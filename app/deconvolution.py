from __future__ import annotations

import argparse
from dataclasses import dataclass, asdict
from functools import lru_cache
import io
import json
from pathlib import Path
import re

import numpy as np
from scipy.ndimage import percentile_filter
from scipy.optimize import nnls
from scipy.signal import find_peaks, peak_widths

VERSION = '1.0.0'
ELECTRON = 0.000548579909
ATOM = {'C': 12.0, 'H': 1.00782503223, 'N': 14.00307400443, 'O': 15.99491461957, 'Cl': 34.968852682, 'Li': 7.0160034366, 'Na': 22.9897692820, 'K': 38.9637064864, 'Rb': 84.9117897379, 'Cs': 132.9054519610}
HIS = 6 * ATOM['C'] + 9 * ATOM['H'] + 3 * ATOM['N'] + 2 * ATOM['O']
ISOTOPES = {
    'C': ((0, 0.0, 0.9893), (1, 1.00335483507, 0.0107)),
    'H': ((0, 0.0, 0.999885), (1, 1.00627674589, 0.000115)),
    'N': ((0, 0.0, 0.99636), (1, 0.99703489445, 0.00364)),
    'O': ((0, 0.0, 0.99757), (1, 1.00421713693, 0.00038), (2, 2.00424499329, 0.00205)),
    'Cl': ((0, 0.0, 0.7576), (2, 1.997049920, 0.2424)),
    'Li': ((0, 0.0, 0.9241), (-1, -1.0008805492, 0.0759)),
    'Na': ((0, 0.0, 1.0),),
    'K': ((0, 0.0, 0.932581), (1, 1.0002916796, 0.000117), (2, 1.9981187715, 0.067302)),
    'Rb': ((0, 0.0, 0.7217), (2, 1.9973907891, 0.2783)),
    'Cs': ((0, 0.0, 1.0),),
}


@dataclass(frozen=True)
class Composition:
    histidine: int
    metal: str
    metals: int
    chloride: int
    hydrogen: int
    charge: int

    @property
    def reference_mz(self):
        return (self.histidine * HIS + self.metals * ATOM[self.metal] + self.chloride * ATOM['Cl'] + self.hydrogen * ATOM['H']) / self.charge - ELECTRON

    @property
    def label(self):
        parts = []
        if self.histidine:
            parts.append(f'{self.histidine}His')
        if self.metals:
            parts.append(f'{self.metals}{self.metal}')
        if self.chloride:
            parts.append(f'{self.chloride}Cl')
        body = '+'.join(parts)
        if self.hydrogen:
            body += f'{self.hydrogen:+d}H'
        return f'[{body}]{self.charge}+'


@lru_cache(maxsize=512)
def element_distribution(element, count):
    distribution = {0: (1.0, 0.0)}
    for unused in range(count):
        updated = {}
        for nominal, (probability, moment) in distribution.items():
            for isotope_nominal, delta, abundance in ISOTOPES[element]:
                key = nominal + isotope_nominal
                old_probability, old_moment = updated.get(key, (0.0, 0.0))
                updated[key] = (old_probability + probability * abundance, old_moment + abundance * (moment + probability * delta))
        distribution = {key: value for key, value in updated.items() if value[0] > 1e-12}
    return distribution


@lru_cache(maxsize=18000)
def envelope(composition):
    counts = {'C': 6 * composition.histidine, 'H': 9 * composition.histidine + composition.hydrogen, 'N': 3 * composition.histidine, 'O': 2 * composition.histidine, 'Cl': composition.chloride, composition.metal: composition.metals}
    distribution = {0: (1.0, 0.0)}
    for element, count in counts.items():
        if not count:
            continue
        updated = {}
        for first_nominal, (first_probability, first_moment) in distribution.items():
            for second_nominal, (second_probability, second_moment) in element_distribution(element, count).items():
                key = first_nominal + second_nominal
                old_probability, old_moment = updated.get(key, (0.0, 0.0))
                probability = first_probability * second_probability
                updated[key] = (old_probability + probability, old_moment + first_moment * second_probability + second_moment * first_probability)
        distribution = {key: value for key, value in updated.items() if value[0] > 1e-10}
    positions, probabilities = [], []
    for nominal in sorted(distribution):
        probability, moment = distribution[nominal]
        if probability >= 0.0005:
            positions.append(composition.reference_mz + moment / probability / composition.charge)
            probabilities.append(probability)
    return np.array(positions), np.array(probabilities)


def read_spectrum(text):
    rows = []
    skipped = 0
    for line_number, line in enumerate(text.lstrip('\ufeff').splitlines(), 1):
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            continue
        fields = re.split(r'[,;\s]+', stripped)
        if len(fields) != 2:
            raise ValueError(f'第 {line_number} 行不是两列 m/z、强度。')
        try:
            rows.append([float(value) for value in fields])
        except ValueError as error:
            if not rows and skipped == 0:
                skipped += 1
                continue
            raise ValueError(f'第 {line_number} 行含有非数值内容。') from error
    if len(rows) < 10:
        raise ValueError('至少需要 10 个数据点；本程序读取连续轮廓谱。')
    data = np.array(rows, dtype=float)
    if not np.isfinite(data).all() or np.any(data[:, 0] <= 0):
        raise ValueError('数据中有非有限数值或非正 m/z。')
    if np.any(np.diff(data[:, 0]) <= 0):
        raise ValueError('m/z 必须严格递增，不能有重复坐标。')
    if np.max(data[:, 1]) <= 0:
        raise ValueError('谱图中没有正强度信号。')
    return data[:, 0], data[:, 1]


def centroid_profile(positions, intensity, minimum_percent):
    baseline = percentile_filter(np.maximum(intensity, 0), percentile=15, size=min(301, len(intensity) // 2 * 2 - 1), mode='nearest')
    corrected = np.maximum(intensity - baseline, 0)
    difference = np.diff(intensity)
    noise = float(np.median(np.abs(difference - np.median(difference))) / 0.9538726)
    threshold = max(float(corrected.max()) * minimum_percent / 100, noise * 5)
    indices, properties = find_peaks(corrected, height=threshold, prominence=threshold)
    if not len(indices):
        raise ValueError('当前阈值下未检出峰；请降低最小峰高，或检查数据。')
    widths, heights, left, right = peak_widths(corrected, indices, rel_height=0.5)
    width_da = np.interp(right, np.arange(len(positions)), positions) - np.interp(left, np.arange(len(positions)), positions)
    return indices, corrected[indices], width_da, baseline, corrected, noise, threshold


def library(metal, maximum_charge, maximum_histidine, include_chloride):
    for charge in range(1, maximum_charge + 1):
        for histidine in range(maximum_histidine + 1):
            for chloride in range(7 if include_chloride else 1):
                for metals in range(min(16, histidine + charge + chloride) + 1):
                    hydrogen = charge - metals + chloride
                    if hydrogen < -histidine or hydrogen > charge:
                        continue
                    if histidine == 0 and (hydrogen != 0 or metals == 0):
                        continue
                    yield Composition(histidine, metal, metals, chloride, hydrogen, charge)


def nearest_index(positions, target):
    insertion = int(np.searchsorted(positions, target))
    possible = [index for index in (insertion - 1, insertion) if 0 <= index < len(positions)]
    return min(possible, key=lambda index: abs(positions[index] - target))


def compress_profile(positions, intensity, budget=7000):
    if len(positions) <= budget:
        return np.column_stack((positions, intensity)).round(6).tolist()
    chosen = {0, len(positions) - 1}
    for group in np.array_split(np.arange(len(positions)), budget // 2):
        chosen.add(int(group[np.argmax(intensity[group])]))
        chosen.add(int(group[np.argmin(intensity[group])]))
    indices = sorted(chosen)
    return np.column_stack((positions[indices], intensity[indices])).round(6).tolist()


def analyze(text, metal='Na', minimum_percent=0.8, tolerance_ppm=100, maximum_charge=2, maximum_histidine=12, include_chloride=True, mass_offset_ppm=0.0, filename='spectrum.txt'):
    if metal not in ('Li', 'Na', 'K', 'Rb', 'Cs'):
        raise ValueError('请选择 Li、Na、K、Rb 或 Cs。')
    if not 0.05 <= minimum_percent <= 20 or not 5 <= tolerance_ppm <= 300:
        raise ValueError('峰高阈值须为 0.05–20%，质量窗口须为 5–300 ppm。')
    if maximum_charge not in (1, 2, 3) or not 1 <= maximum_histidine <= 16 or abs(mass_offset_ppm) > 500:
        raise ValueError('电荷须为 1–3，His 数须为 1–16，质量偏移须在 ±500 ppm 内。')
    positions, intensity = read_spectrum(text)
    indices, peak_height, width_da, baseline, corrected, noise, threshold = centroid_profile(positions, intensity, minimum_percent)
    observed = positions[indices]
    candidates = []
    missing_positions = {}
    for composition in library(metal, maximum_charge, maximum_histidine, include_chloride):
        if not positions[0] - 12 < composition.reference_mz < positions[-1] + 12:
            continue
        theoretical, probabilities = envelope(composition)
        mode = int(np.argmax(probabilities))
        adjusted = theoretical * (1 + mass_offset_ppm * 1e-6)
        if not positions[0] <= adjusted[mode] <= positions[-1]:
            continue
        peak_index = nearest_index(observed, adjusted[mode])
        shift = observed[peak_index] - adjusted[mode]
        if abs(shift) > max(0.015, adjusted[mode] * tolerance_ppm * 1e-6):
            continue
        if peak_height[peak_index] < corrected.max() * minimum_percent / 100:
            continue
        rows = {}
        actual_vector = []
        model_vector = []
        matched = 0
        visible_probability = 0.0
        alignment_window = min(0.05, max(0.012, float(width_da[peak_index]) * 0.75))
        for isotope_position, probability in zip(adjusted + shift, probabilities):
            if not positions[0] <= isotope_position <= positions[-1]:
                continue
            visible_probability += probability
            isotope_index = nearest_index(observed, isotope_position)
            if abs(observed[isotope_index] - isotope_position) <= alignment_window:
                row = int(isotope_index)
                actual = peak_height[row]
                if probability > probabilities.max() * 0.03:
                    matched += 1
            else:
                key = round(float(isotope_position) / 0.02)
                row = missing_positions.setdefault(key, len(observed) + len(missing_positions))
                actual = 0.0
            rows[row] = rows.get(row, 0) + float(probability)
            actual_vector.append(actual)
            model_vector.append(probability)
        if visible_probability < 0.75 or matched < 2:
            continue
        cosine = float(np.dot(actual_vector, model_vector) / max(np.linalg.norm(actual_vector) * np.linalg.norm(model_vector), 1e-20))
        if cosine < 0.65:
            continue
        candidates.append({'composition': composition, 'rows': rows, 'mode_index': peak_index, 'mode_theoretical': float(theoretical[mode]), 'offset_da': float(shift), 'cosine': cosine, 'matched_isotopes': matched})
    full_height = np.zeros(len(observed) + len(missing_positions))
    scale = float(peak_height.max())
    full_height[:len(observed)] = peak_height / scale
    predicted = np.zeros(len(observed))
    components = []
    relative_residual = 1.0
    if candidates:
        used_rows = sorted(set(range(len(observed))).union(*(candidate['rows'] for candidate in candidates)))
        row_lookup = {row: index for index, row in enumerate(used_rows)}
        design = np.zeros((len(used_rows), len(candidates)), dtype=float)
        for column, candidate in enumerate(candidates):
            for row, probability in candidate['rows'].items():
                design[row_lookup[row], column] = probability
        response = full_height[used_rows]
        amplitudes, residual = nnls(design, response, maxiter=max(3000, len(candidates) * 10))
        fitted = design @ amplitudes
        predicted = fitted[:len(observed)] * scale
        relative_residual = float(residual / np.linalg.norm(response))
        total_amplitude = float(amplitudes.sum())
        for candidate, amplitude in zip(candidates, amplitudes):
            if amplitude <= max(total_amplitude * 0.001, 1e-6):
                continue
            composition = candidate['composition']
            mode_observed = float(observed[candidate['mode_index']])
            components.append({**asdict(composition), 'label': composition.label, 'reference_mz': composition.reference_mz, 'mode_theoretical_mz': candidate['mode_theoretical'], 'mode_observed_mz': mode_observed, 'raw_error_ppm': (mode_observed / candidate['mode_theoretical'] - 1) * 1e6, 'local_alignment_da': candidate['offset_da'], 'envelope_cosine': candidate['cosine'], 'matched_isotopes': candidate['matched_isotopes'], 'fitted_height_sum': float(amplitude * scale), 'fitted_share_percent': float(amplitude / total_amplitude * 100), 'neutral_histidine_mass': composition.histidine * HIS, 'ion_mass': composition.reference_mz * composition.charge, 'status': '模型条件下的候选，非结构确证'})
        components.sort(key=lambda component: component['fitted_height_sum'], reverse=True)
    peaks = [{'mz': float(position), 'raw_height': float(intensity[indices[index]]), 'corrected_height': float(peak_height[index]), 'relative_height': float(peak_height[index] / scale * 100), 'fwhm_da': float(width_da[index]), 'fitted_height': float(predicted[index]), 'residual_height': float(peak_height[index] - predicted[index])} for index, position in enumerate(observed)]
    warnings = ['峰归属依赖候选库；同分异构体、配位结构及溶液结合常数不能由本结果确定。', '拟合份额是同位素包络峰高和的模型比例，不是浓度、摩尔分数或总离子流面积。', '按局部峰位对齐模板；原始质量误差仍保留。宽搜索窗口不能替代外部标准校准。', '同一 m/z 可含单电荷与成比例组成的多电荷簇；需用同位素间距及残差审查。']
    if positions[-1] < 900:
        warnings.insert(0, f'数据仅覆盖 m/z {positions[0]:.4f}–{positions[-1]:.4f}；不能推断范围之外的峰，也不能与完整谱直接比较全谱份额。')
    if include_chloride:
        warnings.append('含 Cl 候选是待验证假设；请核对金属盐、组氨酸盐形式、试剂空白与同位素包络。')
    if len(components) == 0:
        warnings.insert(0, '没有满足质量与至少两条同位素线条件的组分；不强行归属。')
    matched_fraction = float(np.minimum(peak_height, predicted).sum() / peak_height.sum())
    if matched_fraction < 0.2:
        warnings.insert(0, '模型仅解释少量已检出信号：当前组分表属于低覆盖探索结果，不宜用于化学结论或比较组成比例。')
        for component in components:
            component['status'] = '低覆盖探索候选，不可作确认归属'
    harmonic_candidates = []
    for child_index, child in enumerate(observed):
        if child > 350:
            continue
        for order in (2, 3):
            if child * order > positions[-1]:
                continue
            parent_index = nearest_index(observed, child * order)
            error_ppm = (child * order / observed[parent_index] - 1) * 1e6
            if abs(error_ppm) < 100 and peak_height[parent_index] > 5 * peak_height[child_index]:
                harmonic_candidates.append({'child_mz': float(child), 'parent_mz': float(observed[parent_index]), 'frequency_multiple': order, 'error_ppm': float(error_ppm), 'status': '可能为倍频伪峰或重叠电荷态，未自动删除'})
    if harmonic_candidates:
        warnings.append('发现低质量峰与强峰成 1/2 或 1/3 关系；可能为倍频伪峰或重叠电荷态，详见 JSON 的 harmonic_candidates。未自动删除。')
    return {'version': VERSION, 'filename': filename, 'settings': {'metal': metal, 'minimum_percent': minimum_percent, 'tolerance_ppm': tolerance_ppm, 'maximum_charge': maximum_charge, 'maximum_histidine': maximum_histidine, 'include_chloride': include_chloride, 'mass_offset_ppm': mass_offset_ppm}, 'summary': {'point_count': len(positions), 'mz_min': float(positions[0]), 'mz_max': float(positions[-1]), 'base_peak_mz': float(positions[np.argmax(intensity)]), 'base_peak_height': float(intensity.max()), 'noise_estimate': noise, 'detection_threshold': threshold, 'detected_peaks': len(indices), 'candidate_count': len(candidates), 'reported_components': len(components), 'matched_peak_height_fraction': matched_fraction, 'relative_l2_residual_including_missing_isotopes': relative_residual}, 'warnings': warnings, 'components': components, 'peaks': peaks, 'harmonic_candidates': harmonic_candidates, 'profile': compress_profile(positions, intensity), 'baseline': compress_profile(positions, baseline)}


def components_csv(result):
    import csv
    output = io.StringIO(newline='')
    columns = ['label', 'histidine', 'metal', 'metals', 'chloride', 'hydrogen', 'charge', 'reference_mz', 'mode_theoretical_mz', 'mode_observed_mz', 'raw_error_ppm', 'local_alignment_da', 'envelope_cosine', 'matched_isotopes', 'fitted_height_sum', 'fitted_share_percent', 'neutral_histidine_mass', 'ion_mass', 'status']
    writer = csv.DictWriter(output, fieldnames=columns)
    writer.writeheader()
    writer.writerows(result['components'])
    return output.getvalue()


def main():
    parser = argparse.ArgumentParser(description='His / alkali-metal positive ESI isotope-envelope NNLS deconvolution')
    parser.add_argument('input', type=Path)
    parser.add_argument('--metal', choices=['Li', 'Na', 'K', 'Rb', 'Cs'], default='Na')
    parser.add_argument('--output', type=Path, default=Path('result'))
    parser.add_argument('--no-chloride', action='store_true')
    parser.add_argument('--max-charge', type=int, choices=[1, 2, 3], default=2)
    parser.add_argument('--min-percent', type=float, default=0.8)
    parser.add_argument('--tolerance-ppm', type=float, default=100)
    parser.add_argument('--mass-offset-ppm', type=float, default=0)
    args = parser.parse_args()
    result = analyze(args.input.read_text(encoding='utf-8-sig'), metal=args.metal, minimum_percent=args.min_percent, tolerance_ppm=args.tolerance_ppm, maximum_charge=args.max_charge, include_chloride=not args.no_chloride, mass_offset_ppm=args.mass_offset_ppm, filename=args.input.name)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'result.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    # newline='' 防止 Windows 文本模式把 components_csv 已有的 CRLF 再翻译成 CRCRLF。
    with (args.output / 'components.csv').open('w', encoding='utf-8-sig', newline='') as handle:
        handle.write(components_csv(result))
    print(json.dumps(result['summary'], ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
