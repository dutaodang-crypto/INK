"""一键复核：用当前程序重跑六份附带样品，并与参考 summary.json 对照。

用法（在本文件夹内）：
    python verify.py
    python verify.py --summary "..\\data\\summary.json"

默认自动查找 ../data/summary.json（仓库布局）或 ../analysis_data/summary.json（打包布局）。

判定：基峰 m/z ±0.01 Th、报告组分数完全相等、已匹配峰高比例 ±0.002。
退出码：0 全部通过；1 存在不一致或参考中缺少样品记录；2 找不到或无法读取参考文件。
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from deconvolution import analyze  # noqa: E402

SAMPLES = [
    ("HisNa1.5mM", "Na"),
    ("HisNa5mM", "Na"),
    ("HisLi5mM", "Li"),
    ("HisK5mM", "K"),
    ("HisRb5mM", "Rb"),
    ("HisCs5mM", "Cs"),
]

# 与参考 summary.json 生成时一致的设置
SETTINGS = dict(
    minimum_percent=0.8,
    tolerance_ppm=100,
    maximum_charge=2,
    maximum_histidine=12,
    include_chloride=True,
    mass_offset_ppm=0.0,
)

TOL_MZ = 0.01
TOL_MATCH = 0.002

REFERENCE_CANDIDATES = (
    HERE.parent / "data" / "summary.json",           # 仓库布局
    HERE.parent / "analysis_data" / "summary.json",  # 打包布局
)


def default_reference() -> Path:
    for candidate in REFERENCE_CANDIDATES:
        if candidate.is_file():
            return candidate
    return REFERENCE_CANDIDATES[0]


def load_reference(path: Path) -> tuple[dict | None, str | None]:
    if not path.is_file():
        return None, f"找不到参考文件 {path}"
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        return None, f"无法读取参考文件 {path}（{error}）"
    if not isinstance(data, dict):
        return None, f"参考文件 {path} 不是 JSON 对象。"
    return data, None


def reference_values(row) -> tuple[float, int, float] | None:
    """取出一条参考记录；结构不完整时返回 None。"""
    if not isinstance(row, dict):
        return None
    try:
        return (
            float(row["base_peak_mz"]),
            int(row["reported_components"]),
            float(row["matched_peak_height_fraction"]),
        )
    except (KeyError, TypeError, ValueError):
        return None


def main() -> int:
    parser = argparse.ArgumentParser(description="复核六份样品的解卷积结果")
    parser.add_argument("--summary", default=str(default_reference()), help="参考 summary.json 路径")
    args = parser.parse_args()

    if hasattr(sys.stdout, "reconfigure"):
        try:
            sys.stdout.reconfigure(encoding="utf-8")
        except Exception:
            pass

    ref, error = load_reference(Path(args.summary))
    if error:
        print(f"错误：{error}")
        print("未做判定。请用 --summary 指定参考 summary.json。")
        return 2

    header = f"{'样品':<12}{'基峰 m/z':>10}{'参考':>10}{'组分':>6}{'参考':>6}{'匹配比例':>9}{'参考':>8}  判定"
    print(header)
    print("-" * 78)

    failures = 0
    for key, metal in SAMPLES:
        sample = HERE / "samples" / f"{key}-cluster-calibration-1.txt"
        if not sample.is_file():
            print(f"{key:<12}  缺少文件 {sample}")
            failures += 1
            continue

        result = analyze(sample.read_text(encoding="utf-8-sig"), metal=metal, filename=key, **SETTINGS)
        got = result["summary"]
        values = reference_values(ref.get(key))

        if values is None:
            print(f"{key:<12}{got['base_peak_mz']:>10.4f}{'-':>10}{got['reported_components']:>6}{'-':>6}"
                  f"{got['matched_peak_height_fraction']:>9.3f}{'-':>8}  （参考中缺少该样品）")
            failures += 1
            continue

        row_base, row_components, row_match = values
        same = (
            abs(got["base_peak_mz"] - row_base) <= TOL_MZ
            and got["reported_components"] == row_components
            and abs(got["matched_peak_height_fraction"] - row_match) <= TOL_MATCH
        )
        failures += 0 if same else 1
        print(f"{key:<12}{got['base_peak_mz']:>10.4f}{row_base:>10.4f}"
              f"{got['reported_components']:>6}{row_components:>6}"
              f"{got['matched_peak_height_fraction']:>9.3f}{row_match:>8.3f}"
              f"  {'一致' if same else '不一致'}")

    print()
    if failures:
        print(f"复核失败：{failures} 份样品与参考结果不一致或缺少参考记录。")
        return 1
    print("复核通过：六份样品的基峰、报告组分数与匹配比例均与原分析一致。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
