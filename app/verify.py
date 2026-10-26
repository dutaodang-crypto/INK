"""一键复核：用当前程序重跑六份附带样品，并与 analysis_data/summary.json 对照。

用法（在本文件夹内）：
    python verify.py
    python verify.py --summary "..\\analysis_data\\summary.json"

判定：基峰 m/z ±0.01 Th、报告组分数完全相等、已匹配峰高比例 ±0.002。
全部通过时退出码为 0，否则为 1。
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

# 与 analysis_data/summary.json 生成时一致的设置
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


def load_reference(path: Path) -> dict | None:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    default = HERE.parent / "analysis_data" / "summary.json"
    parser = argparse.ArgumentParser(description="复核六份样品的解卷积结果")
    parser.add_argument("--summary", default=str(default), help="参考 summary.json 路径")
    args = parser.parse_args()

    if hasattr(sys.stdout, "reconfigure"):
        try:
            sys.stdout.reconfigure(encoding="utf-8")
        except Exception:
            pass

    ref = load_reference(Path(args.summary))
    if ref is None:
        print(f"警告：找不到参考文件 {args.summary}，只打印本次结果，不做判定。")
        print()

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
        row = ref.get(key) if ref else None

        if row is None:
            print(f"{key:<12}{got['base_peak_mz']:>10.4f}{'-':>10}{got['reported_components']:>6}{'-':>6}"
                  f"{got['matched_peak_height_fraction']:>9.3f}{'-':>8}  （无参考）")
            continue

        same = (
            abs(got["base_peak_mz"] - row["base_peak_mz"]) <= TOL_MZ
            and got["reported_components"] == row["reported_components"]
            and abs(got["matched_peak_height_fraction"] - row["matched_peak_height_fraction"]) <= TOL_MATCH
        )
        failures += 0 if same else 1
        print(f"{key:<12}{got['base_peak_mz']:>10.4f}{row['base_peak_mz']:>10.4f}"
              f"{got['reported_components']:>6}{row['reported_components']:>6}"
              f"{got['matched_peak_height_fraction']:>9.3f}{row['matched_peak_height_fraction']:>8.3f}"
              f"  {'一致' if same else '不一致'}")

    print()
    if ref is None:
        print("未做判定：缺少参考文件。")
        return 0
    if failures:
        print(f"复核失败：{failures} 份样品与参考结果不一致。")
        return 1
    print("复核通过：六份样品的基峰、报告组分数与匹配比例均与原分析一致。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
