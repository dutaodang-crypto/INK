# His–碱金属 ESI 质谱：解卷积、峰归属与浓度比较

组氨酸（L-His）与碱金属氯化物（LiCl / NaCl / KCl / RbCl / CsCl）混合溶液的 ESI 质谱数据分析。
本仓库包含**一个可运行的解卷积小程序**、**六份原始文本谱图**、**图文报告**与**全部结果数据**。

## 三个任务

| 任务 | 交付内容 | 位置 |
| --- | --- | --- |
| Task 1 开发解卷积小程序 | 网页界面 + 命令行 + CSV/JSON/SVG 导出 | `app/` |
| Task 2 峰是什么峰；与金属离子浓度是否有关 | 四类信号归属表、两份钠谱图对比、结论与边界 | `report/His_碱金属_ESI分析报告.html` |
| Task 3 不同碱金属如何影响峰值 | 五个金属的基峰位移、同位素包络、主导峰族对比 | 同上 + `data/major_peak_annotations.csv` |

## 结论速览

**Task 2：峰的构成。** 谱图里的峰不是一个化合物一个峰，而是四类信号叠加：

1. **质子化组氨酸聚集体**：[4His+H]⁺ (621.28518)、[8His+2H]²⁺；
2. **金属加合 / 金属取代簇**：[3His+Na]⁺、[3His+2Na−H]⁺，每换一个金属净增 m(M)−m(H) = 21.98194 Th（Na）；
3. **含氯混合簇**：[3His+5Na+Cl−3H]⁺ (612.10211) = 554 系列 + 一个中性 NaCl（+57.95862 Da）；
4. **纯金属氯化物簇**：[4Cs+3Cl]⁺ (636.4958)，不含组氨酸。

**Task 2：与浓度有关，但不能换算成浓度。** NaCl 由 1.5 → 5 mmol/L，基峰从 621.281 变为 554.136；
I(554)/I(621) 由 0.357 升到 5.77（约 16.2 倍），I(554)/I(510) 由 0.676 升到 3.34（约 4.9 倍）。
变化的是**相对分布**：高 [Na] 时富钠/含氯簇的相对贡献上升。由于 NaCl 同时改变 Na⁺、Cl⁻ 与离子强度，
且无重复实验与内标，**不能说成"Na⁺ 的结合常数"**，也不能把谱峰高直接当浓度。

**Task 3：碱金属换了，三类特征同时变。**
理论 m/z 平移（3His+M 系列：Li 472.2239 → Na 488.1977 → K 504.1716 → Rb 550.1197 → Cs 598.1133）；
天然同位素包络改变（Li 有 ⁶Li −1.0009/z、K 有 ⁴¹K +2/z、Rb 85/87 = 72.17/27.83、Cl 35/37 = 75.76/24.24）；
主导峰族改变（Li、Na 为富金属 His 簇，K、Rb 的 621 峰最强，Cs 以含氯混合簇与纯盐簇为主）。
现有数据（单浓度、无重复、无标定）**不足以给出可靠的结合能力排序**。

## 目录结构

```
app/                        解卷积程序
  deconvolution.py            算法本体（候选模型 + 同位素包络 + NNLS）
  server.py                   本地网页服务（只监听 127.0.0.1）
  index.html                  网页界面
  verify.py                   一键复核（6 个样品与报告数字逐项比对）
  start.bat                   启动脚本（Windows）
  requirements.txt            numpy>=2.0,<3 / scipy>=1.14,<2
  README_使用说明.md           参数、模型边界与验证情况
  samples/                    6 份原始文本谱图（第一列 m/z，第二列强度）
report/
  His_碱金属_ESI分析报告.html   图文报告（离线可看，图片内嵌）
  HisMS_界面预览.png           界面截图
  figures/                    报告中的四张图（PNG + SVG）
data/
  summary.json                六份样品的解卷积汇总
  His*_components.csv         每个样品报告出的组分表
  His*_result.json            完整结果（含曲线、基线、残差、候选峰）
  major_peak_annotations.csv  主要峰的候选归属
  sodium_PDF_peak_heights.csv 两份钠谱图 PDF 的读峰表（相对峰高）
  algorithm_validation.json   合成重叠峰恢复测试
  interface_validation.json   界面/接口自检记录
  provenance.json             软件版本、参数、原始文件 SHA-256
  先读我.txt                   数据来源与实验条件的最终确认
```

## 如何运行

需要 Python 3.10+ 与 numpy、scipy：

```bash
pip install -r app/requirements.txt
python app/server.py          # 浏览器打开 http://127.0.0.1:8765/
```

或用命令行直接分析单个文件：

```bash
python app/deconvolution.py app/samples/HisNa5mM-cluster-calibration-1.txt --metal Na
python app/verify.py --summary data/summary.json        # 复核六份样品
```

另外提供**免安装一键包**（内置 Python + numpy + scipy，解压双击即用，约 71 MB），
因为超过仓库常规体积，单独放在 Release / 网盘中，未纳入本仓库。

## 数据说明

- 所有样品 His 均为 5 mmol/L；金属盐均为氯化物；采集参数一致。
- 金属浓度：NaCl 两组 1.5 与 5 mmol/L，LiCl / KCl / RbCl / CsCl 各 5 mmol/L。
- `HisNa1.5mM` 的 TXT 只覆盖 m/z 200.0001–262.5101（原始文件本身不完整），
  报告中两组钠的高质量区相对峰高取自两份钠谱图 PDF 的矢量谱线读数，已在正文注明。
- 采集日期 2022-01-30 与 2022-02-07，正离子模式，Scans 1。

## 方法边界

这是"候选组成模型 + 同位素包络 + 非负最小二乘"的解卷积，不是数据库匹配，也不是无先验通用反演；
拟合份额是包络峰高的模型比例，不等于摩尔分数或离子流积分。同一 m/z 可同时对应单电荷与双电荷组成
（如 [3His+4Na−3H]⁺ 与 [6His+8Na−6H]²⁺ 理论 m/z 相同），程序允许两者并存，分量依赖模型。
低质量区若出现与强峰成 1/2、1/3 关系的峰，可能是倍频伪峰，程序仅在 JSON 中提示，不作定性结论。
