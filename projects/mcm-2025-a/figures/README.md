# 2501909 踏面磨损图：热图 + 横向/纵向投影

本目录只包含**科研绘图**交付物，不修改 WDM 主程序与主测试模块。首轮使用**明确标记的合成磨损格**，用于把坐标、单位、投影定义和渲染流程固定下来；拿到真实矩阵后，同一命令可直接重绘并做论文图对比。

> **合成声明**：默认 `--demo` 数据为闭式解析场景，**不是论文实测数据**。图内、文件名、manifest 三处均带 `SYNTHETIC` 标记，禁止当作论文结果引用。

正式渲染产物。 [PNG](artifacts/synthetic/fig_2501909_tread_wear_synthetic.png) · [PDF](artifacts/synthetic/fig_2501909_tread_wear_synthetic.pdf) · [SVG](artifacts/synthetic/fig_2501909_tread_wear_synthetic.svg) · [生成记录](artifacts/synthetic/figure_manifest.json)

## 1. 坐标与单位约定

与 MCM-ICM-Lab PR #2 的 WDM 输出一致：

| 项 | 约定 |
|---|---|
| `depth_mm[y][x]` | 磨损深度，正值为材料去除；形状 `(ny, nx)`，行主序 |
| `x` | 台阶**左右方向**；`x = 0` 为左端，向右递增 |
| `y` | 踏面**前后方向**；`y = 0` 为**前缘（nosing）**，向踏面内部递增 |
| 原点 `O(0,0)` | 前缘—左端角点；图中用空心圆和文字标出 |
| `dx_mm`, `dy_mm` | 单元格物理尺寸，mm；单元面积 `dx·dy`，mm² |
| 深度单位 | mm |
| 投影单位 | mm³（每条网格列/行上的磨损体积） |

两个边际投影的严格定义（与程序输出同名）：

```text
x_wear_volume_mm3[i] = dx·dy · Σ_j depth_mm[j, i]      # 沿 y 求和，横向剖面
y_wear_volume_mm3[j] = dx·dy · Σ_i depth_mm[j, i]      # 沿 x 求和，纵向剖面
```

总磨损体积 `V = dx·dy · Σ depth_mm`，因此必有
`Σ x_wear_volume_mm3 = Σ y_wear_volume_mm3 = V`。这是图中和测试里的第一道守恒检查。

## 2. 合成场景（闭式、可手算）

`nx = 24`, `ny = 16`, `dx = dy = 10 mm`，踏面 240 mm × 160 mm，单元面积 100 mm²。

\[
d[y,x] = 3.0 \cdot \frac{g_x}{12} \cdot \frac{h_y}{16}\ \text{mm},
\qquad
g_x = 1,2,\ldots,12,12,\ldots,2,1,
\qquad
h_y = 16-j.
\]

- 横向 `g`：居中行走带的三角权重（最两侧 1，中心 12）。
- 纵向 `h`：自前缘 `h=16` 向内部线性衰减到 `h=1`。
- 峰值 3.000 mm，位于前缘中央两列（x = 115、125 mm 的单元中心）。

**手算核对表**（精确值，测试与 `--self-test` 使用同一组数）：

| 量 | 计算 | 期望值 |
|---|---|---|
| 最左列横向投影 | `100·3·(1/12)·(136/16)` | 212.5 mm³ |
| 中心列横向投影 | `212.5 × 12` | 2550 mm³ |
| 前缘行纵向投影 | `100·3·(156/12)·(16/16)` | 3900 mm³ |
| 最内行纵向投影 | `3900 × (1/16)` | 243.75 mm³ |
| 总磨损体积 | `100·3·13·8.5` | 33150 mm³ |
| 前缘中心深度 | `3.0·1·1` | 3.000 mm |
| 最内左角深度 | `3.0·(1/12)·(1/16)` | 0.015625 mm |

两投影求和都应回到 33150 mm³；`--self-test` 和 `pytest` 同时断言这些点值与守恒式。

## 3. 输入格式

1. **CSV/TXT 矩阵**（推荐）：逗号或空白分隔，**无表头**；第 `j` 行是 `y` 的第 `j` 个单元（行 0 = 前缘），第 `i` 列是 `x` 的第 `i` 个单元（列 0 = 左端）。允许首行为非数值表头（自动跳过）。
2. **NumPy**：`.npy` 二维数组；`.npz` 中键名 `depth_mm`（缺省取第一个数组）。
3. 可选的程序输出核对文件：`x_wear_volume_mm3.csv`、`y_wear_volume_mm3.csv`，单列数值或 `(坐标, 体积)` 两列；用 `--check-x-proj/--check-y-proj` 比对。

格式示例见 [`data/example_depth_5x4.csv`](data/example_depth_5x4.csv)（4 行 `y` × 5 列 `x`，单位 mm，仅演示格式，不是实验数据）。

## 4. 输出

默认输出到 `--outdir`（默认 `figures/out/2501909`）：

```text
fig_2501909_tread_wear_synthetic.png   300 dpi 位图
fig_2501909_tread_wear_synthetic.pdf   矢量（字体嵌入，fonttype 42）
fig_2501909_tread_wear_synthetic.svg   矢量（默认文字转路径）
x_wear_volume_mm3.csv                  x_mm, x_wear_volume_mm3
y_wear_volume_mm3.csv                  y_mm, y_wear_volume_mm3
figure_manifest.json                   输入哈希、网格、单位、总量、核对结果、软件版本
```

真实矩阵（非 `--demo`）默认文件名是 `fig_2501909_tread_wear.*`。

图版布局：

- **(a)** 踏面深度热图；横轴 `x`（左→右），纵轴 `y`（前缘→内部），原点/前缘在左下角；色标为深度 mm。
- **(b)** 热图正下方：横向投影 `x_wear_volume_mm3`，与 (a) 共用 x 轴。
- **(c)** 热图右侧：纵向投影 `y_wear_volume_mm3`，与 (a) 共用 y 轴（前缘在下）。
- 只在峰值与总量处直接标注；底部一行给出坐标约定、网格尺寸与合成标记。

## 5. 重绘命令

```bash
# 合成演示 + 自检（不写图）
python plot_stair_wear.py --demo --self-test

# 合成演示出图
python plot_stair_wear.py --demo --outdir figures/out/2501909

# 真实 WDM 矩阵出图，并与程序输出交叉核对
python plot_stair_wear.py \
  --input /path/to/depth_mm.csv \
  --dx-mm 10 --dy-mm 10 \
  --check-x-proj /path/to/x_wear_volume_mm3.csv \
  --check-y-proj /path/to/y_wear_volume_mm3.csv \
  --source-label "WDM depth_mm[y][x] — 真实矩阵" \
  --outdir figures/out/2501909_real
```

一条命令版本（本目录内、合成演示）：

```bash
python plot_stair_wear.py --demo --outdir figures/out/2501909
```

## 6. 测试与验证

```bash
python -m pytest -q test_plot_stair_wear.py
```

覆盖：闭式场点值、两投影手算、均匀场手算、随机场守恒、非法输入（1 行/列、NaN/Inf、负深度）、CSV/表头读取、程序输出交叉核对、`--self-test`、完整渲染 smoke test。

**当前验证状态**

- [x] 合成场的投影与总量已由闭式表达式手算（上表）。
- [x] `python plot_stair_wear.py --demo --self-test` 通过。
- [x] `python -m pytest -q test_plot_stair_wear.py` 共 17 项通过。
- [x] PNG、PDF、SVG 均成功生成，已目视检查 PDF。
- [x] 同一合成网格输入 [WDM 实现](https://github.com/JunQianNuo/MCM-ICM-Lab/blob/7b33a076be7e55c743b7bdce4a6d81eef1f54ec9/projects/mcm-2025-a/wdm.py)，再以 `--check-x-proj` 和 `--check-y-proj` 核对。两轴均逐格一致。
- [ ] 真实矩阵（WDM 输出）重绘并与原论文图对比。

浏览器预览 `preview/tread_wear_preview.html` 可实时重算同一合成格。正式图使用 Matplotlib 的 cividis 色图。

## 7. 环境与字体

- 依赖：Python ≥ 3.10，`numpy`、`matplotlib`（Agg 后端）。
- 中文字体：自动查找 `Noto Sans CJK SC` / `Source Han Sans SC` / `WenQuanYi Zen Hei` / `SimHei` / `Microsoft YaHei`；找不到时标签自动退回英文，避免方框。
- SVG 默认 `svg.fonttype = path`（跨系统显示稳定）；需要改字时加 `--editable-svg`，并在打开机器上安装同一 CJK 字体。
- PDF 使用 `fonttype 42` 嵌入 TrueType；热图 `rasterized=True`，文字保持矢量。

## 8. 落地到仓库

目标位置：**`projects/mcm-2025-a/figures/`**（本包不触碰 WDM 主代码）。

```text
projects/mcm-2025-a/figures/
├── plot_stair_wear.py
├── test_plot_stair_wear.py
├── README.md
├── data/example_depth_5x4.csv
└── preview/tread_wear_preview.html
```

`figures/out/` 用于本地重绘。首轮选定的图、两轴投影 CSV 与 `figure_manifest.json` 保存在 `artifacts/synthetic/`。

## 9. 局限与下一步

- 合成场是**演示/校验**数据，不能用于论文结论；真实矩阵到达后必须替换并保留 `figure_manifest.json` 的输入 SHA-256。
- 当前投影是**未平滑的边际求和**（对应 WDM 输出）；若论文图使用了平滑、插值或面积归一化，需要按论文口径另加归一化选项，并在图中注明。
- 未包含测量不确定性、网格配准误差或不同踏面的对齐；这些属于后续验证设计，不在本绘图模块内。
- 与论文图对比时，先核对坐标原点、轴向、深度单位、网格尺寸和投影量。列体积 `V_x` 除以整列面积 `dx·dy·ny`，得到该列平均深度。
