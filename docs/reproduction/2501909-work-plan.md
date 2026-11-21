# 2501909 首轮复现分工

原文为《Stair Wear: Traces of History》（Team 2501909）。已在项目机器上找到原论文 PDF。原文用于定义模型和核对报告值，测量矩阵等缺项在证据记录中跟踪。 [在线归档索引](https://wxj630.github.io/Math-Modeling-BAO/best_practie/bao-mcm-2025-a-stair-wear.html)

| 工作 | 主责 | 首轮交付 | 交叉审查 |
| --- | --- | --- | --- |
| 原文证据与案例输入 | Codex 论文研究员 | `docs/reproduction/2501909-paper-evidence.md` 汇总原文公式、坐标、表图和案例数据缺项 | 实验验收核对来源和可复算性 |
| 扩展问题与可辨识性 | Claude 论文研究员 | `docs/reproduction/2501909-identifiability.md` 记录所需证据、反例和实验判据 | 原文证据负责人核对论文口径 |
| 磨损图表 | 编码与绘图 | `projects/mcm-2025-a/figures/` 提供热图、横纵投影和重绘命令 | 数值实现核对坐标及体积守恒 |
| 独立实验验收 | 验收协调者 | `docs/reproduction/2501909-validation.md` 验证 WVM、核对案例回代，记录 WDM 验收清单 | 数值实现复核问题和复现命令 |
| WDM 数值实现 | 开发协调者 | `projects/mcm-2025-a/wdm.py` 与测试。先做深度栅格的横纵投影，随后核对混合分布拟合口径 | 原文证据负责人核对式 (3.17) 和坐标，实验验收独立运行 |

交付使用独立工作树和分支。合成网格的结果用于检验代码，爱丁堡案例的复算依赖单独核实的测量输入。
