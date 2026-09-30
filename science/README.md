# TMC1 科研公开协作资料

这里保留第二轮（R2）科研证据、三物种受体修订、第三轮（R3）P0库与筛选证据，以及R3阶段零USH1B公开资料。R2和R3分目录保存，便于三人协作引用同一份数据。

## 从哪里开始

- [R2核心科研报告](R2/evidence/06_report/第二轮核心总报告.md)：方法、结果与讨论；同目录附可编辑Word及PDF。
- [R2一页摘要](R2/evidence/06_report/一页摘要_阅读说明.md)：快速阅读入口。
- [三物种修订](R2/three_species_revision/08_三物种构象修订/第二轮三物种修订说明.md)：人、小鼠、线虫结构的来源与适用边界。
- [R3 P0库证据](R3/P0/01_library_evidence/stage2_library_evidence.md)：当前库的来源和分母。
- [R3匹配口径结果](R3/P0/02_screening_validation/matched_MW350_metrics.md)：与R2口径分开阅读。
- [R3阶段零公开资料](R3/stage0/stage0_usher1b_public_evidence.md)：USH1B / MYO7A研究主线及科学来源。
- [复现入口与缺口](REPRODUCIBILITY.md)：现有脚本与待补输入。

对接分数与几何/接触指标属于计算排序和方法控制证据，尚不能证明生化活性、特异性或临床疗效。R3疾病及临床资料是2026-09-30的公开来源记录，应按需要更新。

## 现有代码（14个）

代码来自R2匿名科学包。完整重算仍缺历史R1输入、部分R2库构建/排序/对接/生成方法及外部工具环境；因此本仓库交付已有证据和代码，不声称完整复现已完成。以下按用途提供入口，未运行科研计算或实现测试。

| 脚本 | 用途 |
|---|---|
| [assemble_incremental_pharmacology.py](R2/evidence/05_independent_audit/formal319_pharmacology/assemble_incremental_pharmacology.py) | 整理增量药理来源与结果 |
| [fetch_incremental_chembl.py](R2/evidence/05_independent_audit/formal319_pharmacology/fetch_incremental_chembl.py) | 获取公开ChEMBL药理来源 |
| [build_verified_stage_figures.py](R2/evidence/06_figures/build_verified_stage_figures.py) | 生成阶段结果图 |
| [plot_forward_three_arms.py](R2/evidence/06_figures/plot_forward_three_arms.py) | 绘制前瞻三臂比较图 |
| [plot_from_verified_csv.py](R2/evidence/06_figures/plot_from_verified_csv.py) | 从结果CSV绘图 |
| [register_reviewed_figures.py](R2/evidence/06_figures/register_reviewed_figures.py) | 登记图件及来源 |
| [compose_pose_shift_figure.py](R2/evidence/09_structure/compose_pose_shift_figure.py) | 组合位姿迁移图 |
| [compose_structure_figure.py](R2/evidence/09_structure/compose_structure_figure.py) | 组合结构图 |
| [launch_render.py](R2/evidence/09_structure/launch_render.py) | 调用外部结构渲染工具 |
| [render_pose_shift.py](R2/evidence/09_structure/render_pose_shift.py) | 渲染位姿迁移图 |
| [render_structure_panels.py](R2/evidence/09_structure/render_structure_panels.py) | 渲染结构面板 |
| [verify_delivery.py](R2/evidence/09_structure/verify_delivery.py) | 原交付内容检查入口 |
| [prepare_portable.py](R2/evidence/replay/prepare_portable.py) | 便携重放准备入口，需要完整R1/R2输入 |
| [redraw_f3_f4.py](R2/evidence/replay/redraw_f3_f4.py) | 依据已有CSV重绘F3/F4，需要相应设计输入表 |

## 公开副本说明

原始材料未被修改。选定公开副本清除了作者属性、本机路径或报告中的内部联络说明；科学数值与结果保留。公开文件的完整清单及SHA256见[内容清单](PUBLIC_CONTENT_MANIFEST.json)和[校验列表](SHA256SUMS.txt)。部分原嵌入清单描述原始文件，公开副本以这里的新清单为准。筛选范围、排除数量及检查限制见[公开内容检查报告](PUBLIC_CONTENT_INSPECTION.json)。
