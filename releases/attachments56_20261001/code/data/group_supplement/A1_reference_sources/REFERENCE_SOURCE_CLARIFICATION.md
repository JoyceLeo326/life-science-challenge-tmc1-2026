# 参考化合物来源补核

## 本次获得的实际进展

已重新读取论文指定作者仓库的完整文件树与 10 份活性配体 SDF，并通过 ZINC 团队的 CartBlanche 接口补回 3 个候选的精确数据库结构。3 个新结构均执行了原冻结制备流程，均因原规则不支持质子化后非环状叔胺的未指定氮构型而停止，原冻结流程未生成新的对接分数。另立两种 N 几何的敏感性分析已完成 6 次对接，结果见 `n_geometry_sensitivity/RESULTS.md`，不并入原冻结统计。

## 必须分开的三个集合

| 集合 | 当前可核实内容 | 使用边界 |
|---|---|---|
| 论文药效团训练集合 | 论文与作者 README 列出 13 个训练分子；10 个模型的 active SDF 合计 78 条记录，标准 InChIKey 去重后为 12 个分子 | 78 是跨模型重复记录数，不能计作 78 个独立化合物；也不能改写为 68 个独立阻断剂 |
| 论文外植体染料实验中降低染料摄取的 12 个候选 | 补充数据 Fig8_stats 已逐行映射；本次补回其中 3 个 ZINC 的数据库结构，另 2 个尚未解析 | 与上行训练配体是不同集合，两个集合刚好同为 12 不代表身份相同；染料摄取是功能相关代理读数，不等于已证实直接结合人 TMC1 |
| 项目曾提到的 68 条参考库 | 本次核查的论文、补充数据及作者冻结仓库未找到可逐条追溯的 68 个独立分子清单 | 不能制作“68 个已确认阻断剂”的分布结论；目前仍应标为来源未证实 |

作者仓库冻结提交为 `74274dad76201258adbd76c7058b705170c4b857`，完整树共 41 个文件。首个 GitHub API 请求返回 403 后，正常 Git 读取成功，已保存文件树、10 份 SDF 的 SHA256 与全部 78 条记录。该提交与本次检查到的远端 main 相同。[作者仓库](https://github.com/ramirezlab/Drug-design-targeting-TMC1/tree/74274dad76201258adbd76c7058b705170c4b857)

12 个训练配体记录名为 UoS-7692、UoS-3607、UoS-3606、UoS-962、Proto-1、E6-berbamine、hexamethyleneamiloride、amsacrine、phenoxybenzamine、carvedilol derivative 13、ORC-13661、FM1-43。论文另列的 UoS-5247 未出现在这 10 份 active SDF 中；这一文件覆盖缺口不等于论文没有使用它。另列的 benzamil、tubocurarine、DHS 是额外参考分子，不能混入药效团训练计数。[作者训练集合说明](https://github.com/ramirezlab/Drug-design-targeting-TMC1/blob/74274dad76201258adbd76c7058b705170c4b857/3D_Pharmacophores_Models/README.md)

## 五个 ZINC 候选的精确结构状态

| 论文编号 | 官方返回编号 | 来源核对 | 原冻结制备流程 |
|---|---|---|---|
| ZINC24739924 | ZINC000024739924 | 已取得 SMILES；重算 InChIKey 与官方值一致 | 非环状叔胺质子化后氮构型不受旧规则支持，停止 |
| ZINC58438263 | ZINC000058438263 | 已取得 SMILES；重算 InChIKey 与官方值一致 | 同上，停止 |
| ZINC12986242 | ZINC000012986242 | 已取得 SMILES；重算 InChIKey 与官方值一致 | 同上，停止 |
| ZINC06530230 | 未返回记录 | ZINC15 验证页；CartBlanche 批量与单独查询均无记录；PubChem 别名检索无记录 | 未运行 |
| ZINC12430014 | 未返回记录 | 同上 | 未运行 |

使用官方 CartBlanche 当前代码支持的 `POST /substances.json`，查询原 ZINC 数字编号，并读取返回任务结果。旧 wiki 中的 GET 示例本次只返回网页，未误记为结构数据。[官方接口源代码](https://github.com/docking-org/cartblanche22/blob/94dafee595e0f02b66b1a19cdf6d73ebf2693a11/backend/cartblanche/main/search.py) [官方检索说明](https://wiki.docking.org/index.php/Zinc22:Searching)

三份成功匹配均按请求编号逐条取交集，并核对结构重算 InChIKey。接口对两份未返回编号仍给出 `missing=[]`，所以未把这个字段当作完整性证明。原始响应里的供应商价格字段未纳入结论。`resolved_ZINC_primary_identity_2D.sdf` 是依据精确数据库 SMILES 生成的二维结构文件，可用于身份核对；对接输入与 pH 处理另有记录。

## Ceforanide 制备问题的具体原因

原始 ChEMBL 母体的两个碳手性均有指定，没有缺失碳构型。旧流程的 pH 7.4 规则给原子 15 加质子，产生未指定氮构型；该 N 与羰基 C 相连，且属于四元与六元两个稠合环。现有规则只允许特定单环叔胺 N 的可交换构型，因此这个双环 N 被排除。阻断点在质子化与立体规则衔接，不能写成“来源结构缺失”，也不能随意指定一个 N 构型后算作原协议复现。[ChEMBL 母体记录](https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL1201046.json)

可执行后续动作是另立版本审查这些位点的质子化与构型处理，再按同一新版本重算候选与参照；本次保留原规则和既有结果。3 个 ZINC 的非环状叔胺问题与 ceforanide 的双环酰胺问题分别记录。

## 可直接使用的文件

- `author_tree.tsv`：作者冻结树，共 41 文件。
- `author_training_unique12.csv`：12 个训练配体身份、重复次数与原文件定位。
- `paper_13_plus_3_author_source_map.csv`：论文 13 + 3 名称与作者结构文件覆盖映射。
- `ZINC_primary_identity_verification.csv`：5 个请求编号、3 个成功结构与 2 个空缺。
- `resolved_ZINC_primary_identity_2D.sdf`：3 个已解析数据库分子的二维 SDF。
- `resolved_ZINC_reference_inputs.csv`：固定流程输入。
- `raw_runs/zinc_reference_followup3/preparation_ledger.csv`：3 次真实制备记录，全部停止，无新对接分数。
- `ceforanide_atom_diagnosis.json`：原子 15 的电荷、邻接与环归属。

原 R1/R2/R3 科研文件、旧对接结果及阈值均未改动。本次结构来源补全不代表 68 条参考库已经补齐，也不代表项目已有新的实测活性。

论文来源：[Identification of druggable binding sites and small molecules as modulators of TMC1](https://pmc.ncbi.nlm.nih.gov/articles/PMC12075566/)，DOI 10.1038/s42003-025-07943-x。
