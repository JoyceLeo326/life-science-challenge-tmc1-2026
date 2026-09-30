# 冻结9098分子结构提示诊断

已实际执行RDKit 2025.09.6官方FilterCatalog：PAINS 480条（A=16、B=55、C=409），BRENK 105条。反应性提示为冻结JSON指定的10条BRENK子集。全部585条运行条目与冻结官方SMARTS重建后的完整序列化逐一一致。

输入来自冻结R2的16698行预测表，严格取 `docking_eligibility == eligible` 的9098行。compound_id与parent_key均唯一，顺序与原输入逐一一致。主表示为原冻结的 `rule_smiles`；敏感性表示为 `standardized_isomeric_smiles`；同时报告两表示命中并集，不改变库或原排名。

| 表示 | PAINS命中/9098 | BRENK命中/9098 | 反应性子集命中/9098 |
|---|---:|---:|---:|
| 主表示 | 552 | 3498 | 447 |
| 源表示敏感性 | 649 | 3568 | 447 |
| 两表示并集 | 653 | 3626 | 447 |

全部尝试行保留；主表示失败0，源表示失败0。两表示规则集合不同的分子数为311。分母始终为9098，另列可解析分母；失败不能解释为无提示。

| 六候选 | PAINS并集规则数 | BRENK并集规则数 | 反应性并集规则数 |
|---|---:|---:|---:|
| LIB_DNDVGFIEWOHGAH | 0 | 0 | 0 |
| LIB_FCRHJSNRFAAYQD | 0 | 0 | 0 |
| LIB_MYAFDHCCQBIZLC | 0 | 0 | 0 |
| LIB_YWDIXVZCRMMQQY | 0 | 0 | 0 |
| LIB_ZNDXPRCORYDEAP | 0 | 0 | 0 |
| LIB_ZVHBDYFMHCQADD | 0 | 0 | 0 |

具体规则见 `chemical_alerts_9098.csv` 与 `six_candidate_alerts.csv`。每条规则频次（含零命中）见 `rule_frequencies.csv`；所有比例与失败分母见 `alert_summary.csv`；互斥并集模式见 `overlap_patterns.csv`。图 `chemical_alerts_9098.png` / `.pdf` 使用独立条形展示重叠提示，并另列合计9098的互斥模式；可复算图源见 `figure_source.csv`。

自检从产出CSV独立重读重算汇总、规则并集与全部身份。固定种子随机12分子加六候选的两个表示用全部585条官方SMARTS（mergeHs=True）直接核对目录，共21060次比较，零差异。条目检查见 `catalog_runtime_checks.csv`，抽查见 `direct_smarts_spotchecks.csv`。输入、规则、脚本SHA256及真实软件、时间和执行回执见 `execution_receipt.json`。

结构命中是复核提示，不能当作实验反应性、毒性、测定干扰、TMC1生物学活性或疗效证据。BRENK子集包含广义三元杂环，不能仅解释为环氧化物。各提示集合有重叠，禁止直接相加得到分子总数。规则来源、许可与原始引用完整保留于 `../rules/README.md`、`../rules/RULES_FROZEN.json` 和 `../rules/official_sources/`。

复现：使用RDKit 2025.09.6及matplotlib运行 `../scripts/run_chemical_alerts.py`。程序按冻结SHA校验输入和规则；不会写回R2源文件。
