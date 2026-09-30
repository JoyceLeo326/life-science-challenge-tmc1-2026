# Amiloride 准备阻断独立只读核查

日期：2026-09-30。范围：只读既有源SMILES、准备日志、源码及公开数据库记录；未运行分子准备、枚举、对接、MD或模型训练，未修改原结果。

**结论：本次是来源SMILES互变异构表示中的潜在 C=N 双键触发保守准备门限。不是未指定四面体手性中心，也不是Dimorphite枚举后产生的问题；现有证据不足以把它定性为药物本身存在独立、稳定、未明确的立体身份。** 应报告为管线技术阻断，同时保留该分子未准备/未对接的事实。

## 本地直接证据

路径均相对 `01b_new_experiments/release_dir/`：

- `raw_runs/reference_panel/logs/REF_CHEMBL945_preparation.json`：CHEMBL945，amiloride；源和规范化输入SMILES均为 `N=C(N)NC(=O)c1nc(Cl)c(N)nc1N`；InChIKey为 `XSDQTOBWRPYKKA-UHFFFAOYSA-N`。
- 日志第17行明确记录 `source_unspecified_potential_stereo = [["Bond_Double", 0]]`，没有 `Atom_Tetrahedral`；第18–20行记录blocked及异常；准备耗时0.001秒。日志没有rule_variant_count、rule_smiles或构象生成字段。
- `protocol_inputs/chemistry_rules.py:11-14` 直接将 RDKit `Chem.FindPotentialStereo` 中标记为 `Unspecified` 的元素返回；没有检验互变异构、E/Z互变速率、溶液稳定性或能垒。因此“potential”不是“已确证稳定立体异构”。
- `protocol_inputs/dock_batch.py:124-127` 对源输入潜在立体元素直接阻断；**第128行以后才调用 `protonate_smiles`**。所以此例没有走到Dimorphite pH7.4规则或枚举，不能归因于Dimorphite误枚举。
- 按该SMILES书写顺序，第0键是开头的 `N=C`，位于酰基胍部分；这是潜在双键立体元素，不是“未定手性碳”。阻断也不是Vina分数或QC失败。

日志仍另有 `database_parent_identity_article_salt_and_batch_not_confirmed` 范围说明：数据库母体身份不等于原文章盐型/批次身份。这个独立来源限制应保留，不能与此处的潜在双键门限混为同一问题。

## 公开数据库交叉证据与判断边界

[PubChem Amiloride，CID16231](https://pubchem.ncbi.nlm.nih.gov/compound/Amiloride) 的当前网页检索记录给出相同InChIKey，规范SMILES为 `C1(=C(N=C(C(=N1)Cl)N)N)C(=O)N=C(N)N`，已定义/未定义原子立体中心和键立体中心计数均为0。

该规范表示将双键放在酰基胍另一位置，其双键碳接两个相同NH2支链；与ChEMBL输入中终端 `N=C(N)N–C(=O)` 的形式不同，但数据库标准InChIKey一致。这支持“检测依赖互变异构/形式表示”的解释，不支持“amiloride存在未指明四面体手性中心”的解释。

注意：PubChem字段来自本轮主来源网页检索；直接动态网页打开只能返回JavaScript提示，PUG REST经网页工具读取失败，未伪称已取得实时原始API JSON。本核查没有独立计算任何互变异构平衡或旋转能垒，不能进一步断言该形式C=N在所有条件下没有E/Z问题；也不能在冻结协议中静默换互变异构形式或放宽门限。

## 建议报告措辞

中文：

> Amiloride（CHEMBL945）在来源SMILES的酰基胍形式中被RDKit标记为潜在未指定C=N双键立体元素（Bond_Double 0），触发冻结准备流程的保守来源检查，因此本批未完成配体准备与对接。该技术阻断发生在质子化枚举之前，不代表该药物具有未明确的四面体手性中心或已证实的稳定立体身份歧义。数据库母体与文献盐型/批次的对应关系仍待独立核实。

English:

> Amiloride (CHEMBL945) was not prepared or docked because the frozen source-input gate flagged a potential unspecified C=N double-bond stereo element in its supplied acylguanidine tautomer representation (Bond_Double 0). The gate stopped before protonation enumeration. This is a representation-dependent preparation limitation, not evidence of an unspecified tetrahedral center or an experimentally established stable stereochemical ambiguity of amiloride. Article salt/batch equivalence remains independently unconfirmed.

原README的“unspecified potential source stereochemistry”可以逐字对应日志，但需要补充上述发生阶段和表示依赖性，避免读者理解为药物天然身份不明确。原blocked记录、零新增对接及全部冻结数据应保留。本核查不提出重跑或门限变更。

## 本地证据SHA256

- `REF_CHEMBL945_preparation.json`：`80e2639d987316465cd2741fb077b87657e66b133cb8222d3026d90df9f26a5b`
- `dock_batch.py`：`b5cc1fd410f85d7251ed4343307ac10f3aab91590052f2de4aef4e7c49a6f543`
- `chemistry_rules.py`：`f595cdc49a82ef6a0305677f2c92b60f4ee38e1ac21a97f54b93fd04a9e2a535`
