# 六候选的类似物与已知药理

官方 ChEMBL 两档结构检索共得到 1458 条候选—邻居配对，按 ChEMBL ID 去重为 847 个记录。每个候选均已取尽分页；≥0.6数量另向官方60%接口交叉核对。阈值来自官方相似度百分数除以100，不等同本项目原Morgan1024分数。

|候选|≥0.4|≥0.6|
|---|---:|---:|
|CHEMBL1079964|427|26|
|CHEMBL1080362|258|63|
|CHEMBL1081740|150|62|
|CHEMBL1080669|158|59|
|CHEMBL1076497|292|32|
|CHEMBL1080205|173|34|

六候选输入完整 InChIKey 与官方分子记录一致。本次相似接口未返回同ID、同完整 InChIKey 或同规范SMILES的自结构，不人为补入相似度1的自匹配。`including_self`与`excluding_self`列相同描述的是实际返回集合；未把这一观察写成已查证的接口默认规则。

为回答已知药理问题，按预先登记规则查看每个候选≥0.6前三个非正式候选邻居、0.4至0.6前两个邻居以及全部带临床阶段记录的邻居，连同6候选共33个分子；实际取得622条实验活动、422个assay、48篇来源文献。原始值、关系符、单位、assay、靶标与DOI均保留。完整847记录中未进入这个短名单的行明确标为未查药理，不能理解为无活性。相似结构的活性不能直接赋予候选。

## 对当前候选的解释

五个候选已有EphB4研究背景：CHEMBL1081740、CHEMBL1080669、CHEMBL1080362的记录分别为IC50 850、10500、50000 nM；CHEMBL1079964和CHEMBL1080205是10µM下36%和32%抑制，不能当作IC50。来源为[2009年EphB4原始论文](https://pubmed.ncbi.nlm.nih.gov/19879134/)。

CHEMBL1076497存在B-Raf IC50 44 nM记录，同时有HT-29和WM266.4细胞毒性IC50 740与310 nM。来源为[2009年B-Raf原始论文](https://pubmed.ncbi.nlm.nih.gov/19864136/)。这提示后续需要关注激酶相关作用和细胞安全性，不能仅凭TMC1对接分数推断选择性。

六个候选的ChEMBL首批上市年份、最高临床阶段均为空；现有材料不足以把它们统称为获批老药。建议定位为“已知研究化合物的新用途候选筛选”，不称发现全新分子，也不声称六个都是老药重定位。

数据来源：[ChEMBL API](https://www.ebi.ac.uk/chembl/api/data/docs)；各请求的获取时间、URL、SHA256在raw子目录回执中。ChEMBL来源及本目录派生表按 [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/) 归属与相同方式共享，见[官方许可说明](https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/general-questions)。
