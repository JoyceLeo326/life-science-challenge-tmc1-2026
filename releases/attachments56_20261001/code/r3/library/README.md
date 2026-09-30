# 第三轮分子库完整来源与结构警示附件

这里补齐可追溯的公共分子库及9098条冻结候选的结构警示结果，第二轮原库和候选排名保持原版。

## 使用顺序

1. `FILTERS_AND_DENOMINATORS.md`：先看每一步筛什么、分母是多少。
2. `library/CURRENT_MANIFEST.json`：定位当前CSV、正式二维SDF和13条隔离记录。
3. `diagnostics/README.md`与`diagnostics/six_candidate_alerts.csv`：看结构提示及六个候选的逐项结果。
4. `provenance/REPLAY_LIBRARY.md`：查来源与离线重建；`provenance/ARCHIVE_SCOPE.md`解释历史审计对应版本。
5. `independent_review/LIBRARY_ROOT_REVIEW.json`：本附件的文件和分母独立核验。

当前正式SDF有17632条，另13条保持CSV身份并列入隔离。正式SDF的SHA256为46d1073b15e93aca106c5c8ba03ce5d8d1be3537522323c1148671b9712c691a，不使用旧版17619条文件。

9098条的主表示PAINS/BRENK/反应性子集提示为552/3498/447；源表示为649/3568/447；两表示并集为653/3626/447。两种表示均无解析或计算失败，六候选全部未命中这些冻结规则。提示是后续人工复核信息，不能当作实验药效或安全结论。

## 重现及保持原件

本包既有原始化学来源、处理脚本、来源到公开副本的哈希映射，也有本次实际结构警示脚本、冻结规则、输入和运行回执。原库重建入口只做过文件/语法/输入检查，没有再次执行整库重建；其既有独立化学审计与正式SDF SHA完全匹配。

如需重跑结构警示，请先复制整个本包到新的独立工作目录，然后在该副本运行`scripts/run_chemical_alerts.py`，因为原冻结脚本向自身相邻diagnostics目录写结果。不要在唯一归档原件上重跑。不要把旧审计的历史分母当作新诊断分母。

使用`python verify_manifest.py`检查本附件所有公开文件。实验验证、湿实验数据和临床效应不在本次计算附件内。

## 数据署名

ChEMBL data is from https://www.ebi.ac.uk/chembl ; the version of ChEMBL is ChEMBL_37. 数据许可为[CC BY-SA 3.0](https://chembl.github.io/chembl-licensing/)，保留来源与相同许可要求。RDKit规则及源码保留各自原许可，见rules/official_sources和rules/README.md。本次官方许可页面复核日期为2026-09-30。
