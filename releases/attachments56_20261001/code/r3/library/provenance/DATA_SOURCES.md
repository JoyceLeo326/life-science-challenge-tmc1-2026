# 公共分子库来源与使用边界

## ChEMBL

- 官方数据服务：<https://www.ebi.ac.uk/chembl/api/data/molecule.json>
- 官方接口与分页文档：<https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services>
- 官方许可说明：<https://chembl.github.io/chembl-licensing/>；数据许可 CC BY-SA 3.0。下游再分发须保留署名及相同许可条件。
- 下载时官方 `/status.json` 报告 `chembl_db_version=ChEMBL_37`、`chembl_release_date=2026-05-01`；原始响应存 `provenance/chembl_status.json`，SHA-256 `29dc1fb09487a8d2962253a4d0d92c4ded0aa0d831effab42445d0c3483d4740`。
- 本轮请求选择 `molecule_type=Small molecule`，按 `molecule_chembl_id` 排序，每页 1000 条；每页原始 JSON、URL、UTC 获取时间、SHA-256、记录数与失败尝试保存在 `provenance/raw/expanded/`，归档后的相对路径和原字节哈希见 `provenance/summary.json`。该排序可能偏重早期 ID，不代表随机抽样或全 ChEMBL 化学空间。
- `build_library.py` 默认使用上述原始缓存和第一轮只读缓存；移机时可通过环境变量 `TMC1_R2_RAW_CACHE` 与 `TMC1_R1_APPROVED_CACHE` 指向对应目录，再运行 `prepare`。原始页 SHA 必须与冻结 `summary.json` 一致。
- 第一轮 ChEMBL 37 `max_phase=4` 的 3475 条原始记录从只读路径引用，原始文件与哈希已载于本包 `provenance/legacy_acquisition_manifest.json`。本轮合并时按 ChEMBL ID 去重，再按标准化母体键去重，明确区分来源记录数与独立母体数。
- ChEMBL 化合物身份和性质不表示 TMC1 活性。此库仅供计算筛选候选池和无标签化学描述，不以 ChEMBL 其他靶点活性作为 TMC1 标签。

## 处理定义

RDKit 2025.09.6：解析原始 canonical SMILES、标准化、取最大有机片段、去电荷、互变异构体规范化，保留标准化 isomeric SMILES；移除立体信息后生成 InChIKey 首段作为母体键。盐、立体异构形式和互变异构体可在来源行分辨，但不额外增加母体计数。此规则会合并部分实验上可能不同的立体异构体，所以下游需要核查具体立体身份。

性质过滤只是可计算性和常见小分子空间的预筛，阈值及逐项原因见 `build_library.py` 与 `provenance/summary.json`；不能从过滤通过推断药效、耳蜗递送、安全性或通道调节方向。SDF 为二维结构，不包含对接构象。
