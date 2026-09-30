# 恢复的 R1 历史输入

本目录恢复冻结的历史科研输入，供 R2 便携准备脚本读取。文件保持 R1 相对目录结构和原始字节；没有重新运行对接、模型训练、排序或实验。

## 文件与核对结果

R2 `science/R2/evidence/replay/prepare_portable.py` 所需的六个 R1 输入全部恢复，SHA256 六项全部相符；额外核对该脚本指定哈希的 `M_PUB_ALL` 受体，一项相符；受体不包含在此公开补丁中。核对合计 7/7，暂存科学输入文件为六项。逐项哈希、字节数和行数见 `R1_INPUT_MANIFEST.json`。这里的来源匹配数量指核实的原文件，不表示已统计所有重复副本。

| 输入 | 数据行数 | 用途 |
|---|---:|---|
| `coarse_v1/analysis/job_audit_join.csv` | 971 | 历史作业标签与质控字段 |
| `coarse_v1/inputs/library/selected_candidates.csv` | 971 | 历史化合物库 |
| `confirmation_v1/inputs/library/selected_candidates.csv` | 12 | 设计亲本输入 |
| `confirmation_v1/analysis/parent_setting_complete_denominators.csv` | 72 | 亲本条件完整分母 |
| `screen_v2/inputs/library/eligible_pool.csv` | 971 | 旧 eligible 排除池 |
| `coarse_v1/config/vina_box.txt` | — | 冻结对接盒 |
| `screen_v2/inputs/receptors/M_PUB_ALL.pdbqt`（外部输入，未包含） | — | 原文件已核对；固定来源与处理身份见 `SOURCE_RIGHTS/needs_reference_fetch.json` |

`--r1` 应指向本目录根。R2 便携脚本也要求其 R2 根中的 `03_compute/receptors/M_PUB_ALL.pdbqt`，且要求同一哈希；是否提供该 R2 文件应另行核对。本次恢复没有证明完整 R1/R2 重算输入、软件环境或平台复现全部齐备。

## 已有科学来源与许可证据

化合物表保留 ChEMBL 标识符及公开 `source_url`。`provenance/chembl_approved_pool_source_manifest.json` 是原始下载登记，记录公开 API、2026-09-27 的读取时刻和原下载文件哈希。它不代表重新下载或新增药理验证。

`provenance/receptors.json` 记录 `M_PUB_ALL` 为作者提供的小鼠六蛋白 AF2/MD 坐标模型的处理版本，原坐标 SHA256 为 `14b5d4010eb2220ebddca0ef5cd0be346e6ba9bcc6db14311745b684e9807d82`。它属于模型条件敏感性比较，没有已实验指定的开/闭状态；原输入没有膜脂、离子或水。

`provenance/author_model_provenance.json` 给出 [作者模型仓库](https://github.com/ramirezlab/Drug-design-targeting-TMC1)、提交 `74274dad76201258adbd76c7058b705170c4b857`、源下载入口及原坐标哈希；论文 DOI 为 `10.1038/s42003-025-07943-x`。原登记明确表示未找到仓库 LICENSE，且该模型未收入当时的公开结构资源 ZIP。因此公开补丁不包含 `M_PUB_ALL.pdbqt`。`SOURCE_RIGHTS/needs_reference_fetch.json` 保存固定提交下载 URL、原坐标与目标受体哈希、原处理步骤及依赖哈希，供新设备从原作者入口获取后重建核对。六个 CSV/盒输入可独立使用；原坐标和派生受体的许可限制继续保留。ChEMBL 字段和其他科研结果的来源归属也应保留；本次恢复未作新的法律判断。

## 隐私与实际缺项

本目录及清单不记录本机绝对路径或个人身份。完整文本经本机路径、身份及凭据关键字筛查，并检查 CSV 字段；没有复制私人聊天、凭据、环境目录或通信导出。

本任务指定的六个 R1 输入没有缺失或哈希不匹配项。受体真实原文件已找到且哈希相符，但公开补丁将它列为需从固定来源取得的外部输入。受体上游许可仍属未确立事项，完整重算环境和 R2 其他输入应由总体复现清单独立确认。恢复输入并不等于科研计算重跑完成。

## 补充对接盒元数据依赖

另已恢复原始 `coarse_v1/config/common_box.json`，SHA256 `2330b07cf13cb289f845344898fcda828c2f6f6024be55af8dc75e45c3e52dca`。这是原 `dock_batch.py` 的默认 `--box-meta` 输入，补充到前述六项冻结输入之外。本目录现有七个科学输入文件；原受体仍作为外部输入记录。新增元数据只作原始文件身份恢复，没有重跑计算。
