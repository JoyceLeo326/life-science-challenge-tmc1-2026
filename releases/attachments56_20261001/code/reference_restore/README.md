# 本轮完整四输入恢复（2026年10月1日）

先在已安装科研依赖的 Python 3.12.10 环境补充结构恢复依赖，再依次恢复鼠源与人源固定输入：

```powershell
python -m pip install -r reference_restore/requirements-restore.txt
python reference_restore/restore_m_pub_all.py --output-root ../reference_restore_local
python reference_restore/restore_human_af.py --output-root ../reference_restore_local
python replay_integrated.py --inputs ../reference_restore_local/EXTERNAL_INPUTS.local.json --work ../full_replay
```

前两条恢复命令从代码包根目录执行。两来源本轮已实际联网下载并核对 SHA256，鼠源与人源转换输出均与历史冻结 SHA256 相同；人源入口在确认四输入齐全后生成相对路径的 EXTERNAL_INPUTS.local.json。受体坐标和本地回执留在代码包外。下方保留早期恢复方法与当时网络失败记录；最新验收见上级交付目录。

新增的人源入口调用相同冻结对齐函数和原 Meeko 参数，不产生新的科学模型、对接或 MD 结果。

### 网络不可用时：使用已下载的固定源文件

将依法取得的两份原始文件放在代码包外的 `../source_inputs/`。鼠源使用作者固定提交的 [MET-channel_complex.pdb](https://raw.githubusercontent.com/ramirezlab/Drug-design-targeting-TMC1/74274dad76201258adbd76c7058b705170c4b857/Mechano_electrical_transducer_channels_complex/MET-channel_complex.pdb)；人源使用 [AF-Q8TDI8-F1-model_v6.pdb](https://alphafold.ebi.ac.uk/files/AF-Q8TDI8-F1-model_v6.pdb)。

| 原始文件 | 必须匹配的 SHA256 |
|---|---|
| `MET-channel_complex.pdb` | `14b5d4010eb2220ebddca0ef5cd0be346e6ba9bcc6db14311745b684e9807d82` |
| `AF-Q8TDI8-F1-model_v6.pdb` | `ba6b438b036b9386f4f04534cfa1dda3b6a35fbe05dd7fb2ec1b415bc6768d3b` |

从代码包根目录运行：

```powershell
python -m pip install -r reference_restore/requirements-restore.txt
python reference_restore/restore_m_pub_all.py --source-pdb ../source_inputs/MET-channel_complex.pdb --human-pdb ../source_inputs/AF-Q8TDI8-F1-model_v6.pdb --output-root ../reference_restore_offline
python reference_restore/restore_human_af.py --human-pdb ../source_inputs/AF-Q8TDI8-F1-model_v6.pdb --output-root ../reference_restore_offline
python replay_integrated.py --inputs ../reference_restore_offline/EXTERNAL_INPUTS.local.json --work ../offline_replay
```

此备用方式已在新的独立输出目录实际完成两次恢复，分别用时21.219秒与8.937秒，四个受体输入均匹配原SHA；见 `../validation/offline_source_restore/OFFLINE_RESTORE_VALIDATION.json`。这是结构输入恢复验证；24阶段整合流程的新环境复现见前述独立验收记录。

---

# M_PUB_ALL 新设备本地恢复



该目录只交付脚本、冻结参数和验收记录。作者原始坐标的明确仓库许可证未在原记录中找到；因此这里不包含作者 PDB、对齐 PDB、加氢 PDB、Meeko 参数化 JSON 或派生 PDBQT。脚本在用户设备的独立本地目录下载和生成这些文件。



## 已完成的验收



2026-09-30，在 Windows 11 x64、Python 3.12.10 和原科研环境中，从已有且 SHA256 匹配的两份输入运行本脚本。新生成的对齐 PDB 和 M_PUB_ALL.pdbqt 均与历史文件逐字节一致。该运行未执行对接、MD 或模型推理。详见 `conversion_acceptance_receipt.json`。



首次固定 URL 在线下载尝试失败，错误为 `SSL: UNEXPECTED_EOF_WHILE_READING`；没有关闭 TLS 验证、替换来源或把在线下载记为成功。详见 `download_attempt_receipt.json`。该次历史验收完成“已有输入 → 原流程结构转换”；本轮在线与离线恢复均已通过，具体证据见本页开头。



## 运行



建议使用 Windows x64 和 Python 3.12.10。以下命令在该交付目录运行。虚拟环境和坐标输出均在独立本地目录；该目录不要加入 Git 或公共压缩包。



```powershell

py -3.12 -m venv ../reference_restore_local\.venv

../reference_restore_local\.venv\Scripts\python.exe -m pip install -r requirements-restore.txt

../reference_restore_local\.venv\Scripts\python.exe restore_m_pub_all.py --output-root ../reference_restore_local

```



`py -3.12` 应选择 **3.12.10**；脚本会检查实际 Python 和依赖版本，版本不匹配会停止。只有取得最终历史 SHA256 后才能将本地生成的 `screen_v2/inputs/receptors/M_PUB_ALL.pdbqt` 用于冻结结果重放。脚本不会自动写入公开仓库。



只下载和核对输入，不需要安装科研依赖：



```powershell

py -3.12 restore_m_pub_all.py --fetch-only --output-root ../reference_restore_local

```



如网络失败而已有合法取得的固定源输入，可传入 `--source-pdb <作者原始PDB>` 和 `--human-pdb <AFDB-v6-PDB>`。两份输入的 SHA256 仍必须完全匹配。缓存不匹配会停止，不会用新来源覆盖它。



脚本默认本地目录是当前用户主目录下的 `TMC1_reference_restore_local`，也接受其他独立路径。它拒绝把输出放入自身交付目录、`completion_staging` 或已有 Git 工作区；会生成 `NOT_FOR_UPLOAD.txt` 和本地 `restore_receipt.json`。



## 固定输入与哈希



作者模型原始来源：



[固定提交的作者 MET-channel_complex.pdb](https://raw.githubusercontent.com/ramirezlab/Drug-design-targeting-TMC1/74274dad76201258adbd76c7058b705170c4b857/Mechano_electrical_transducer_channels_complex/MET-channel_complex.pdb)



提交：`74274dad76201258adbd76c7058b705170c4b857`。



作者原始 PDB SHA256：`14b5d4010eb2220ebddca0ef5cd0be346e6ba9bcc6db14311745b684e9807d82`。



人类对齐参考：[AFDB Q8TDI8 v6 PDB](https://alphafold.ebi.ac.uk/files/AF-Q8TDI8-F1-model_v6.pdb)。SHA256：`ba6b438b036b9386f4f04534cfa1dda3b6a35fbe05dd7fb2ec1b415bc6768d3b`。



历史对齐 PDB SHA256：`a0c2c7b07eb3a4939397ffb104d0d8b8bdc5ef4ab9b7af732deef03944385baa`。



历史 M_PUB_ALL.pdbqt SHA256：`c024cc9efaf910ac11642cea4a4c701a0a88a335d7deb0d040320516d5a2eda1`。



URL 若发生变更、内容改变或失效，脚本会失败。不会悄悄改用最新版模型。



## 原科学方法与实现来源



`original_prepare_receptors_v2.py` 是原 `src/prepare_receptors_v2.py` 的原样源码副本，SHA256 为 `76a534122dd35c0593ceea26fda721439634ef8d91aa53e671d00edc1b9f512a`。它作为方法来源保存；单独执行它需要原项目的其他模块。恢复脚本核对该源码哈希，再通过 Python AST 只加载其中的 `atoms`、`ca_map`、`align` 三个函数，避免执行原项目的其他受体、搜索盒和后续流程。



1. 原函数读取 ATOM/HETATM，去除 H/D；保留原重原子和原顺序。

2. 按 A 链 Cα 残基序列，用 Biopython 全局比对：BLOSUM62，开隙 -10，延伸 -0.5；取首个比对结果。

3. 取人类 UniProt 注释跨膜区内、残基类型相同的 Cα 对，使用原 NumPy SVD/Kabsch 算法对齐到人类 AFDB 参考。将同一刚体变换应用于作者模型全部链，不修改内部几何；PDB 写入 0.001 Å 精度。

4. 冻结的跨膜区间和15个 middle 位点的映射来自原 `data/raw/uniprot_Q8TDI8.json` 和 `results/mouse_human_residue_map.csv`；其原文件哈希与准确数值均在 `frozen_alignment_inputs.json`。原 M_AF 输入和全残基映射用于原项目的其他受体、搜索盒和报告；单独恢复 M_PUB_ALL 不使用它们的坐标，不重新生成公共搜索盒。

5. 调用与原日志完全一致的 Meeko 参数；以同一 Python 的模块入口替代 Windows `.exe` 包装器，当前验收得到相同产物：



```text

python -m meeko.cli.mk_prepare_receptor --read_pdb screen_v2/inputs/receptors/M_PUB_ALL_aligned_heavy.pdb -o screen_v2/inputs/receptors/M_PUB_ALL -p -j --write_pdb screen_v2/inputs/receptors/M_PUB_ALL_H.pdb

```



原日志显示 `Template padding will be used` 和 `gasteiger charges will be read from template file`。没有额外的 pH 参数、质子化状态覆盖、用户模板、外部质子化工具、优化或 `--compute_charges`。补氢与残基化学状态由该版 Meeko 的默认残基模板/RDKit 处理；Gasteiger 部分电荷读取模板；输出为刚性受体。不要把此默认模板过程描述为已在特定 pH 下重新质子化或重新计算电荷。验收回执记录了三份默认模板的 SHA256。



## 依赖与范围



实际验收版本：Python **3.12.10**、Biopython **1.88**、NumPy **2.5.3**、Meeko **0.8.0**、RDKit **2025.9.6**、Gemmi **0.7.5**、SciPy **1.18.1**。精确版本在 `requirements-restore.txt`；urllib 下载仅使用 Python 标准库。



这份历史记录对应原环境运行。本轮另已完成全新 Windows 虚拟环境的在线与离线恢复，Linux/macOS 尚未实测。脚本保留历史 Windows CRLF 对齐 PDB 字节格式，并对最终 PDBQT 作严格哈希判定；跨平台差异若导致不匹配会明确失败。



输出目录还包含 `M_PUB_ALL.json`、`M_PUB_ALL_H.pdb` 和日志。这些参数化/坐标产物也只用于本地研究，不属于公开交付。验收结果不等同于结合验证或门控状态判定。

