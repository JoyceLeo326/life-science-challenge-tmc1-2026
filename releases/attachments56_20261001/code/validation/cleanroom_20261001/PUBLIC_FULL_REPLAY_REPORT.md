# 公开代码包全链路独立环境验收

## 结论

公开派生代码在新建独立 Python 3.12.10 环境中，完成全部 **24 个不同阶段** 的实际运行。所有预设结果比对通过，原六候选和冻结科研结果未改变。

源包 SHA256：`55100dd627da724f0dda6b1135ef26c03511bd29a7223b9d1ce4dae13289e7d1`。

本次在原物理电脑上的全新软件环境验收，使用新下载、新生成的外部输入；不借用旧受体文件，不称为另一台物理电脑实测。

## 外部输入完整恢复

1. 从原作者固定 Git 提交链接和 AFDB 固定 v6 链接重新联网下载。两份原始 PDB 的 SHA256 均与冻结值一致；标准 TLS 验证开启。
2. 按包内原恢复入口重建鼠源模型对齐 PDB 与 PDBQT，二者均与历史 SHA256 完全一致。
3. 补充 `reference_restore/restore_human_af.py`。其读取原始脚本的同一组函数，以相同参数处理人源 AFDB 模型。人源对齐 PDB 与 PDBQT 均与历史 SHA256 完全一致。
4. 新脚本生成四输入相对路径配置 `EXTERNAL_INPUTS.local.json`。原始和派生坐标仅保存在包外本地目录。

AFDB 官方[FAQ](https://alphafold.ebi.ac.uk/faq)和[许可说明](https://alphafold.ebi.ac.uk/assets/License-Disclaimer.pdf)说明预测数据适用 CC-BY-4.0。鼠源作者模型的明确再分发授权仍未获得，因此此公开包交付恢复入口，不打包作者坐标。

## 运行顺序

以下从公开代码根目录执行，先进入新建的独立 Python 3.12.10 环境：

```powershell
python -m pip install -r requirements.txt --extra-index-url https://download.pytorch.org/whl/cpu
python -m pip install -r reference_restore/requirements-restore.txt
python reference_restore/restore_m_pub_all.py --output-root ../reference_restore_local
python reference_restore/restore_human_af.py --output-root ../reference_restore_local
python replay_integrated.py --inputs ../reference_restore_local/EXTERNAL_INPUTS.local.json --work ../full_replay
```

本轮实际使用相同入口，将输出目录写在隔离测试目录的包外相对路径；逐项真实命令保存在 PUBLIC_ACTUAL_COMMANDS.json 和 PUBLIC_ACTUAL_REPLAY_RECEIPTS.json。受体恢复依赖额外包含 Biopython 1.88 和 Gemmi 0.7.5，其余锁定依赖与代码包一致。安装前 numpy/torch 均不可见，用户和系统 site-packages 未被继承。安装后两次 pip check 均通过。ENVIRONMENT_FULL_LOCK.txt 记录此次 Windows 环境全部 35 个已安装发行包的版本，便于再次恢复。

## 覆盖范围

- 正式 6 候选与 R2/R3 结果读回。
- 16698 项前向排序、主动选择 60 项、原始批次/续跑账本合并、360 项评价。
- SMILES/SELFIES 各 1 epoch 的真实小例训练和采样、7 项教学设计。
- 三组既有结果科研图重新生成。
- 9098 分子双表示结构提示重新计算。
- 6 个 ANM 网络与 24 个 Cα 位移重新计算、ANM 图重新生成。
- 26 个静态姿势记录和 17 个既有 R3 对接任务的独立几何/结果读回。

本次未新增 319 次对接、正式生成模型 10 epoch 重训、20000 样本生成、膜环境 MD 或湿实验。评分、结构提示和粗粒度位移不构成生物活性验证。

7 项教学设计来自独立规则枚举，其 15 列结果与冻结小例及新环境短链路逐字节一致。1 epoch 神经生成小例不承担正式筛选结论：SMILES 32 个采样均无效；SELFIES 32 个采样中 31 个有效、29 个唯一有效，但理化约束后合格数为 0。两个真实训练权重文件已产生；这些生成小例未贡献正式六候选，也未被计为新的有效药物。

## 发现并修复的问题

首次运行在 R3 来源校验处失败：公开版将 ANM 复现说明中的本机路径改成相对路径后，内层整合来源账本仍保存旧版说明的集成哈希。已明确记录新公开文件的集成哈希，并保留原 source_sha256 / copied_sha256。同步内层 ANM 清单和外层清单后，原校验脚本不作改动，全链路重新通过。

首次失败的退出码、子进程回执和日志保留。之后使用 `--resume` 重跑失败阶段并完成后续阶段；没有跳过失败校验，没有改变科学数值。

## 每个阶段的实际耗时

| 阶段 | 秒 | 退出码 |
|---|---:|---:|
| official_candidates | 0.846 | 1 |
| official_candidates | 2.103 | 0 |
| prepare_check | 0.252 | 0 |
| prepare_actual | 0.656 | 0 |
| rank_all16698 | 306.452 | 0 |
| rank_preflight | 0.619 | 0 |
| active_acquire60 | 26.706 | 0 |
| active_physical | 2.935 | 0 |
| initial_merge | 3.139 | 0 |
| active_resume | 2.999 | 0 |
| active_merge | 3.026 | 0 |
| evaluate_full360 | 4.481 | 0 |
| small_example | 33.306 | 0 |
| plot_F8 | 2.714 | 0 |
| plot_stages | 6.513 | 0 |
| plot_confirmation | 2.33 | 0 |
| chemical_domain9098 | 0.746 | 0 |
| chemical_alerts9098 | 97.58 | 0 |
| ANM_check | 0.189 | 0 |
| ANM_actual | 12.203 | 0 |
| ANM_figures | 6.585 | 0 |
| structure_raw26 | 1.738 | 0 |
| reference_top10 | 0.815 | 0 |
| reference_reference_panel | 0.864 | 0 |
| reference_denatonium | 0.578 | 0 |

成功的不同阶段共 24，保留早期失败子进程 1 次。

验收时公开源树 4018 个文件在重放前后逐字节不变。代码根清单 SHA256：`2f6545857c6874dd0ad252e17280fca4a9a92b5fabbd14444fe37b5ace17e271`。

## 证据入口

PUBLIC_FULL_STATUS.json、PUBLIC_INTEGRATED_REPLAY_REPORT.json、PUBLIC_ACTUAL_REPLAY_RECEIPTS.json、PUBLIC_METADATA_REPAIR.json、FRESH_SOURCE_FETCH.json、FRESH_MOUSE_RESTORE.json 和 FRESH_HUMAN_RESTORE.json。所有对外记录使用相对路径或占位符，未包含坐标、凭据、私聊和本机账户信息。
