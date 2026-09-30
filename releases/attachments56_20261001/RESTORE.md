# 换设备继续运行

## 1 获取完整项目

```sh
git -c core.longpaths=true clone https://github.com/JoyceLeo326/life-science-challenge-tmc1-2026.git tmc1
cd tmc1
git config core.longpaths true
cd releases/attachments56_20261001/code
python verify_manifest.py
python make_results.py --output ../my_results/results.csv
```

Windows建议选择较短的克隆目录。最后一个入口使用Python标准库核对已有证据并导出六候选，不会重跑全部对接。数据、脚本、正式生成模型权重及历史运行记录已在代码目录中。

## 2 安装科研环境

实际验收环境为Windows、Python3.12.10、CPU。请使用独立环境；不需要CUDA。

```powershell
python -m venv ../tmc1_env
../tmc1_env/Scripts/python.exe -m pip install -r requirements.txt --extra-index-url https://download.pytorch.org/whl/cpu
../tmc1_env/Scripts/python.exe -m pip check
../tmc1_env/Scripts/python.exe run_example.py --train-small --output ../my_example
```

Linux或macOS把环境解释器路径换为 `../tmc1_env/bin/python`。完整科学环境的实测记录是Windows，不把标准库跨平台检查等同于所有科学依赖跨平台实测。

## 3 运行Notebook与新增科学表

```powershell
../tmc1_env/Scripts/python.exe -m pip install -r notebooks/requirements-notebook.txt
../tmc1_env/Scripts/python.exe notebooks/run_notebook.py
../tmc1_env/Scripts/python.exe data/group_supplement/verify_science_supplement.py
../tmc1_env/Scripts/python.exe data/group_supplement/rebuild_science_tables.py --output ../my_science_tables
```

Notebook执行训练小例、正式结果入口、617项补充核查和四张表重建，输出写入包外 `../notebook_outputs/`。也可在Jupyter或VS Code中打开 `notebooks/TMC1_walkthrough.ipynb`。

11分子的ADMET-AI预测有单独锁定环境，重跑时依照 `data/group_supplement/REPRODUCE.md`；不要将其依赖混装进本节科研环境。已保存的原始预测和网站输出可直接离线核对。

## 4 恢复外部结构并完整重放

上游受体坐标通过固定公开来源在本机依法取得，不装入公开仓库。执行前先阅读 `reference_restore/README.md`。

```powershell
../tmc1_env/Scripts/python.exe -m pip install -r reference_restore/requirements-restore.txt
../tmc1_env/Scripts/python.exe reference_restore/restore_m_pub_all.py --output-root ../reference_restore_local
../tmc1_env/Scripts/python.exe reference_restore/restore_human_af.py --output-root ../reference_restore_local
../tmc1_env/Scripts/python.exe replay_integrated.py --inputs ../reference_restore_local/EXTERNAL_INPUTS.local.json --work ../my_full_replay
```

如果已有两份固定原始结构，代码README提供 `--source-pdb` 和 `--human-pdb`离线方式；该方式已实际恢复四个输入且SHA匹配。完整重放包含登记的24阶段，具体重新计算与读回范围见 `REPLAY.md`。若需要六个完整复合物，使用 `assemble_candidate_complexes.py`及 `submission/complex_provenance/`，输出保留在新的本地目录。

## 5 接着修改并共享

从最新main建立任务分支。新结果另建日期目录并记录输入版本、命令、运行环境和输出哈希。每批可用成果完成后及时推送，再通过Pull Request合并，由另一名成员核对。原R1/R2/R3和本次已验收目录作为冻结基线保留。

账户配置、聊天原文、签名及本机路径留在本地。正式成员分工与签章由负责同学处理；科研和材料工作可以继续推进。
