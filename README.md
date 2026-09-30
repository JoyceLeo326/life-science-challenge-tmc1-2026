# TMC1 生命科学挑战赛 · 科研与协作仓库

这里保存可运行代码、科学证据、结果、中文材料和协作规范，供团队直接接手。

## 最新交付 · 2026年10月1日

**靶向TMC1机械转导通道的小分子变构调节剂AI辅助筛选——用于Usher综合征1B的听力代偿保护策略**

- [完整成果入口](releases/attachments56_20261001/README.md)：最新代码、中文数据、报告、PPT、图件和视频。
- [完整代码与材料下载](https://github.com/JoyceLeo326/life-science-challenge-tmc1-2026/releases/tag/attachments56-20261001)：两份独立ZIP；代码包68.08MB。
- [换设备运行](releases/attachments56_20261001/RESTORE.md) · [代码目录](releases/attachments56_20261001/code/README.md) · [已执行Notebook](releases/attachments56_20261001/code/notebooks/TMC1_walkthrough.ipynb)。
- [中文科学数据表](releases/attachments56_20261001/data_chinese/A1_A3_科学补充数据.xlsx) · [当前进展和待办](collaboration/CURRENT_STATUS.md)。

最新材料已整合群内A1—A3要求，包含18页中文配音PPT（约7分29秒）和独立真实代码视频（2分23.71秒，8.49MB）。旧R1/R2/R3研究证据保留，新版本以 `releases/LATEST.json`为准。

## 现在如何接手

```sh
git -c core.longpaths=true clone https://github.com/JoyceLeo326/life-science-challenge-tmc1-2026.git tmc1
cd tmc1
git config core.longpaths true
git switch -c work/my-task
cd releases/attachments56_20261001/code
python verify_manifest.py
python make_results.py --output ../my_results/results.csv
```

本轮交付直接使用下方最新文件。如需继续修改，从最新main建立分支，按[协作细则](collaboration/TEAM_COLLABORATION.md)及时推送并互审。科学环境安装、Notebook和完整重放步骤见[接手说明](releases/attachments56_20261001/RESTORE.md)。

公开读取无需邀请；直接推送需要写入权限。尚未获得权限时可以先克隆或Fork。

## 后续协作

已有[三人协作细则](collaboration/TEAM_COLLABORATION.md)和[A科学](work/a_science/README.md)、[B计算](work/b_compute/README.md)、[C材料](work/c_delivery/README.md)工作区保留供后续使用。当前已交付成果直接复用，新增改动通过任务分支和Pull Request合并。

## 文件组织

- `releases/attachments56_20261001/`：最新完整代码、材料、中文数据和运行视频。
- `science/`：此前已发布的R1/R2/R3科研记录，保留用于追溯。
- `collaboration/`、`work/`：三人分工、交接标准、当前进展与各自工作入口。
- `releases/LATEST.json`：当前版本及代码路径。
- `PUBLIC_FILE_MANIFEST.json`与各版本清单：公开文件大小及SHA256。

克隆后可以用 `python scripts/verify_files.py`核对公开文件，或用 `python scripts/check_latest_release.py`在临时副本核对最新正式结果。仓库的自动检查覆盖文件完整性和标准库结果入口；完整科学环境实跑范围见代码测试报告。

## 科学范围

项目完成的是TMC1虚拟筛选、计算分析和候选整理。原六候选、历史分母、失败项与负结果均保留。新增公开药理记录、SwissADME/ADMET-AI预测及SA结果都有来源；预测不属于团队湿实验。68条参考库的完整成员来源尚未建立，完整膜MD、功能方向及听力代偿效果仍需后续验证。

成员真实分工、签章和官网正式提交由负责同学办理。文件发布、代码运行和比赛正式提交分别记录。

## 跨设备与隐私

公开仓库不放原始聊天、私聊、个人信息、凭据和设备配置。外部受体按固定来源与哈希在各自设备恢复。每次工作前拉取最新提交，完成后及时推送，避免成果只保存在一台电脑上。
