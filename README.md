# TMC1 生命科学挑战赛 · 三人协作仓库

这里集中保存项目代码、计算结果、科研材料与协作进度，供团队在不同设备接手。

## 现在如何接手

1. 克隆本仓库，或从绿色 **Code → Download ZIP** 下载。
2. 阅读 [三人协作细则](collaboration/TEAM_COLLABORATION.md) 和 [科研文件入口](science/README.md)。
3. 在 Issues 认领工作，从最新 `main` 新建自己的分支。
4. 每批成果完成后 5 分钟内提交并推送；用 Pull Request 汇合，由另一位成员复核。

```sh
git clone https://github.com/JoyceLeo326/life-science-challenge-tmc1-2026.git
cd life-science-challenge-tmc1-2026
git switch -c work/my-task
```

公开读取无需邀请。直接向本仓库推送需要仓库负责人添加协作者；添加前可以先下载、克隆或通过 Fork 提交修改。

## 三人分工

| 负责人 | 主要工作 | 完成时提交 |
|---|---|---|
| A · 科学依据与项目逻辑 | 核对文献、明确研究问题、修订科学正文与结论 | 来源表、科学正文、规则要求与证据对应表 |
| B · 计算与复现 | 整理代码与依赖、核对输入输出、复现计算结果 | 运行入口、参数、日志、结果与图表数据 |
| C · 展示与交付 | 按规则排版报告、制作展示材料、整理交付包 | 报告成稿、PPT、讲稿、交付清单与下载包 |

最终整合与发布负责人由三人稍后认领。A、B、C 的具体职责、互审和交接要求见 [协作细则](collaboration/TEAM_COLLABORATION.md)。

直接领取任务：[A 科学](https://github.com/JoyceLeo326/life-science-challenge-tmc1-2026/issues/1) · [B 计算](https://github.com/JoyceLeo326/life-science-challenge-tmc1-2026/issues/2) · [C 展示与交付](https://github.com/JoyceLeo326/life-science-challenge-tmc1-2026/issues/3) · [最终整合待认领](https://github.com/JoyceLeo326/life-science-challenge-tmc1-2026/issues/4)。

每人的具体文件清单：[A 工作区](work/a_science/README.md) · [B 工作区](work/b_compute/README.md) · [C 工作区](work/c_delivery/README.md)。

## 文件组织

- `science/`：按轮次组织的科研材料、计算证据与复现代码。
- `collaboration/`：分工、进展、交付要求。
- `.github/`：提交说明模板。
- `PUBLIC_FILE_MANIFEST.json`：本次公开文件的大小与 SHA256 校验值。

已收录 3,165 份科研文件，包括 14 个 Python 脚本、3 份 Word、4 份 PPT、10 份 PDF，以及数据、结构、计算日志和图件。R2、R3 分目录保存；原始聊天及个人资料不在公开包内。

完整重算所需的部分历史 R1 输入和方法尚未收入公开包，详见 [复现说明](science/REPRODUCIBILITY.md)。B 的任务包含补齐可公开的缺失输入与运行入口。

克隆后可运行 `python scripts/verify_files.py` 核对交付文件哈希；这是文件完整性核对，不会启动科研计算。

## 结果口径

项目以 TMC1 虚拟筛选与计算分析为基础。计算评分、结构预测和候选排序应按各报告的证据范围解释。湿实验属于后续验证工作；正式赛事提交以实际回执为准。

## 跨设备与隐私

项目公开内容通过 Git 同步。原始聊天、个人信息、账号凭据、密钥和设备私有配置不进入本仓库。原有历史迁移档案继续保存在私有库中。公开文件的下载与运行说明随材料一并提供。

在另一台设备继续工作时先拉取最新提交，完成一批内容便及时推送，避免成果只留在本机。
