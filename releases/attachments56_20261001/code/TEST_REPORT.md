# 本轮新增验收（2026-10-01）

本节与下方历史测试分别记录。正式6候选CSV和核心运行入口沿用已验证字节。

- 公布基线在独立Python3.12.10环境完成24个不同阶段，累计520.64秒；4018个源文件前后未改。完整证据见 `validation/cleanroom_20261001/PUBLIC_FULL_REPLAY_REPORT.md`。曾出现的来源账本错误及修复已留痕。
- `--source-pdb` 与 `--human-pdb` 备用方式在新的私有输出目录真实完成鼠源、人源恢复，分别21.219秒、8.937秒；两次进程均退出0，四个输入SHA与历史完全一致，坐标仍保存在包外。证据见 `validation/offline_source_restore/OFFLINE_RESTORE_VALIDATION.json`。
- 新版代码运行录屏为143.71秒、8,494,463字节，真实小样本训练与正式结果读回均退出0，6个核心入口/输入/结果哈希与当前工作树一致；证据见 `validation/recording_20261001/`。视频独立提交。
- 当前代码小样本训练重新实跑31.484秒，规则枚举的7行结果与冻结CSV完全一致；见 `validation/group_revision_tests/SMALL_EXAMPLE_RECEIPT.json`。两类神经小例与规则演示、正式六候选分别解释。
- A类126个文件纳入 `data/group_supplement/`。新增只读校验入口实跑1.047秒、离线重建实跑1.297秒，均退出0；622条药理活动、1458类似物配对、11分子成药性总表、6候选交集表共享字段全部0差。原始来源响应、网站输出及451个预训练模型预测值分别保存，异常值未删；见 `validation/group_revision_tests/SCIENCE_ENTRY_RECEIPTS.json` 和 `OFFLINE_REBUILD_RECEIPT.json`。
- A1来源/版本、A2药理、A3及公开副本独立核验合计85/85通过。A1仅在已确认的作者12个训练母体范围内完成比较，68条完整来源仍未建立，报告如实保留此范围。回执见 `validation/group_revision_tests/PUBLIC_SCIENCE_INDEPENDENT_ACCEPTANCE.json` 和 `A2_PHARMACOLOGY_A3_INDEPENDENT_CHECK.json`。

- 新增 Notebook 两次执行全部8个代码单元均通过。首次完整46.354秒；命令入口复跑中小例训练32.135秒、正式读回2.180秒、617项补充检查0.992秒、四表重建1.320秒。正式六候选CSV哈希不变，所有新输出在包外；见 `notebooks/NOTEBOOK_EXECUTION_RECEIPT.json` 和 `NOTEBOOK_RUNNER_REPLAY.json`。

---

## 以下为2026-09-30的原测试记录

# 真实复现测试报告（2026-09-30）

本次实际重跑完整16698排序、60项主动选择和360项既有结果评价。九份科研CSV与原成果逐字节SHA相同，见CORE_REPLAY_CONSISTENCY.json。六候选清单与独立QC、四条件确认、SMILES和六个配体SDF逐项绑定；统计另列summary_metrics.csv。

可移植输入的准备、初段154+124合并、主动27+14合并、14项缺口冻结、41新物理ID分区、完整评价均实际通过。公开副本路径编辑采用固定哈希映射，未去掉门禁。初次测试发现恢复支持文件被旧副本覆盖、续跑job账与QC/接触表缺失，均补齐后真实重测。失败记录保留在ACTUAL_TEST_RECEIPTS.json，不能把失败次与最终成功次混算。

原名单LIB_RNNSZIOJHWINFW实际重新制备/对接/导出：-5.659 kcal/mol与原分数相同，姿势SHA也完全相同。没有重复319项全部对接。24确认结果是独立QC逐条复算和审阅的既有结果；不是本次重跑24次对接。

三个绘图入口和登记入口实际通过。正式SELFIES权重实际安全加载，结构51763参数、67词表，前向输出有限；原20000采样行数确认。正式10epoch训练和20000采样没有重复。小例实际2048条输入、SMILES/SELFIES各1epoch/32采样，性质合格数均0，未放宽门槛。原40生成联合达标0/40保留。

原匿名清单实核3099/3099，当前科研公开清单实核3163/3163。群里另一版本3121/3121没有对应的本地同版清单，本文不沿用该数字。最终交付包采用本包重建MANIFEST与SHA256SUMS，实际项目范围以该清单为准。

脚本、命令参数、开始时间、实际退出码、日志原SHA与关键输出SHA见ACTUAL_TEST_RECEIPTS.json；完整原始含本机路径的内部日志不外发。路径在发布测试回执中用具名变量映射，科学数值保持。ENTRY_TEST_MATRIX.csv列每个本轮声明可用的入口和测试范围。

录屏为匿名网页执行真实Python子进程，展示候选输出；未使用动画冒充运行。16.88秒，H.264 MP4，340100字节。录屏展示已有计算评价和候选结果生成，不展示全对接或正式训练。
