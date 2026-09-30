# 附件5/6与科学补充 · 最新完整交付

**靶向TMC1机械转导通道的小分子变构调节剂AI辅助筛选——用于Usher综合征1B的听力代偿保护策略**

2026年10月1日修订。这个目录同时包含可运行代码、正式材料、图件来源、中文数据和实际运行视频。下载整仓库即可接手，独立代码包和材料包也可从[本轮 Release](https://github.com/JoyceLeo326/life-science-challenge-tmc1-2026/releases/tag/attachments56-20261001)获取。

## 直接使用

| 需要做什么 | 打开哪里 |
|---|---|
| 快速理解成果和下一步 | [当前进展](../../collaboration/CURRENT_STATUS.md) |
| 运行代码或换设备 | [接手步骤](RESTORE.md) · [完整代码说明](code/README.md) |
| 一步步查看实际计算 | [已执行 Notebook](code/notebooks/TMC1_walkthrough.ipynb) |
| 阅读完整科研报告 | [Word](materials/documents/科研报告_附件56口径修订版.docx) · [PDF](materials/documents/科研报告_附件56口径修订版.pdf) · [Markdown](materials/documents/科研报告_附件56口径修订版.md) |
| 查药理与成药性 | [补充报告](materials/documents/候选身份、既有药理与ADMET预测补充报告.md) · [中文数据工作簿](data_chinese/A1_A3_科学补充数据.xlsx) |
| 用正式PPT汇报 | [18页配音PPT](materials/documents/靶向TMC1机械转导通道的小分子变构调节剂AI辅助筛选——用于Usher综合征1B的听力代偿保护策略项目介绍（区赛）.pptx) · [PDF](materials/documents/靶向TMC1机械转导通道的小分子变构调节剂AI辅助筛选——用于Usher综合征1B的听力代偿保护策略项目介绍（区赛）.pdf) · [无配音编辑版](materials/documents/靶向TMC1机械转导通道的小分子变构调节剂AI辅助筛选——用于Usher综合征1B的听力代偿保护策略项目介绍（区赛）_无配音编辑版.pptx) |
| 调整讲解或图件 | [逐页讲稿](materials/documents/逐页中文讲稿_附件56修订版.md) · [科研成图](materials/figures) · [数据来源](materials/source_data) · [配音分段](materials/narration) |
| 下载代码运行视频 | [143.71秒MP4](video/靶向TMC1机械转导通道的小分子变构调节剂AI辅助筛选——用于Usher综合征1B的听力代偿保护策略代码运行录屏（区赛）.mp4) · [录屏说明](video/录屏说明.md) |
| 核对六个正式候选 | [冻结候选CSV](code/results.csv) · [提交材料目录](code/submission) |
| 完善正式表单 | [摘要Word](materials/documents/靶向TMC1机械转导通道的小分子变构调节剂AI辅助筛选——用于Usher综合征1B的听力代偿保护策略摘要.docx) · [匿名评审表Word](materials/documents/靶向TMC1机械转导通道的小分子变构调节剂AI辅助筛选——用于Usher综合征1B的听力代偿保护策略评审表（区赛）.docx) |

配音PPT共18页，自动播放约7分29秒。代码视频为2分23.71秒、8.49MB，展示真实小样本训练和正式结果读回；新增 Notebook 补充展示 A1—A3 校验和离线表格重建。视频和代码压缩包分别提交。

## 科学内容

- 原R2仍为360次尝试、319个物理任务；六候选四条件确认共24项。R3结构分析、正常模分析、补充筛选及其负结果原样保留。
- A1核清作者12个训练母体与论文活性集合的区别，并保留三参照物六次N几何敏感性计算。原“68条库”完整成员来源仍未建立。
- A2包括相似度≥0.4的1458对、≥0.6的276对、847个不同ChEMBL ID，以及33个代表分子的622条公开药理活动。六个正式候选属于已知研究化合物，未据此认定为已上市药物。
- A3包括11个分子的真实SwissADME输出、两套SA评分及ADMET-AI的451个端点预测。10个超出范围的原始预测值单独标明；ADMETlab没有形成可用结果。
- 变构调节、TMC1功能方向、USH1B听力代偿和湿实验效果仍是后续验证内容。正常模分析没有被写成完整膜环境MD。

## 已完成的验证

代码基线在同一电脑的新隔离环境完成24阶段实跑，科学CSV保持一致；本轮又实际执行本地源文件恢复、六候选读回、小训练和A类离线重建。Notebook两次执行8个代码单元均通过。A类资料617项本地检查、85项独立科学检查通过。各项范围、用时和失败修复记录都在 [代码测试报告](code/TEST_REPORT.md)。

## 来源与版本

[来源和公开处理](SOURCES_AND_RIGHTS.md)说明再分发范围；[覆盖清单](COVERAGE.md)对应本轮交付内容。`MANIFEST.json` 和 `SHA256SUMS.txt`核对本目录文件，`ASSETS.json`记录完整下载包。科学代码与正式ZIP逐字节对应；材料公开副本只去除评审表的私有Office元数据，正文、表格和图均未改。

仓库已有 `science/` 历史内容仍保留。最新入口由 `releases/LATEST.json`指定，接续工作从本目录开始。
