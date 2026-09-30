# 靶向TMC1机械转导通道的小分子变构调节剂AI辅助筛选——用于Usher综合征1B的听力代偿保护策略

这是可单独提交的代码与模型附件。解压后进入本目录，根目录具有唯一候选结果入口。整合时间2026-10-01，科研结果日期2026-09-30。R1、R2及R3的冻结原件均保持独立；原六候选CSV与数值未改变。代码与录屏分开提交；本代码附件不含音视频。

## 附件5建议目录与本包实际位置

本包按功能对应官方建议目录，保留科研证据和脚本之间的原有相对路径。

| 官方建议位置或功能 | 本包实际位置与入口 |
|---|---|
| `README.md` / 环境说明 | 本文件、`ENVIRONMENT_ACTUAL.json`、`requirements.txt`；外部受体恢复依赖见 `reference_restore/requirements-restore.txt` |
| `data/` / 数据来源与处理 | `data/`、`r3/library/`、`SOURCE_INDEX.json`、`reproduction_bundle/science/R1/provenance/`；受限结构获取见 `reference_restore/README.md` |
| `src/` / 核心算法 | `reproduction_bundle/science/R2/evidence/03_compute/method/`、`reproduction_bundle/reproduction/` 与 `r3/` 下对应计算脚本 |
| `models/` / 正式权重 | `reproduction_bundle/science/R2/evidence/04_design/model_artifacts/final_selfies_17645/`；模型适用范围见 `MODEL_CARD.md` |
| `train.py` / 训练入口 | `run_example.py --train-small`；底层训练脚本见 `reproduction_bundle/science/R2/evidence/04_design/method_results/train_smiles_gru.py` 与 `train_selfies_gru.py` |
| `design.py`、`screen.py` / 设计筛选 | `run_example.py --train-small` 展示短流程；`replay_integrated.py` 复现已登记完整计算阶段；原设计入口为 `reproduction_bundle/science/R2/evidence/04_design/method_results/run_design.py` |
| `predict.py` / 最终清单入口 | `python make_results.py`，生成正式六候选与汇总文件 |
| `results/` / 结果示例 | 根目录 `results.csv`、`summary_metrics.csv`、`r3_summary.csv`，以及 `example_actual_outputs/` |
| `logs/` / 日志与随机种子 | `logs/`、`validation/` 和各科学模块原执行回执 |
| `notebooks/` / Notebook | `notebooks/TMC1_walkthrough.ipynb`，包含真实输出；`python notebooks/run_notebook.py` 可重跑全部单元 |

代码附件和运行视频独立提交。本包中的三份 `.csv.gz` 是原始采样数据的压缩表示，供原脚本读取；不含嵌套交付 ZIP 或音视频。

## 一键结果（Python 3.10以上，标准库）

```powershell
python make_results.py
```

输出根目录`results.csv`（6个冻结候选）、`summary_metrics.csv`（R2固定预算评价）、`r3_summary.csv`和`R3_READBACK_REPORT.json`（R3单独读回）。入口真实核算360次尝试、319个物理任务、MW≥350匹配评价、四条件确认24次计算，并验证9098分子结构提示、三静态结构矩阵、正常模结果及177补筛的原始证据。候选CSV中模型/运行版本、标准化与规则态SMILES、来源ID、配体SDF、确认条件和备注齐全。`make_results_r2.py`是原冻结入口的原字节副本。

## 完整既有计算复现

本次实际验证主机：AMD Ryzen 7 8845H（8核16线程），可见物理内存27.81 GiB；科学运行使用Python 3.12.10、Windows 11（精确版本见ENVIRONMENT_ACTUAL.json），CPU运行，不需要CUDA。完整复现入口将OMP/MKL/OpenBLAS/NumExpr线程设为1。最低内存需求未单独测量，未声明最低配置。

先在独立Python环境安装依赖：

```powershell
python -m pip install -r requirements.txt --extra-index-url https://download.pytorch.org/whl/cpu
```

2026年10月1日在同一物理电脑新建独立 Python 3.12.10 环境，完成依赖安装、两次 pip check、固定来源在线恢复和24个不同阶段的实际复现，累计520.64秒；三份正式结果CSV逐字节一致。报告和失败修复记录见 `validation/cleanroom_20261001/PUBLIC_FULL_REPLAY_REPORT.md`。这项验收覆盖加入本轮A类资料前的科学基线；本轮新增资料的测试另记在 `TEST_REPORT.md`，不混为一次运行。实际机器配置见 `ENVIRONMENT_ACTUAL.json`。

复现受体为外部输入：本包未获得上游作者受体模型的再分发授权，因此不装入该受体或派生全原子复合物。按`reference_restore/README.md`取得固定来源并恢复；填写`EXTERNAL_INPUTS.template.json`的四个外部文件。程序严格校验SHA，不能换结构绕过校验。

```powershell
python replay_integrated.py --inputs EXTERNAL_INPUTS.local.json --work ../TMC1_replay_new
```

输出目录须新建或为空且在包外。此入口实际重算16698前向排序、60项主动选择、初段/续跑账本合并、360项完整评价、小例训练及三幅图；实际重跑9098 PAINS/BRENK双表示诊断、6个ANM网络与24个Cα位移，以及26项静态姿势和17项此前R3姿势（9生成物、7文献参照、1苯甲地那铵）的独立读回。既有对接物理任务、正式10epoch生成模型训练和20000采样没有在本次整合中全部重新运行。代表性真实对接入口及既有执行回执见`REPLAY.md`、`REPRESENTATIVE_DOCKING_COMPARISON.json`。

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

此备用方式已在新的独立输出目录实际完成两次恢复，分别用时21.219秒与8.937秒，四个受体输入均匹配原SHA；见 `validation/offline_source_restore/OFFLINE_RESTORE_VALIDATION.json`。这是结构输入恢复验证；24阶段整合流程的新环境复现见前述独立验收记录。

## 真实训练、生成和筛选小例

```powershell
python run_example.py --train-small --output example_results
python make_results.py
```

小例为2048条库输入、SMILES/SELFIES各1epoch和32采样；透明规则设计与7项先验排序是独立教学范围。不得把它说成正式模型重训、校准亲和力预测或生成了冻结六候选。正式模型是17645库、10epochs、20000采样、0微调；权重、词表、原始采样和方法位于`reproduction_bundle/science/R2/evidence/04_design/`。原40分子联合0/40和R3补筛10项联合0/10的负结果保留。经典ExtraTrees/Ridge代理使用固定种子、同一935/995标签和冻结特征重建，版本、参数与训练/前向隔离规则见MODEL_CARD；未声称已有历史代理序列化权重。正式SELFIES保存权重、词表与训练配置另有完整快照。六个正式候选均为ChEMBL已有记录的化合物，属于已知化合物的再定位筛选；是否属于已上市药按逐项药品状态说明，不称发现新分子。

## 可直接运行的 Notebook

`notebooks/TMC1_walkthrough.ipynb` 已保存真实执行结果。它调用相同的小例训练、正式六候选读回和新增 A 类证据入口，并生成四条件确认图。已完成两次全单元执行，8 个代码单元均通过；617 项补充核对通过，四张科学表离线重建共享字段零差。执行回执见 `notebooks/NOTEBOOK_EXECUTION_RECEIPT.json` 和 `NOTEBOOK_RUNNER_REPLAY.json`。

在上述科研环境额外安装 Notebook 依赖后，可从代码包根目录直接运行：

```powershell
python -m pip install -r notebooks/requirements-notebook.txt
python notebooks/run_notebook.py
```

也可用 Jupyter 或 VS Code 打开 Notebook 并执行全部单元。运行结果和新执行的 Notebook 保存在代码包外的 `../notebook_outputs/`；源代码和冻结结果不改写。这个演示覆盖小例训练、正式结果核对和 A 类表格重建，不把它作为重新运行全部对接的证据。

## 完整来源及新增科学附件

| 包内路径 | 内容 |
|---|---|
| `reproduction_bundle/science/R1` | 历史输入、971→24证据与固定盒 |
| `reproduction_bundle/science/R2/evidence` | 冻结库、模型、生成样本、方法、原始任务与独立QC |
| `r3/library` | ChEMBL双批次来源/许可/获取记录、安全SDF17632加13隔离、17645理化合格表及9098分子结构提示原件 |
| `r3/anm` | 6个网络、24个Cα位移、全部数值和独立验收、复现入口及ProDy参考源码许可 |
| `r3/structure` | 7USX来源/预处理、24格矩阵、26个原始job、12新增计算及姿势/省略HET碰撞证据 |
| `r3/reference_generation` | 冻结177候选、162可选、top10、9次实际对接及参照物来源、原始结果和复现入口 |
| `submission/` | 与冻结CSV一致的官方六候选Excel和材料来源候选CSV |

`INTEGRATED_SOURCE_SNAPSHOT_INDEX.json`逐项列原项目相对来源、原SHA、包内位置及纳入/排除/改名；`INTEGRATION_TRANSFORMATIONS.json`记录新增入口与元数据整合。短R3目录用于避免Windows深路径，组件内部相对路径保持原样。`MANIFEST.json`与`SHA256SUMS.txt`覆盖本新包；用`python verify_manifest.py`核验。

## 群内A1—A3补充证据

本轮新增内容在 `data/group_supplement/`，与原R1/R2/R3计算分母分别保留。

| 内容 | 直接入口 | 本轮实际范围 |
|---|---|---|
| A1 参照来源与交集核查 | `data/group_supplement/A1_reference_sources/A1_README.md` | 核清作者12个训练母体与论文12个染料活性物的区别；68条完整库的来源仍未建立，未制造68行成员表 |
| A2 六候选类似物与药理 | `data/group_supplement/A2_chembl_neighbors/A2_README.md` | ChEMBL ≥0.4 的1458配对、≥0.6的276配对，按ChEMBL ID去重847项；33个代表分子622条公开实验活动 |
| A3 成药性预测与合成可及性 | `data/group_supplement/A3_ADMET_SA/A3_README.md` | 6候选加5探索/对照，共11分子的SwissADME原输出、SA评分和ADMET-AI 41端点预测；10个越出物理范围的原输出单列保留 |
| 来源、模型版本及复现 | `data/group_supplement/REPRODUCE.md` | A类来源时间、API原响应、网站实际输出、工具与权重哈希；ADMET模型使用单独环境 |

A2的相似度来自ChEMBL官方接口，与原筛选Morgan1024描述符分别解释。六个正式候选均缺少已上市证据，定位为已知研究化合物的新用途候选。A3得到的是多参数成药性预测：六候选中5个胃肠吸收预测为Low，六个BBB预测均为No，不能据此推断血迷路屏障通透性。

新增结果可在原科研环境离线重建：

```powershell
python data/group_supplement/verify_science_supplement.py
python data/group_supplement/rebuild_science_tables.py --output ../group_supplement_rebuilt
```

新入口已在本次工作树真实运行。622条药理活动、1458条类似物配对、11个成药性总表及6候选与作者训练12结构交集的共享字段全部重建一致；见 `validation/group_revision_tests/OFFLINE_REBUILD_RECEIPT.json`。网站实际预测结果从已保存的官方输出读取，不把离线读回称为再次调用网站。

## 结构文件与证据范围

`structures/*.sdf`是规则态配体坐标，单位Å；它们不是完整小分子—靶点复合物。复现者需先依法取得并恢复外部受体，结合`REPLAY.md`、候选CSV和原始姿势文件恢复对应结构条件；候选Excel的“结构与设计”和新增“复合物索引”页逐候选引用六个PDB，实际PDB作为09反馈修订团队私有复合物附件另交；本代码包包含assemble_candidate_complexes.py及submission/complex_provenance/以便用有权使用的外部受体重建，同包不分发作者受体坐标。Excel“冻结源CSV”6×17原值不变。

PAINS/BRENK和冻结反应性子集只是结构提示，不构成毒性、药效或实验反应性证据。ANM是Cα弹性网络，不是MD；24帧全原子动态对接受体资格为0。三物种24格矩阵中蛋白QC通过23/24；矩阵里的线虫8/8首姿势与省略HET配置存在小于1.5Å重叠，扩展到全部线虫原始作业为10/10；保留失败与完整分母，不作跨物种亲和力排序或复合物兼容结论。旧760行全局表保留为历史，局部15位点已纠正，不称全表已校正。68参考库成员仍未核实，不造清单或零重叠结论。功能方向、人源效应及USH1B补偿待真实实验。

## 历史整合测试（保留原运行记录）

24个不同阶段已完成，保留26次子进程记录，其中2次早期绘图目录适配失败后成功；更早的环境、证据路径比较及输出生成顺序错误一并公开说明于validation/INTEGRATION_FAILURE_HISTORY.json。已完成16698排序和主动60选择没有重复运行。小例20.813秒、9098诊断67.570秒、ANM8.036秒。除已披露的证据路径元数据之外，科学单元格严格一致；全部其他指定科研CSV逐字节匹配。
