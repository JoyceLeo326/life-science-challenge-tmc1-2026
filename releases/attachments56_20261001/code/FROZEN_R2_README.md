# AI辅助的TMC1小分子候选调节剂筛选与验证

研究副题：TMC1筛选计算证据与USH1B下游补偿的后续假设。

本包保留六个冻结候选、模型源代码与权重、训练记录、输入快照、原始计算账本、独立质检及真实复现测试。对接分数和几何接触属于静态计算指标，TMC1功能方向、人源效应与USH1B补偿尚待实验验证。

## 一键候选结果

Python 3.10或以上；此入口仅用标准库。

```powershell
python make_results.py
```

输出 `results.csv`：六个候选的编号、赛道描述、标准化与规则态SMILES、来源ID、Vina分数、ExtraTrees预测均值与树间标准差、四条件确认读数、模型/运行版本、备注及结构文件名。`structures/*.sdf`为本项目配体坐标，单位Å，采用对接规则态质子化表示。结构不包含上游第三方受体。赛道字段为本项目方法描述，若赛事候选模板规定专用枚举，应据正式模板映射。

`summary_metrics.csv`另存360臂尝试/319唯一物理任务、每臂120、MW≥350匹配分析和24确认任务。程序逐行核分数阈值、独立QC、接触和联合端点，并保留生成初段0/40联合达标的负结果。

## 完整小例与训练

真实测试环境：Windows 64位、Python 3.12.10、CPU；RDKit 2025.09.6、NumPy 2.5.3、SciPy 1.18.1、scikit-learn 1.7.2、Dimorphite-DL 2.0.2、Meeko 0.8.0、PyTorch 2.5.1+cpu、SELFIES发行元数据2.2.0（运行时报告2.1.1）、joblib1.5.2、threadpoolctl3.6.0。CUDA不需要。原历史记录未冻结NumPy/SciPy精确版本，以上为本次实测环境，不能据此声称历史完整锁定。

```powershell
python -m pip install -r requirements.txt --extra-index-url https://download.pytorch.org/whl/cpu
python run_example.py --train-small
```

小例实际运行：两亲本局部规则设计→理化/结构警示约束→固定7项选择→透明的亲本邻域先验排序，并以2048条真实库输入演示SMILES/SELFIES各1epoch、32采样的训练与生成。结果写入独立 `example_results`。小例的先验不属于已校准亲和力预测。1epoch的小规模示例出现0个性质合格生成分子时保留此结果，不替换或放宽标准。正式生成模型使用17645库、10epochs、20000样本、无微调；权重、词表、原始采样与完整训练记录在 `reproduction_bundle/science/R2/evidence/04_design/model_artifacts/final_selfies_17645`。正式模型没有用TMC1活性/对接标签训练。

正式训练入口是 `reproduction_bundle/science/R2/evidence/04_design/method_results/train_selfies_gru.py`；正式前瞻排序、主动获取、对接及评价源代码位于 `03_compute/method`。`run_example.py`在输出目录复制方法，避免改写归档输入。

## 绘图

```powershell
python plots/plot_forward_three_arms.py
python plots/build_verified_stage_figures.py
python plots/plot_from_verified_csv.py F4 data/confirmation_plot.csv figures/confirmation_four_conditions
python plots/register_reviewed_figures.py
```

前三项实际绘图已通过。登记入口只记录实际文件完整性，默认不冒称已获得人工视觉批准。

## 完整科学重放

见 `REPLAY.md` 和 `reproduction_bundle/reproduction/README.md`。准备入口严格校验六项历史输入、冻结方法/支持文件、盒和外部受体。新设备需按固定来源取得受体并恢复；本包不分发原作者受体坐标。正式Windows Vina1.2.7二进制和Meeko导出入口均需显式指定。原子坐标与随机种子控制不代表跨操作系统字节相同。

`TEST_REPORT.md`逐项说明本次实际运行范围。完整排序、主动获取和既有作业评价已重算；319项全对接、正式10epoch模型训练与20000采样没有重复运行。另实际重跑一个原名单任务，并实际加载正式保存权重执行前向计算。既有结果、代表性重跑、小例结果分别报告。

`design_evidence`和`design_independent_audit`补齐04设计账本与05设计复核。`DESIGN_EVIDENCE_SOURCE_INDEX.json`保留原始与发布SHA对应关系；元数据中的机器路径改为逻辑来源，不改科学数值。R1的971→24和R2的6×4属于不同阶段，见 `R1_971_TO_24.md`。

校验入口 `python verify_manifest.py`核对本包MANIFEST/SHA256SUMS。环境安装时间不计入计算运行时间；轻量候选结果通常数秒，小例数十秒，16698全排序本次约213秒。程序结果不构成湿实验或临床结论。

## 正式代码附件与独立录屏

本目录是正式代码与模型提交包。运行录屏作为独立文件提交，媒体及录屏截图不在代码ZIP内；文字回执保留。见FORMAL_CODE_PACKAGE.md和MEDIA_SEPARATION.json。
