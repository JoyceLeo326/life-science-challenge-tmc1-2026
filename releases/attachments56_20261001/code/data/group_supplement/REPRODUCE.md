# 补充证据的离线重建与模型重算

## 1 离线复核和总表重建

使用Python 3.12及`requirements_offline_audit.txt`，与本次SA复算版本一致。也可在原代码包已具备同版RDKit的环境中运行。

```text
python -m pip install -r requirements_offline_audit.txt
python verify_science_supplement.py
python rebuild_science_tables.py --output rebuilt_tables
```

重建程序从包内ChEMBL原始JSON重建622条实验记录和1458条候选—邻居配对，从两份SwissADME官方CSV和本地实际ADMET-AI输出重建11分子成药性总表，并重新计算6候选与12训练母体的交集。输出与交付表共享字段逐项比对，收到PASS和`OFFLINE_REBUILD_RECEIPT.json`为通过。只读原证据，输出到指定新目录；不是新做对接。

## 2 ADMET-AI独立可选重算

此模型环境与正式筛选环境分开。实际使用Python 3.11、ADMET-AI 2.0.1、PyTorch 2.8.0 CPU及列明依赖。版本和10个权重的SHA256保存在`A3_ADMET_SA/ADMET_AI_v2_inference_manifest.json`。先复制A3目录到新的工作副本，避免覆写原始结果。

```text
python -m venv .venv_admet
# 使用新虚拟环境中的Python执行以下命令
python -m pip install -r A3_ADMET_SA/requirements_admet_ai_exploratory.txt
python A3_ADMET_SA/run_admet_ai11.py
python A3_ADMET_SA/range_qc11.py
```

官方包和模型来源为 https://pypi.org/project/admet-ai/2.0.1/ 与 https://github.com/swansonk14/admet_ai 。模型资源应位于安装包`admet_ai/resources/models/`，运行清单会对其10个`.pt`文件计算SHA256。缺模型时应从官方包/对应版本恢复并核对清单，不使用其他权重冒充本次结果。入口在导入PyTorch前隐藏CUDA并核对模型实际device，CPU单线程；新设备有GPU也可以执行。本次新版入口已实际重放451值，与首轮最大绝对差0，见`ADMET_CPU_REPLAY.json`。

## 3 网站计算的复核

SwissADME输入为`swissadme_input11.smi`，分别保存了5条与6条完整原始CSV、返回结构身份核验和官方job URL。网站模型未在本项目本地重实现；重新请求可能遇到版本/服务变化，应另存并记录日期，不覆盖冻结值。按官方要求串行提交，不开多个同时任务。ADMETlab本次官方接口404，没有其预测结果。

## 4 来源与分母

6正式候选和5探索/对照角色见`ADMET_SA_panel11.csv`。ChEMBL两档阈值使用API的百分比结果除以100，不能当作原项目Morgan1024相似度。A1中的68来源仍未成立；离线重建只检查可复核的作者12训练结构，不补造68成员。
