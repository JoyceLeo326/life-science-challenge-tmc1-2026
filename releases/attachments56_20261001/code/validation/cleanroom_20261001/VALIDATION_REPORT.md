# 附件5独立环境安装与运行核验

结论：**PASS**。冻结代码包无需修改。

## 测试对象与隔离方式

- 代码 ZIP SHA256：`55100dd627da724f0dda6b1135ef26c03511bd29a7223b9d1ce4dae13289e7d1`
- Python 3.12.10，Windows CPU 环境。
- 新建独立 venv，`include-system-site-packages=false`，用户 site-packages 禁用；安装前 numpy 与 torch 均不可见。
- 按包内 requirements.txt 的锁定版本安装；pip 使用官方 PyTorch CPU 索引作为附加源。部分 wheel 来自本机 pip 下载缓存，但未复用其他环境的已安装库。
- 安装耗时 184.095 秒，`pip check` 通过。
- 在同一台物理电脑上的全新软件环境进行。已验证本包安装与相对路径运行；未称为另一台物理电脑实测。

## 实际运行

从解压后的代码根目录，使用新环境执行：

```powershell
python run_example.py --train-small --output ../../outputs/example_results
python make_results.py --output ../../outputs/formal/results.csv
```

- 小例耗时 22.348 秒，包含规则设计、SMILES 1 epoch 训练、SELFIES 1 epoch 训练与小规模采样；3 个真实子进程均退出 0，输出 7 个教学设计。
- 正式结果入口耗时 1.378 秒，输出 6 个冻结候选，并核算已有 R2/R3 证据。
- results.csv、summary_metrics.csv、r3_summary.csv 与包内冻结基线逐字节一致。
- 解压代码树全部 4016 个文件在执行前后逐字节不变，输出写在代码包外。
- 本次未重跑 319 个对接任务、正式 10 epoch 生成训练、完整 16698 前向排序或膜环境 MD。

## 证据

- ACTUAL_COMMANDS.json：实际命令、开始时间、耗时与退出码。
- RESULT_CHECKS.json：逐项结果比较。
- isolation_before_install.log：安装前独立环境状态。
- declared_dependency_install.log、pip_check.log、installed_versions.log：安装与依赖版本。
- training_generation_screening_small.log、formal_six_replay.log：实际入口日志。
- MANIFEST.json：本核验材料哈希。

## 安装版本

| 包 | 实际版本 |
|---|---|
| rdkit | 2025.9.6 |
| numpy | 2.5.3 |
| scipy | 1.18.1 |
| scikit-learn | 1.7.2 |
| dimorphite-dl | 2.0.2 |
| meeko | 0.8.0 |
| torch | 2.5.1+cpu |
| selfies | 2.2.0 |
| joblib | 1.5.2 |
| threadpoolctl | 3.6.0 |
| matplotlib | 3.11.2 |
