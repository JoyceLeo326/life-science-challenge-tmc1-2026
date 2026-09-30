# 最新交付包的自动复现检查

检查入口从 `releases/LATEST.json` 读取当前版本。更新交付版本时，修改该指针及对应清单即可，不必修改检查脚本中的日期。

## 只用 Python 标准库复现正式结果

在仓库根目录运行，建议 Python 3.12：

```bash
python scripts/check_latest_release.py --report-dir verification/latest
```

程序实际执行交付代码的清单校验和 `make_results.py`，再把生成的 `results.csv`、`summary_metrics.csv`、`r3_summary.csv` 与冻结版本逐字节比较。临时结果写到代码包外，并核对代码包运行前后没有变化。

这一步使用已经保存的计算证据重建结果，不需要下载受体或安装对接工具。

## 同时复现 A1—A3 补充表

```bash
python -m pip install -r scripts/requirements-ci.txt
python scripts/check_latest_release.py --with-supplement --report-dir verification/latest
```

补充检查使用 NumPy 和 RDKit，调用代码包原有的两个入口：

- `data/group_supplement/verify_science_supplement.py`：核对保存的数据及来源关系。
- `data/group_supplement/rebuild_science_tables.py`：从包内来源响应重新生成药理记录、相似度配对、成药性汇总和作者训练集交集表，并比较全部共享字段。

依赖安装完成后，这两种检查都可离线运行。它们不执行新的分子对接、模型训练、ADMET 网络请求或分子动力学模拟。完整计算复放与 Notebook 实际执行的既有证据保留在交付包的 `validation/` 和 `notebooks/` 中。

## 自动检查和回执

仓库的 `Validate public delivery` 工作流在 Windows 和 Ubuntu 上执行相同检查，并上传 `latest-release-check-<操作系统>` 回执。克隆仓库、安装依赖和上传 CI 回执需要联网。

回执包含当前版本、代码包来源 SHA256、实际命令、退出码、耗时、结果一致性和源文件是否改变。日志中的临时目录和环境路径使用通用标签替代。

检查工具本身的测试可单独运行：

```bash
python -m unittest discover -s scripts/tests -p "test_*.py" -v
```

测试用的微型数据只验证检查器行为，科研结论以正式交付数据的实际复现回执为准。

## 本次实际验证范围

`docs/reproduction_evidence/` 保存了 Windows Python 3.12.10 上对当前冻结代码的实际执行回执。Linux 工作流已经配置，是否在 GitHub 的当前提交上通过，应以该提交的 Actions 记录为准。

这些检查验证的是交付内容的完整性、已有结果的可重建性与补充来源表的一致性；不把代码通过自动检查等同于新的实验发现或湿实验验证。
