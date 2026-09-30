# 可运行 Notebook

将本目录内容放在代码包的 `notebooks/` 中。

先在同一 Python 环境安装代码根目录的科学依赖与本目录的 Notebook 依赖：

```bash
python -m pip install -r requirements.txt --extra-index-url https://download.pytorch.org/whl/cpu
python -m pip install -r notebooks/requirements-notebook.txt
```

在 Jupyter 或 VS Code 中打开 `notebooks/TMC1_walkthrough.ipynb`，选择上述 Python 环境并执行全部单元。也可以从代码包根目录直接运行：

```bash
python notebooks/run_notebook.py
```

这个入口会执行整个 Notebook，并把新执行副本保存在代码包外的 `notebook_outputs/`。每次计算输出使用独立时间目录，原始代码和冻结结果不改写。从 `notebooks/` 目录启动 Notebook 也能找到代码根目录。

## 内容

1. 读取 2048 条小样本输入。
2. 调用现有入口实际训练、生成并展示 7 项规则设计示例。
3. 调用正式入口核对六候选，生成实际四条件确认图。
4. 调用现有入口检查 A1—A3，再从保存的官方响应和预测文件重建表格。
5. 保存此次命令、运行时间、结果校验值和执行回执。

7 项规则示例与正式六候选分别保存；两条单轮生成支路的实际合格数直接显示。A1 的 68 库来源缺口、A3 的 10 个范围异常及 ADMETlab 未完成情况保留在说明中。

已执行的 Notebook 自带真实输出，结果和运行耗时见 `NOTEBOOK_EXECUTION_RECEIPT.json`。Notebook 使用包内现有主入口，不复制筛选算法。
