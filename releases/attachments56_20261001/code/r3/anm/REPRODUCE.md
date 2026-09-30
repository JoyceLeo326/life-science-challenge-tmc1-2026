# 实际复现说明

## 环境

本次实际为Python3.12.10、NumPy2.5.3、SciPy1.18.1、Matplotlib3.11.2、threadpoolctl3.6.0、Windows CPU。SciPy ANM源码公开，执行时限制每个BLAS库为2线程。ProDy整包未安装；已审阅的v2.6.1参考Hessian函数与许可在reference_code。

可在独立Python环境安装`requirements.txt`。该文件是兼容最低版本范围；当前实测版本见ACTUAL_RUN_RECEIPT。跨平台低模的符号/近简并基可变，不能要求所有二进制和时间回执逐字节相同；应核数值精度、子空间和几何指标。

## 外部输入

复制`INPUT_PATHS.template.json`为自己的配置。H_AF_domain和H_AF_full_chainA填写同一个已配准人源PDB；M_PUB_ALL_domain与M_PUB_ALL_full_chainA填写同一个已配准原作者小鼠PDB。程序强制验证与METHOD_FREEZE匹配的SHA；不能换原始未配准文件并仍声称是本次同版输入。

固定作者源模型和恢复方法参见此前代码附件中的reference_restore；源SHA/区域见本包SOURCE_INDEX.csv。外部坐标不在本公开包内。prody_reference_source由安全入口自动定位，用户无需另行配置。

```powershell
python reproduce_anm.py --inputs INPUT_PATHS.local.json --out ../ANM_replay --check-only
python reproduce_anm.py --inputs INPUT_PATHS.local.json --out ../ANM_replay
python plot_anm_results.py --release ../ANM_replay/release
```

输出目录必须独立、新建或为空。入口会核冻结协议、科学脚本、参考源码及两源PDB，再启动真实ANM；原包和源坐标均不改写。新输出的内部路径配置/日志可能含本机路径，交付时应另做匿名副本。

## 新输出

- release：模态NPZ、全部特征值、节点索引、MSF、边界/截断敏感性、24帧位移及几何表、实际回执。
- internal：原参考及24帧Cα PDB，未作全原子重建，不可直接提交为可对接受体。

完整单链对照只有参考结构，不生成展示帧。冻结科学脚本的通用参考文件名含static_domain，实际区域以协议、节点表和作用域ID为准。

公开包校验：

```powershell
python verify_manifest.py
```

本次实际安全入口重新运行退出0；科学表和模态数值与首轮一致，见REPLAY_ENTRY_TEST.json。图件入口也实际生成四图并完成视觉检查。独立检查没有新增对接或MD。
