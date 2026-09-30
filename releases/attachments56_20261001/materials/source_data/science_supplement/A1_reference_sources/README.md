# 参考来源与补充敏感性证据包

先读 REFERENCE_SOURCE_CLARIFICATION.md，再读 n_geometry_sensitivity/RESULTS.md。

这是 3 个精确 ZINC 数据库结构的来源补回，以及保留原失败后的独立 N 几何敏感性分析。6 个姿势不是 6 个新母体，不并入原 R2/R3 命中率。

原作者 10 个模型的 active SDF 合计 78 重复行、12 个训练分子；与论文 12 个染料活性候选是不同集合。68 个独立参考库的来源仍未证实。

未附原作者整套模型和商业供应商价格。作者结构来源保留文件哈希、冻结提交与逐条链接。内部绝对路径已替换为通用占位符，原分子/姿势字节和科学参数保持。

重新执行使用同版依赖、Vina 1.2.7、原冻结受体与框。设置 TMC1_R1_ROOT、TMC1_WORK_ROOT，向 protocol_inputs/dock_batch.py 显式传入 input、batch、scratch、receptor-dir、vina-bin、export-bin、box、box-meta。sensitivity_inputs.csv 对应敏感性版，resolved_ZINC_reference_inputs.csv 对应 frozen_protocol_failure 原版。不要交叉使用两套脚本。

共同参数：--receptor M_PUB_ALL --seed 20260927 --exhaustiveness 8 --workers 1 --timeout 600。每次 Vina 使用 2 CPU。完整数据用于证据复核，不能作为已验证的生物活性。
