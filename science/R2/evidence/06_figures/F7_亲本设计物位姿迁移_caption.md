# 图 F7｜Trilaciclib 去甲基设计物的更负对接分数伴随位姿迁移

在相同小鼠六链受体、中区搜索盒和 Vina e8 协议下，比较 Trilaciclib 亲本与去甲基设计物 `DES_TWTKMJWAXOPXSX` 的三个离散随机种子（各两条作业；六条均通过独立技术质控）。A 为种子 20260929 的首姿势叠加：亲本接触四个预设 Middle 残基，设计物不接触 Middle 且距搜索盒边界仅 0.637 Å。B–C 分别给出三个种子的 Vina 分数和距 Middle 最近重原子距离。设计物在种子 20260927/20260929 分数更负，但距离 Middle 为 6.854/6.863 Å，且靠近盒边；种子 20260928 回到 Middle 区域时，分数比同种子亲本高 0.046 kcal/mol。三个种子是离散重复，散点不表示连续时间变化或剂量趋势。该计算不支持稳健的同区域改进，也不证明真实亲和力或功能方向。源数据、首姿势坐标哈希、相机参数和渲染脚本见 `09_结构图/图09_区域迁移数据.csv`、`pose_shift_render_manifest.json`、`render_pose_shift.py` 与 `compose_pose_shift_figure.py`。

