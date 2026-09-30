# 模型说明

正式生成模型是从本轮公开库结构从零训练的无条件SELFIES单层GRU，embedding48、hidden96、51763参数、67词表、10epochs、20000采样、未微调。按Murcko骨架分训练/验证：15952训练结构、1520留出结构；17472进入模型范围，17645是原库输入分母。权重、词表、原采样、训练损失与配置均保存。模型学习化学序列分布，没有使用TMC1对接或生物活性标签作条件化。原40次生成结构尝试联合达标0/40如实保留。

前瞻排序代理：Morgan radius2/1024位+8描述符、ExtraTrees160棵树、max_features0.5、min_samples_leaf2、单线程、seed20260929；对照是相同历史935标签的描述符Ridge alpha100与随机选择。主动阶段只纳入60个AI初段终局独立QC通过分数，与935历史标签合计995；不读取Random/Ridge新分数来训练。候选树间标准差是模型内部分歧，不是生物学误差条或校准置信区间。

所有标签为相同静态受体协议的Vina代理，不能当结合活性、激活/抑制方向、人源效果或疾病补偿。MW、骨架和电荷分布影响臂间比较；MW≥350匹配子集单列。小规模1epoch示例用于运行演示，不取代正式模型，性质合格数0同样保留。

方法来源：RDKit、scikit-learn、AutoDock Vina1.2.7、Meeko0.8.0、Dimorphite-DL2.0.2、PyTorch2.5.1+cpu、SELFIES2.2.0。数据源、处理和许可记录在 `reproduction_bundle/science/R1/provenance`、`SOURCE_RIGHTS`和R2冻结输入/方法清单。第三方上游受体坐标由固定源本地取得，不随公开包分发。
