# 成药性与合成可及性补充

## 做了什么

6个正式候选加5个探索/对照共11个分子，均已取得SwissADME真实网页输出。首批连接55秒超时后仅保留完整5行，核对结果未再增长后只续交剩余6个；没有重复或并发提交。返回结构11/11与输入完整InChIKey相符。ADMETlab官方接口本次404，未把失败请求写成成功预测。

另以官方ADMET-AI 2.0.1已发布模型在本地CPU单线程实际推理11×41=451个端点值，10个模型文件与版本哈希已登记。两套模型结果分列，不混合平均，也不把理化描述符当作ADMET预测。

5个探索/对照分别为Lumacaftor、其甲基酰胺探索物JMW、Trilaciclib历史对照、LDAG生成探索物，以及LDAG原冻结分子量/电荷匹配对照GEN_YUOYQOBVJWRGGY。Trilaciclib不是已确认的LDAG直接母体。未新增筛选或改变冻结分母。

## 结果如何使用

六正式候选中5个SwissADME胃肠吸收预测为Low，只有CHEMBL1076497为High；全部六个BBB预测为No，多数存在多种CYP抑制预测及较低水溶性。它们并非已经证实成药性优秀的候选，应将这些信号用于后续优先次序和结构优化。BBB预测不能代替血迷路屏障通透性证据，也不能据此声称能进入耳蜗。

SwissADME合成可及性和RDKit Ertl–Schuffenhauer实现均真实计算，两个版本分别保留，不能互相替换；数值小通常表示结构层面估计更易合成，不是路线验证、采购证明或成本报价。不同软件的原子计数、疏水性、转动键和SA实现会不同，原值不做手工对齐。

ADMET-AI有10个超出元数据物理范围的输出：蛋白结合率大于100%，或分布容积/肝细胞清除率小于0。异常值完整保留并单列，不能当作合理绝对值，不裁剪伪装成合格。分类端点为模型输出概率，不是个体真实风险概率；不把某个阈值当作安全批准。

|分子|SwissADME SA|RDKit SA|胃肠吸收|BBB|ADMET-AI hERG输出|ADMET-AI AMES输出|
|---|---:|---:|---|---|---:|---:|
|LIB_YWDIXVZCRMMQQY|3.71|2.424|Low|No|0.907|0.454|
|LIB_MYAFDHCCQBIZLC|3.63|2.450|Low|No|0.949|0.465|
|LIB_ZNDXPRCORYDEAP|3.61|2.677|Low|No|0.964|0.283|
|LIB_FCRHJSNRFAAYQD|3.93|2.543|Low|No|0.965|0.246|
|LIB_DNDVGFIEWOHGAH|3.95|2.665|High|No|0.977|0.255|
|LIB_ZVHBDYFMHCQADD|3.63|2.580|Low|No|0.941|0.370|
|LIB_CHEMBL2103870|3.28|2.756|High|No|0.462|0.263|
|LIB_CHEMBL3894860|4.11|3.526|High|No|0.893|0.601|
|DES_JMWNSKZCTCNWSL|3.42|2.790|High|No|0.756|0.323|
|GEN_LDAGBLZIWYMUFJ|2.95|2.434|High|No|0.749|0.310|
|GEN_YUOYQOBVJWRGGY|1.91|1.601|High|Yes|0.604|0.326|

建议摘要表述：完成候选的多参数成药性预测与合成可及性评估，识别溶解度、代谢相互作用和潜在安全性限制，为后续优化提供依据。不能改写为完成药代实验、证实安全有效或可以治疗患者。

来源与方法：[SwissADME官方说明](https://www.swissadme.ch/faq.php)；[非营利使用条款](https://www.swissadme.ch/disclaimer.php)；[SwissADME原始论文](https://doi.org/10.1038/srep42717)；[ADMET-AI官方模型与MIT许可](https://github.com/swansonk14/admet_ai)；[RDKit SA实现](https://github.com/rdkit/rdkit/blob/master/Contrib/SA_Score/sascorer.py)。SwissADME官方CSV保留原始值，本文仅作学术解释。
