# F3 规则设计与亲本同协议配对

04线从 Lumacaftor 与 Trilaciclib 两个商业母体出发的7个受约束规则结构，和2个亲本按同一小鼠六链受体、预设中区盒、Vina e8、seed 20260929 实际复算。横轴为对接分数 kcal/mol，数值 ΔVina=设计物−亲本；蓝色为与预定15中区残基有≤4 Å接触的设计物，橙色为未满足该接触条件的区域迁移姿势。7对中5对分数更负，4对兼具中区接触与更负分数；这些都是一次静态计算排序，不能证明结合、功能作用或真实优化。DES_TWT 最负−9.895却无中区接触，三种子复核不支持稳健同区改良。源表和代码：04线 `actual_parent_product_comparison.csv`、`analyze_parent_product.py`；输入哈希与图注见 `actual_parent_product_comparison_manifest.json`。05线已做逐作业独立QC和中区接触复核。本图取自04线已生成版本，06线目视检查，未重复绘制。
