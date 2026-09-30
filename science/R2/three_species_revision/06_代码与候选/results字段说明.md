# results字段说明

- compound_id：候选母体编号。
- conditions_confirmed_n/conditions_total_n：通过条件数/总条件数。
- score_min_kcal_mol、score_max_kcal_mol：静态Vina代理分数范围，单位kcal/mol。
- middle_contacts_n/middle_contacts_total_n：孔道中区4 Å接触次数/总次数。
- min_box_margin_A：盒边界裕度，单位Å；负值表示需谨慎复核。
- qc_status：独立质量控制状态。
- interpretation：证据边界说明。

这些字段不能解释为TMC1实测结合、激活、抑制、选择性或治疗效果。
