# 分子库证据和离线重建

冻结分母：23475来源、23403来源ID、21356母体、17645理化合格、16698前向池；9098是后续固定表示/立体门限后的对接域。库与R2主榜未改。

正式exact二维SDF：library/filtered_for_screening_exact_stereo_fullprops_2d.sdf；17632条、54370979字节、SHA256 46d1073b15e93aca106c5c8ba03ce5d8d1be3537522323c1148671b9712c691a。原26 CTAB/full-key冲突修复13、隔离13；隔离保持CSV权威身份。旧17619 safe SDF未复制为正式结构。

provenance/raw有24页冻结API原字节：旧四页3475条、新二十页20000条，每页SHA已按原summary断言。获取时间、URL、查询式见两acquisition_manifest。chembl_status.json原响应标ChEMBL37、release_date2026-05-01。数据许可CC BY-SA3.0，保留ChEMBL署名与相同许可；许可原页https://chembl.github.io/chembl-licensing/。DATA_SOURCES.md保留原获取说明，本次没有新下载。

原清洗和过滤代码见source_scripts/build_library.py：RDKit2025.09.6，LargestFragmentChooser preferOrganic/Uncharger/canonical tautomer，isomeric SMILES保留立体而母体计数移除立体。允许C/H/N/O/S/P/F/Cl/Br/I/B、有碳、MW150–650、cLogP−2–7、TPSA≤180、可旋转≤12、|形式电荷|≤2。summary.json记录首失败原因分母；不是每个独立过滤器全量命中率，也不是完整Lipinski/Veber/ADMET评估。新增9098结构警示实算另见diagnostics，本节仅描述原库处理。

科学值保持：all_source_records.csv只改raw_page机器路径元数据；各CSV其余字段逐格相同，其他科学CSV和SDF保持原字节。SOURCE_PUBLIC_HASH_MAP.csv分original_sha256/public_sha256；summary/snapshot/audit中的历史hash保持原来源身份，它们不是本发布总清单。用新公开清单核当前字节；不将来源旧hash冒称匿名副本hash。

离线重建需要原RDKit2025.09.6及依赖，单线程，仅写用户指定全新独立目录：
    python source_scripts/replay_library.py --out NEW_INDEPENDENT_DIRECTORY --check-only
    python source_scripts/replay_library.py --out NEW_INDEPENDENT_DIRECTORY
第二条按原prepare、历史标签身份排除、identity/profile从24页处理，不运行对接。默认禁止写回公开目录。本次只静态解析入口、复制/哈希/元数据核对，没有实际运行RDKit重建；不声称新配置入口的科学重算已通过。重算时新路径/时间会变化，科学身份与分母应比较原表，不要求元数据字节一致。

源脚本是匿名配置副本，算法保留：TMC1_LIBRARY_WORK控制新输出根，TMC1_R1_ROOT指历史输入，TMC1_R2_RAW_CACHE/TMC1_R1_APPROVED_CACHE指两批原页。audit_repaired_sdf.py还需完整外部R2树TMC1_R2_SOURCE_ROOT，其他修复源方法需要已生成的原CTAB链和工作快照，不能直接对冻结发布目录运行。当前exact SDF及既有独立修复审计可直接读取，复制不等于重跑化学审计。
