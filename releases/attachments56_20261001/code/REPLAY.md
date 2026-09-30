# 科学方法真实重放

所有命令从本包根目录执行。六项R1输入及盒元数据、R2库和科研方法均在 `reproduction_bundle`；受体和Windows Vina可执行文件显式从外部供应。上游原作者受体未获得明确再分发许可，本包不包含其坐标。按 `reference_restore/README.md` 从固定源取得并本地恢复；既有本地研究文件仅作本次试运行。

```powershell
$BUNDLE = (Resolve-Path reproduction_bundle).Path
$R1 = "$BUNDLE\science\R1"
$R2 = "$BUNDLE\science\R2\evidence"
$WORK = Read-Host '新的独立可写重放目录'
$RECEPTOR = Read-Host '外部M_PUB_ALL.pdbqt文件绝对路径'
python reproduction_bundle/reproduction/prepare_portable.py --r1 $R1 --r2 $R2 --work $WORK --receptor $RECEPTOR --check-only
python reproduction_bundle/reproduction/prepare_portable.py --r1 $R1 --r2 $R2 --work $WORK --receptor $RECEPTOR
$P03 = "$WORK\portable_method\03_compute"
$P04 = "$WORK\portable_method\04_design"
python "$P03\rank_library.py" --library "$R2\02_library\formal_snapshot\forward_unlabeled_pool.csv" --output "$WORK\ranking_replay" --per-arm 120 --ai-initial 60 --exclude-ids "$P03\runtime_pilot_12.csv" --expected-library-sha256 8fd8424797e441a3921ba3805dd8e9bf87d3421e5fe3440a32211753ec964c1e
python "$P03\active_acquire.py" --freeze "$R2\03_compute\ranking_formal_v2_16698" --ai-attempts "$R2\05_independent_audit\AI_initial_60_independent_outcomes.csv" --output "$WORK\active_replay"
python "$P03\prepare_active_physical.py" --freeze "$R2\03_compute\ranking_formal_v2_16698" --active "$R2\03_compute\ranking_active_v2_after_ai60" --output "$WORK\active_physical"
python "$P03\merge_formal_batches.py" --output "$WORK\merged_initial" --require-complete
python "$P03\freeze_active_resume.py" --output "$WORK\active_resume"
python "$P03\merge_active_batches.py" --output "$WORK\merged_active"
```

预检只对新的完整排名输出运行；它写新排名目录的预检JSON：

```powershell
$env:TMC1_PREFLIGHT_ROOT = "$WORK\ranking_replay"
python "$P03\freeze_preflight_v2.py"
```

完整评价保留原物理批次四个显式路径，以便重算暂停/待机时间；公开账中的历史路径字段为来源元数据，不作为实际文件定位默认值。

```powershell
$RAW = "$R2\03_compute\raw_runs"
python "$P03\evaluate_forward.py" --freeze "$R2\03_compute\ranking_formal_v2_16698" --active "$R2\03_compute\ranking_active_v2_after_ai60" --initial-old-batch "$RAW\formal_forward_v2_initial278" --initial-resumed-batch "$RAW\formal_forward_v2_resume124" --active-old-batch "$RAW\formal_forward_v2_active_new" --active-resumed-batch "$RAW\formal_forward_v2_active_resume14" --initial-batch "$RAW\formal_forward_v2_merged" --active-batch "$RAW\formal_forward_v2_active_merged41" --initial-qc "$R2\05_independent_audit\formal_v2_initial_independent_qc.csv" --active-qc "$R2\05_independent_audit\formal_v2_active_independent_qc.csv" --initial-contacts "$R2\03_compute\method\formal_v2_initial_middle_contacts.csv" --active-contacts "$R2\03_compute\method\formal_v2_active_middle_contacts.csv" --output "$WORK\evaluation_replay"
```

接触输入从R2证据根读取；准备入口只复制冻结支持清单中的方法文件。

代表性重新对接入口：

```powershell
$VINA = Read-Host 'Windows Vina1.2.7二进制路径'
$EXPORT = Read-Host 'Meeko0.8.0 mk_export可执行文件路径'
$RECEPTOR_DIR = Split-Path $RECEPTOR
python "$P03\dock_batch.py" --input "$R2\03_compute\ranking_formal_v2_16698\physical_union_docking_list.csv" --batch representative1 --scratch "$WORK\docking" --limit 1 --receptor M_PUB_ALL --receptor-dir $RECEPTOR_DIR --vina-bin $VINA --export-bin $EXPORT --box "$R1\coarse_v1\config\vina_box.txt" --box-meta "$R1\coarse_v1\config\common_box.json" --seed 20260927 --exhaustiveness 8 --workers 2 --timeout 600
```

受体SHA `c024cc9efaf910ac11642cea4a4c701a0a88a335d7deb0d040320516d5a2eda1`，VinaSHA `e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5`，盒和元数据也有冻结SHA。更换受体、模型或条件属于新研究，不混入原固定预算结果。

`run_example.py --train-small`为实际可运行的小例。正式生成若确需重训，应按原 `training_manifest.json` 使用完整17645库、10epochs、20000样本与0微调，不把小例的1epoch当正式模型。新训练输出的耗时、文件路径和时间戳会变化；正式已保存权重有独立SHA，模型和样本的旧权威身份由该SHA确定。

原始中断的四个非终局尝试保留未知资源成本；2×清醒墙钟属于请求槽位容量估计，未计作实测CPU秒。结果仅支持本次条件下的计算排序与几何判据。
