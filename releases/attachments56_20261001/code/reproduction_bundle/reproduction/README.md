# R2 reproduction supplement

This supplement restores the actual R2 scientific methods and their required frozen support files. It records the source hashes separately from the hashes of the edited publication copies in `SOURCE_MANIFEST.json` and `SOURCE_MANIFEST.csv`. No ranking, docking, evaluation, model training or sample generation was run during this migration.

Use `reproduction/prepare_portable.py` as the current preparation entry. The older `science/R2/evidence/replay/prepare_portable.py` is a historical artifact: its original exact string patches and original input hashes do not match the anonymized source layout. The new entry checks the seven R1 inputs including `common_box.json`, two R2 library snapshots, restored method/support publication hashes, and an explicitly supplied external receptor. It never treats an anonymous path token as a runnable default.

```powershell
$R1 = Read-Host 'Absolute directory containing the seven frozen R1 inputs'
$R2 = Read-Host 'Absolute R2 evidence directory containing 02_library, 03_compute and 04_design'
$WORK = Read-Host 'Absolute new empty replay work directory'
$RECEPTOR = Read-Host 'Absolute path to original local M_PUB_ALL.pdbqt'
python reproduction/prepare_portable.py --r1 $R1 --r2 $R2 --work $WORK --receptor $RECEPTOR --check-only
python reproduction/prepare_portable.py --r1 $R1 --r2 $R2 --work $WORK --receptor $RECEPTOR
```

`--check-only` reads and hashes inputs and makes no work tree. Preparation copies published methods unchanged into `WORK/portable_method/03_compute` and `WORK/portable_method/04_design`, then creates `portable_paths.json` and `PORTABLE_MANIFEST.json`. The latter describes preparation, not a scientific rerun. Output methods use explicit configured roots; originals, input tables and frozen results remain provenance inputs.

No receptor is redistributed in this supplement. The required original receptor SHA256 is `c024cc9efaf910ac11642cea4a4c701a0a88a335d7deb0d040320516d5a2eda1`. The parent repository's source and rights record explains the upstream reference. Supply the original local receptor; a replacement with different bytes fails the check. Docking also requires explicit `--receptor-dir`, `--vina-bin`, and `--export-bin`. Vina must match the original Windows 1.2.7 binary SHA `e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5`; the Meeko export executable must exist. The original box is checked by its frozen SHA.

The evaluation method requires four physical batch paths: initial original154/resumed124 and active original27/resumed14. It uses these paths to read the frozen plans and ledgers without rewriting merged manifests. The wall clock and standby subtraction formulas, failures consuming budget, per-arm counts, independent QC checks, seeds, descriptors, fingerprint lengths and numerical model settings are preserved.

The restored interrupted attempt table has all four original rows. Only the local `partial_log_path` root was mapped to `03_compute/raw_runs/formal_forward_v2_active_new/logs/`. Original SHA `95156f2f5f100bc817df714e98dc31b7ba0a3c41e5b7c8b4741c3db83dd5aead` and the different publication SHA are both recorded. The preexisting publication artifact had only a header and cannot serve as the four-row evidence table.

The actual 20,000-sample formal SELFIES base weights, vocabulary and raw compressed samples are archived under `science/R2/evidence/04_design/model_artifacts/final_selfies_17645`. Weights and samples match the hashes in the original training manifest; vocabulary has an observed source hash with no prior frozen expected hash. The weights are stored as bytes and were never loaded. `training_manifest.json` retains the numerical record while replacing local path fields with logical evidence references. Every archived file is below 90 MiB.

The auxiliary methods require these existing frozen evidence directories: `03_compute/ranking_formal_v2_16698`, `03_compute/ranking_active_v2_after_ai60`, the four physical batches under `03_compute/raw_runs`, and the independent QC/contact tables described by the original `REPLAY.md`. Reconstruction commands must supply explicit writable output paths. Preflight uses `TMC1_PREFLIGHT_ROOT`, which preparation sets to `WORK/ranking_replay`.

Merge and resume methods retain their original frozen ledger hashes. `frozen_publication.py` also accepts the exact three published ledger hashes listed in `publication_hashes.json`, which has its own fixed SHA guard. Before this mapping was emitted, original frozen source identity, row counts, and every nonpath scientific field were checked. Any subsequent change to those published ledger bytes fails. `FROZEN_DEPENDENCY_AUDIT.json` gives the required ranking/QC source locations and reference hash relationships for the compute recovery line.

See `ENVIRONMENT.md`, `MIGRATION_DIFF.md`, and `STATIC_VALIDATION.json` for the environment record, scope of edits and static validation. Preparing copies does not claim reproduction of scores, byte-identical retraining, or current wet-lab evidence.
