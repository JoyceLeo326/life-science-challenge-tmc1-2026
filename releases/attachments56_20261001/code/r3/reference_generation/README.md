# R3 public scientific evidence release

All files in this directory have been screened for internal absolute filesystem paths and the local account name. Upload this directory as a unit; the internal parent directory is not the public release.

## Results and boundaries

- Generated supplement: 177 remaining, 162 identity-eligible, frozen cached ExtraTrees top10; 1 preparation block, 9 actual dockings, independent technical QC9/9, Middle contact5/9, score threshold0/9, joint0/10 selected (0/9 docked). Original40 remains0/40.
- Reference supplement: 11 distinct identity/form rows, 3 preparation blocks, 8 actual dockings, QC8/8. Six dye-reduction drug parents docked, three pass score+Middle: indinavir(-9.102), posaconazole(-8.710), lapatinib(-8.403). This small panel is not a biological model benchmark.
- Amitraz: -6.982, Middle distance9.662 A; its published dye-loading direction is increase. Denatonium: -6.397, specificity comparator only, no direct TMC1 positive label.
- Preparation blocks: amiloride source-tautomer-dependent conservative source gate; ceforanide unsupported unspecified stereo after pH rule; indinavir sulfate multiple components. No rule loosening.
- Amiloride(CHEMBL945) was not prepared or docked because the source gate flagged a potential unspecified C=N double-bond element in its supplied acylguanidine tautomer representation(Bond_Double0). The gate stopped before protonation enumeration. This technical limitation is not evidence of an unspecified tetrahedral center or an experimentally established stable stereochemical ambiguity. Article salt/batch equivalence remains independently unconfirmed. See R3_amiloride_independent_chemistry_review.md; original blocked records and supplied SMILES were not altered.
- One amitraz PDBQT partial-charge warning is retained in the independent QC table; score QC acceptance uses Vina's charge-independent scoring policy. The chemical-state limitation remains.
- 12/15 published dye-reduction members are mapped to original Fig8_stats rows. This is P3 mouse cochlear-explant AM1-43 uptake, not direct TMC1 electrophysiology, binding, selectivity, or therapeutic evidence.
- Original paper confirms13 pharmacophore training references plus3 added MET references (16 names). An asserted68-member library has no complete verified membership source here; no68 distribution was fabricated. Five ZINC active names are mapped but exact structures remain unresolved. Drug database identity does not establish article salt/batch equivalence.

## Raw evidence and hashes

`raw_runs/top10`, `raw_runs/reference_panel`, and `raw_runs/denatonium` contain ligand SDF/PDBQT, full atomic poses, Vina/export/preparation logs, ledgers and frozen plans. Coordinates, scores and molecule graphs are preserved. Metadata paths were normalized to release-relative references. `release_source_to_public_hashes.csv` records original-source SHA256 separately from public-file SHA256; modified metadata never inherits the original hash. `release_files_manifest.csv` verifies all public files.

Original QC summary ledger hashes refer to the original internal ledger, deliberately preserved as provenance; use the hash mapping for the normalized public ledger. No historical ledger hash is relabeled as a public hash.

`protocol_inputs` holds the box metadata and runner. The original author's receptor model is NOT redistributed. Supply the external `M_PUB_ALL.pdbqt` via mandatory `--receptor-dir`. The driver asserts SHA256 `c024cc9efaf910ac11642cea4a4c701a0a88a335d7deb0d040320516d5a2eda1` before work. Original input/receptor/box/Vina hashes are retained in plans. Frozen QC replay also requires that external receptor; the original17/17 readback used an internal research copy. `reference_sources` includes the original paper XML, supplement PDF, text extraction and XLSX. Fixed author repository commit74274dad76201258adbd76c7058b705170c4b857 was supplied by control; a fresh API tree read failed HTTP403 rate-limit, so this release does not claim a new full-repository membership audit.

## Recalculate

Use Python3.12, the recorded requirements, AutoDock Vina1.2.7 and the installed Meeko mk_export executable. The Vina binary is an external prerequisite, with its expected SHA256 in frozen plans. Run from this directory:

```text
python replay_R3.py --batch top10 --receptor-dir EXTERNAL_RECEPTOR_DIRECTORY --vina PATH_TO_VINA_1_2_7 --export PATH_TO_MK_EXPORT --out NEW_OUTPUT_DIRECTORY
python replay_R3.py --batch reference_panel --receptor-dir EXTERNAL_RECEPTOR_DIRECTORY --vina PATH_TO_VINA_1_2_7 --export PATH_TO_MK_EXPORT --out ANOTHER_NEW_OUTPUT_DIRECTORY
python replay_R3.py --batch denatonium --receptor-dir EXTERNAL_RECEPTOR_DIRECTORY --vina PATH_TO_VINA_1_2_7 --export PATH_TO_MK_EXPORT --out ANOTHER_NEW_OUTPUT_DIRECTORY
```

Add `--limit 1` for a one-row entry check; this is not full-batch evidence. The first top10 row is preparation-blocked, so use `--batch denatonium --limit 1` for an actual one-job docking check. The driver creates a new output directory, uses at most2 concurrent docking jobs and independently checks the resulting poses. It does not overwrite frozen evidence. This public path adaptation has been statically checked; no duplicate docking was performed after the actual runs. Floating point/platform differences may change results; exact seeds and original inputs are supplied. The previously issued archive SHA bb1c5c72f8c8fee9ec4131021a1f9cd902de288a1a22d994287eeb1c92ff7fef is withdrawn because it included the restricted receptor; use the replacement archive receipt.

For independent QC of the existing frozen poses without docking, run from this directory: `python audit_real_jobs_R3.py top10 READBACK_top10 raw_runs EXTERNAL_RECEPTOR_DIRECTORY` (replace top10 with reference_panel or denatonium). Supply `TMC1_QC_OUTPUT_DIR` for audit outputs in another directory. The auditor uses the external receptor only after enforcing the original ledger receptor SHA256.

Sources: https://doi.org/10.1038/s42003-025-07943-x ; https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId=12433 ; https://rupress.org/jgp/article/144/1/55/43386/Conductance-and-block-of-hair-cell
