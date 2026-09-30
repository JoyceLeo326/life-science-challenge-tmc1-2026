# Anonymous frozen evidence package

Read `三物种修正与补充对接报告.md`, `FINAL_MACHINE_SUMMARY.json`, the24-cellCSV, and `receptor_and_pose_audit/worm_pose_qc.csv` together. Known failures and all omitted-HET collisions are retained. A completed run is not a QC pass or biological positive.

## Contents and scope

`raw/`:26 actual frozen jobs (12 reusedR2 candidates +2 reusedR3 mouse references +12 newjobs), all ligand inputs, all Vina output modes and path-redacted logs with original/public dual hashes. `protocol/`:unchanged chemical rules and new-condition runners, fixed boxes, official experimental7USX-derived rigid receptor. Human/mouse model coordinates are NOT included; only pose ligand coordinates and plots are included. `mapping_review/` preserves original and anonymized log hashes and documented normalization. Original760 correspondence remains a historical hypothesis.

Source metadata paths were anonymized; scientific fields and atom files are unchanged; Vina logs retain scoring tables while local file paths and text encoding/newlines are anonymized. Original and packaged log hashes are kept separately. `anonymous_transformation_receipts.csv` retains original/public hashes. Original condition freezes' input file hashes refer to original unredacted provenance inputs; public transformed input hashes are in that table and FILES_SHA256.csv. Do not confuse metadata anonymization with a new chemical input.

## Dependencies

Python3.12, RDKit2025.09.6, Meeko0.8.0, Dimorphite-DL2.0.2, NumPy, SciPy; AutoDockVina1.2.7 binary SHA e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5. Use your installed `mk_export` / `mk_prepare_receptor` executables. Supply external authorized receptor files whose hashes match the ledgers: M_PUB_ALL=c024cc9efaf910ac11642cea4a4c701a0a88a335d7deb0d040320516d5a2eda1; H_AF=35b4acccd44a5c2c5c8aa7e756abe26c081018a2c26d0899190751c23392f2c1. Bundled experimental worm receptor=907ef675c3f92bbb3efbec7234cbbbf8a3fb1bd412f4191ce574e238ae9b6651.

## Frozen pose reread (no docking)

`python readback_qc.py --mouse <authorized_M_PUB_ALL.pdbqt> --human <authorized_H_AF.pdbqt> --output <outside_package_qc.csv>`

Reports independent current hash/identity/score/geometry checks and agreement with original adjudication; expected protein geometry result25/26, with worm indinavir failing its unchanged box margin. All worm omitted-context collision flags remain separate. Exact public file reread does not prove biological activity.

## Actual docking replay

`python replay.py --species worm --vina-bin <vina_binary> --export-bin <mk_export> --output <new_output_directory> --limit 1`

Omit `--limit1` to execute11 frozen selections (expectedamilorideblocked,10actual). For two human references use `--species human --receptor-dir <authorized_directory>`. Replay output must be outside frozen package. Same chemistry/pH/preparation seeds ande8/Vina seed,2CPU/job,max2workers. It enforces exact receptor and box hashes. Raw exported SDF andPDBQT contain exchangeableN handling; no silent tautomer correction.

## Experimental receptor rebuild

Obtain [official7USX PDB](https://files.rcsb.org/download/7USX.pdb), check source SHA above. Run `python prepare_worm_receptor.py --source-pdb <7USX.pdb> --prepare-bin <mk_prepare_receptor> --output <new_receptor_directory>`. Strict Meeko missing-residue removal is allowed only beyond fixedbox+8Å; it does not fill missing heavy atoms. Exact frozen PDBQT hash is required. Resolved site126heavyatoms retained; terminal oxygen name swaps outside the region are documented.

No new score threshold, cross-species pharmacology, R2 reranking, MD or NMA. Contracted is a geometric state label, not proof of closed state. Full membrane/competitive displacement modeling is future work.
