# Recorded research environment

The original R2 HANDOFF/REPLAY records specify RDKit 2025.09.6, scikit-learn 1.7.2, Dimorphite-DL 2.0.2 and Meeko 0.8.0. The formal model manifest records PyTorch 2.5.1+cpu and SELFIES distribution metadata 2.2.0, whose runtime reported version was 2.1.1. The restored source is Python 3.10 or later syntax; the research directory also contains Python 3.12 cache evidence. No Python environment was copied.

NumPy and SciPy are required by the actual core imports. Their original exact versions were not established by the allowed frozen records, so the requirements file leaves them unspecified. This is a partial environment record, not a complete lockfile and not a claim of byte-identical reproduction. R2's local sklearn vendor directory reported joblib 1.5.2 and threadpoolctl 3.6.0; these names/versions were read from distribution directory metadata only.

Install dependencies into a separately managed research environment. The original Windows AutoDock Vina executable is external and its frozen SHA is enforced by `dock_batch.py`. Meeko's `mk_export` command is passed explicitly with `--export-bin`. NumPy/SciPy/RDKit/model versions should be recorded for any new rerun.

Published methods no longer prepend the original `vendor` or `pydeps` locations. They import libraries from the selected Python environment. `portable_config.py` reads explicit absolute roots from environment variables or adjacent `portable_paths.json`. Required variables are `TMC1_R1_ROOT`, `TMC1_R2_ROOT`, `TMC1_WORK_ROOT`, `TMC1_RULE_OUTPUT`, and (for preflight only) `TMC1_PREFLIGHT_ROOT`. The preparation entry writes all five for copied methods.

Core source dependencies are documented by actual imports: chemistry rules/RDKit; rank and acquisition/NumPy/scikit-learn/Dimorphite-DL; docking/RDKit/NumPy/SciPy/Meeko/Dimorphite-DL; evaluation/RDKit and ranking helpers; rule enumeration/RDKit including SA Score; SMILES/Torch/RDKit; SELFIES/Torch/RDKit/SMILES helpers.
