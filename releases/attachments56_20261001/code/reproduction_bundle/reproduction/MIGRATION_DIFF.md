# Migration edits

The seven original scientific scripts match the frozen HANDOFF SHA table. Chemistry rules match the prior package manifest. Design constraints, warning IDs, standby intervals, pilot exclusions and original interrupted attempts match the original preparation constants. Formal raw samples and weights match the original training manifest. Source identity and publication identity are separate fields in `SOURCE_MANIFEST.json`.

| File or group | Publication edit |
|---|---|
| chemistry_rules.py | Original bytes retained |
| rank_library.py | Explicit R1 root; normal environment import instead of a local vendor prefix |
| active_acquire.py | Normal environment import instead of a local vendor prefix |
| dock_batch.py | Explicit R1/work roots, receptor/tool arguments and existing frozen Vina/box/receptor identity guards |
| evaluate_forward.py | Four required physical batch paths and existence checks; timing functions accept corresponding explicit roots |
| run_design.py | Explicit R1 parent root and separate writable rule output root |
| train_smiles_gru.py | Explicit work and parent roots; normal environment imports |
| train_selfies_gru.py | Explicit work/legacy roots; normal imports; frozen rule enumeration identity guard |
| freeze/merge auxiliary scripts | Explicit R2 evidence/preflight roots and writable output arguments |
| Three original ledger guards | Original expected SHA retained; fixed, path-only publication hashes accepted through an independently hashed mapping table |
| Interrupted attempts CSV | Four original rows restored; path root only normalized |
| Formal model files | Weights, compressed samples and vocabulary retained as original bytes |
| Training manifest | Logical evidence paths replace local path fields |
| prepare_portable.py | New preparation entry checks the fixed publication table; copies rather than patching anonymized strings |

`migration_diffs/*.diff` show each edited source. Original local root text in the minus lines is rendered as a descriptive audit label to avoid publishing machine paths. These are readable migration records, not patches intended to reconstruct the local original bytes. Original SHA fields refer to the actual source read before editing, not to these redacted diff lines.

Seeds, fingerprint lengths, architecture, training/sample sizes, thresholds, score/charge processing, QC/identity checks, denominator logic and numerical formulas were not changed. Static validation compares the original and publication methods' numeric AST constants and checks syntax/configuration without importing them. It does not execute research calculations or load model weights.
