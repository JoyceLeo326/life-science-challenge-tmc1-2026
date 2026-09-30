# Independent actual-output QA

Date: 2026-09-30 (Beijing). The reviewed run is `ANM_TMdomain_v1_20260930`. No frozen protocol, run script, main output or delivered R2/R3 package was edited. The independent checker writes only this `independent_review` directory, imported no functions from `run_anm.py`, and used at most two threads in each actual BLAS library.

## Result

**PASS for numerical calculation and faithful reporting. No atomistic docking eligibility is established.** The 24 frames include deliberately retained distortion failures; these do not invalidate the ANM algebra and must not be relabeled or removed.

- The actual source hashes match the frozen configuration, and the frozen run script and protocol hashes still match `METHOD_FREEZE.json`.
- Source Cα node identities match the published node table exactly. Human domain has 552 nodes (A183–734); mouse domain 555 (A177–731). Full resolved A chains have 760 and 666 nodes respectively. No nonblank insertion codes, duplicated residue numbers or missing Middle15 indices were found in these selected inputs. Thus the number-only downstream maps have no insertion-code ambiguity for this run.
- All six networks were independently rebuilt using a full inter-node distance matrix and vectorized Hessian-block scatter accumulation, rather than the parent's KD-tree/loop implementation. Edge counts, one-component connectivity and degree extrema match the published QC.
- Every saved low20 NPZ eigenvalue/vector set has the expected shape, source node order, finite values, orthogonality and agreement with the saved complete eigenvalue table. The published spectra have six eigenvalues within the frozen absolute `1e-6` zero tolerance and no negative eigenvalues below `-1e-6`.
- All six rebuilt matrices annihilate the independent six-column rigid translation/rotation basis (maximum relative residual around `1e-17`). Saved positive modes satisfy independent eigenpair residual checks: maximum relative residual scaled by `||Hu||+|lambda|*||u||` is `3.4304e-9` for the very soft human full-chain mode; other domain results are around `1e-12`. Orthogonality maximum error is below `2.7e-15`.
- As a separate eigensolver spot check, human-domain 15 Å network's lowest 26 eigenvalues were recalculated using LAPACK `evr`, compared with the original full `evd` calculation. Maximum absolute difference was `8.063e-15`; exactly six were zero under the frozen tolerance. This is one independently repeated spectrum, not a claim of six repeated complete diagonalizations.
- Low20 inverse-eigenvalue weighted MSF and within-network mean-one normalization were independently reproduced for all nodes and all six networks.
- All 24 recorded per-node displacement tables reproduce the published NPZ mode, sign and declared `amplitude*sqrt(N)` scale. Every frame's actual RMS displacement, fitted RMSD, maximum displacement, Middle15 displacement and pair distances, consecutive Cα distances, coarse clash counts and fixed-box margin was independently recomputed. Maximum metric discrepancy was below `1e-11`. All stored Cα PDB coordinates match their full-precision frame after the expected three-decimal PDB rounding (maximum coordinate error ≤0.0005 Å).

## Frozen coarse-distortion flags verified

| Source domain | Imposed Cα RMS displacement | Frames | Coarse flag pass |
|---|---:|---:|---:|
| Human | 0.25 Å | 4 | 4 |
| Human | 0.50 Å | 4 | 2 |
| Human | 1.00 Å | 4 | 0 |
| Mouse | 0.25 Å | 4 | 4 |
| Mouse | 0.50 Å | 4 | 4 |
| Mouse | 1.00 Å | 4 | 0 |

The two failing primary-amplitude human frames are mode 2, both signs. Their maximum consecutive-Cα distance changes are approximately 0.231 and 0.283 Å, exceeding the predeclared 0.2 Å diagnostic threshold. All eight 1 Å frames fail that coarse neighbor-distortion criterion. These are imposed-amplitude sensitivity results, not receptor-validation or biological findings. All 24 retain `eligible_for_full_atom_docking=False`.

## Scientific interpretation needed in the report

1. **The selected lowest modes are localized.** Participation ratios are 0.02694/0.02388 for human modes 1/2 and 0.02910/0.07864 for mouse modes 1/2. The Middle15 squared-amplitude fractions are only 0.001012/0.001213 and 0.001390/0.003217, respectively. Maximum displacement occurs at human residues 497/326 and mouse 494/715. The outputs support describing low-frequency ANM directions, but not calling these modes pocket-dominated opening/closing motions. Source B-column metadata remains source metadata, not a simulated mobility measure.

2. **Global amplitude hides local magnitude.** At the declared 0.5 Å global RMS scale, per-node maximum displacement ranges from 2.604 to 4.519 Å. Middle15 per-node maxima are much smaller, about 0.129–0.263 Å. Thus an attractive small global RMS number does not certify all backbone geometry.

3. **Boundary dependence is material.** The published restricted-full-chain low20 basis ranks are both 20, so the parent's QR-based subspace comparison does not suffer rank-deficient arbitrary completion in this run. Domain-vs-restricted-full RMSIP values are approximately 0.669 (human) and 0.769 (mouse); relative-MSF Spearman correlations are approximately 0.637 and 0.726. These are method sensitivity measurements. They do not establish a unique physiological motion model.

4. **Full human chain contains a very soft nonzero mode.** The minimum positive eigenvalue is 6.161 times the fixed absolute zero tolerance; mouse full chain is 38.906 times, whereas domain networks are 6,289–15,286 times. The independent rigid-basis and saved-eigenpair checks support treating the human value as numerically positive in this run, while its stiffness contrast with the cropped domain must remain visible. Do not infer temperature-calibrated amplitude or functional gating from it.

5. **No all-atom gate passed.** Coordinates contain Cα only. There is no side-chain/backbone reconstruction, constrained atomistic relaxation, auxiliary-subunit/dimer clash audit, membrane environment, protonation/charge preparation or new docking. Numerical ANM QC and the coarse distortion flag cannot substitute for those checks.

## Minor metadata caveat

The internal file `static_domain_CA.pdb` is also used for each `*_full_chainA` baseline, and its generated REMARK says `chain A TM envelope only`. For those two internal files, the actual atoms and frozen selections are full observed chain A (760/666 nodes), not TM-domain-only. The source and node tables are authoritative; public prose or any selected future export should avoid inheriting that misleading generic REMARK. No coordinate or mathematical result is affected, and this reviewer did not alter the files.

## Review evidence

- `audit_actual_outputs.py`: independent read-only checker.
- `ACTUAL_INDEPENDENT_QA_RECEIPT.json`: actual final exit-zero receipt, thread count, scope and counts.
- `INDEPENDENT_SOURCE_QA.csv`: frozen input hashes and node identities.
- `INDEPENDENT_NETWORK_QA.csv`: all six independent matrix/eigenpair checks and the one repeated spectrum.
- `INDEPENDENT_FRAME_QA.csv`: all 24 independently checked frames, flags, PDB rounding and centroid displacements.
- `MODE_LOCALIZATION_DIAGNOSTIC.csv`: participation ratios, Middle15 amplitude fractions and maximum-amplitude residues.

The checker initially needed the existing project dependency path to locate `threadpoolctl`; a subsequent run completed all numerical assertions but its receipt serialization needed conversion of a NumPy Boolean-derived count to a standard integer. Both helper issues were repaired only in the independent checker. The final actual run exited zero in approximately 0.66 s. They did not modify the frozen scientific method or any parent result.
