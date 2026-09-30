# Independent review: reproducible Cα ANM for TMC1

Reviewed: 2026-09-30, Beijing time. Scope: method review only; no new docking, MD, coordinate completion, or worm residue remapping was performed by this reviewer.

## Requirements and source facts

The supplied `stage1.md` explicitly allows CPU ProDy ANM as an added conformational condition. It requires honest NMA labeling and keeps the frozen R2 receptor, box, library and rankings separate. The existing three-species release uses 238 observed, sequence-checked, homologous transmembrane Cα positions for structural alignment. These 238 points are not a complete protein model. Its worm docking gate is blocked by four unresolved mapped positions and additional large local correspondence errors. ANM cannot resolve absent coordinates or repair sequence correspondence.

The available static input provenance is H_AF (AFDB human Q8TDI8 canonical v6 aligned monomer), M_PUB_ALL (author mouse six-protein model), and W_7USW_A (experimental 7USW expanded protomer chain A). Expanded is not an experimentally validated open-channel label. Preserve exact source SHA and atom identifiers in every new computation. Author model derivatives remain internal until source redistribution permission is established.

## Recommended frozen baseline

1. Use standard isotropic Cα ANM, uniform spring constant `gamma=1.0`, interaction cutoff `15.0 Å`. These are the official ProDy ANM defaults. Report it as a harmonic, coarse-grained fluctuation model around the supplied coordinates, without membrane/lipid or solvent mechanics. ANM eigenvalues and inverse-eigenvalue variances have relative units; they are not physical frequencies, calibrated Å² fluctuations, temperatures, times, free energies, or channel-state probabilities. [Official ANM API](https://www.bahargroup.org/prody/manual/reference/dynamics/anm.html), [original ANM paper](https://doi.org/10.1016/S0006-3495(01)76033-X).

2. Select one explicit protein unit and preserve its full resolved Cα context. For the low-cost baseline, use all observed TMC1 chain-A Cα positions per source as a stated monomer analysis; mouse auxiliary chains are then omitted from dynamics and this limitation must be named. A separately named all-protein-complex mouse ANM is scientifically different and must not be mixed as if it had the same degrees of freedom. If a full complex exceeds memory/time, do not infer its motions from a cropped pocket.

3. Prefer building ANM on the full chosen resolved unit, then reporting target/TM/pocket displacements by slicing its modes. Building an independent ANM directly on the 238 matched TM points removes loop and terminal coupling and imposes artificial free boundaries. If a TM-only model is necessary, state the exact node list, missing spans and resulting boundaries. ProDy distinguishes mode slicing, full-model reduction and connectivity trimming. [Official editing API](https://www.bahargroup.org/prody/manual/reference/dynamics/editing.html).

4. Parse only one explicit structural MODEL and one alternative-location policy. Use `(model, chain, residue number, insertion code, atom name)` as identity; assert one Cα per retained residue. Keep original numbering, chain identities, residue names, occupancies and source confidence/B-factor metadata. Do not add missing Cα, connect sequence gaps with invented peptide bonds, renumber to a fictitious continuous chain, or treat a renumbered USalign core export as intact sequence. Spatial ENM springs may bridge different sequence segments if their observed Cα distance is within the cutoff; that is a distance contact, not a covalent-chain reconstruction.

5. Limit numerical BLAS/OpenMP/MKL threads to 2 and call ProDy `calcModes(..., nproc=2)` if the installed version supports it; record actual versions. With `N` Cα nodes, dense float64 Hessian storage is `8*(3*N)^2` bytes, before eigensolver workspace and copies. Use lowest-eigenpair diagonalization or official sparse ANM if a complex is too large. Do not report a ProDy run if only a separately implemented SciPy equivalent ran; record the implementation actually used and independently verify its Hessian formula/eigenpairs.

## Numerical gates before interpreting modes

Compute and save these without inspecting ligand docking scores:

- Node count, sequence gaps, duplicate Cα identities, finite-coordinate check, very short/coincident Cα pairs, per-chain counts and contact-degree distribution.
- Undirected Cα contact graph using the same `15.0 Å` spring cutoff: number and size of connected components, isolated nodes and membership table. Connectivity is necessary but does not guarantee mechanical rigidity.
- Hessian symmetry, PSD check, eigenvalue spectrum including zero modes, and residual norms `||H*u - lambda*u||`. State the floating-point threshold explicitly relative to matrix scale. A rigid connected 3D ANM has six translational/rotational zero modes; more zero modes indicate disconnected pieces or floppy constraints and require explanation. Never silently discard an arbitrary number of zero modes and call the next values validated collective motion. ProDy's ANM solver expects six zero modes. [Official implementation](https://www.bahargroup.org/prody/_modules/prody/dynamics/anm.html).
- Retain the first 20 positive modes for export/relative mobility analysis if feasible. Report relative mobility as `sqrt(sum_k ||u_ik||²/lambda_k)` normalized to a declared reference, and name the exact mode truncation. A plot from 20 modes is low-frequency mobility, not total experimental RMSF.
- Show how much each low mode is localized to terminal/low-confidence regions versus the mapped TM core. AlphaFold confidence in flexible tails can change the apparent lowest modes. Record core squared-amplitude fraction and participation ratio. Do not assert a pore mode because the largest arrows happen to be drawn near the pore.
- An optional cutoff sensitivity comparison at 12 and 15 Å can be a separate, predeclared method check. Compare low-mode subspaces, not signed eigenvectors or single vectors across near-degenerate eigenvalues. Do not alter cutoff to improve candidate scores or force a desired number of modes.

## Conformation export and amplitude definition

For a transparent small deterministic set, export the reference plus ± displacements along the first three validated nonzero modes. Predeclare Cα RMSD amplitudes of `0.25` and `0.50 Å` as conservative numerical sensitivity levels. For a unit-normalized eigenvector `u` over `3N` entries, use `x(a) = x0 ± a*sqrt(N)*u`, and independently verify the measured Cα displacement RMSD. Store both the requested and measured amplitude, maximum displacement, mode index/eigenvalue, sign and node list. Fix eigenvector sign by making its largest absolute entry positive for reproducibility.

These amplitudes are chosen analysis settings, not estimated physiological fluctuation sizes. They must not be called calibrated dynamics or sampling probabilities. Every ± frame must remain in the evidence; do not select signs/amplitudes according to favorable ligand score or use an increasing amplitude search to manufacture hits. This is deterministic ANM mode traversal. ProDy traversal also uses user-specified maximum RMSD; its random `sampleModes` uses Gaussian coefficients weighted by inverse square-root eigenvalues and scales the ensemble to a requested mean RMSD. Consequently even official ANM sampling requires an explicit, nonphysical amplitude qualification. [Official sampling API](https://www.bahargroup.org/prody/manual/reference/dynamics/sampling.html), [official sampling formulas and source](https://www.bahargroup.org/prody/_modules/prody/dynamics/sampling.html).

The reference plus all 12 perturbed frames makes 13 inspectable conformations per analysed unit. This is a small sensitivity set, not 13 independent MD representatives. If execution cost calls for six perturbed frames at only `0.25 Å`, freeze that choice before generating or screening anything.

## Geometry and pocket preservation

At minimum, each Cα export must have unchanged node identities/counts and finite coordinates. Record, relative to the static parent: Cα displacement profile, RMSD before and after rigid fitting, adjacent observed backbone Cα-distance changes, all spring-pair distance changes, mapped pocket-pair distance changes, pocket Cα centroid displacement and frozen-box margins. A negative box margin or absent required pocket residue is an explicit exclusion from docking readiness. Fixed-box tests stay in the original parent coordinate frame; a fitted copy must carry and publish its actual transform.

Set any accept/reject tolerances before ligand evaluation and label them project QC choices, not literature-validated biological criteria. Save continuous measurements and explicit reasons so a future reviewer can change a threshold without rerunning the modes. A small amplitude alone is insufficient: a localized mode can yield large maximum displacement even at low global RMSD.

Cα coordinates are enough for checking modal direction, amplitude, network integrity and gross pocket shape. They are not a chemically valid docking receptor. ProDy `extendModel` transfers each residue's Cα vector to all atoms in that residue; it does not rebuild or relax the peptide chain, optimize side chains, or certify atomic geometry. [Official extension API](https://www.bahargroup.org/prody/manual/reference/dynamics/editing.html). Before using an all-atom derivative for independent screening, additionally verify residue/atom identities, covalent bond lengths/angles, peptide continuity, chirality, clashes (including any retained static auxiliary chains), pocket heavy-atom coverage and distances, protonation/charge preparation and exact search-box coverage. If geometry reconstruction or constrained relaxation is unavailable or fails, retain the frame as an ANM visualization and do not submit it as a docking receptor.

Do not derive a channel pore radius from Cα alone or label a larger geometric gap as channel opening. An atom-based pore analysis would require a defined pore axis, full protein atoms, consistent probe/radius rules and an explicit static parent comparison; even then geometric radius alone does not establish conduction.

## Minimum audit deliverables

- Frozen method JSON with selection/cutoff/gamma/solver/thread limits, mode count, amplitude policy and numerical/QC thresholds; source SHA and actual software versions.
- Node table, network components/degrees, eigenvalues including zero modes, eigenvectors, numerical residuals/orthogonality and relative-mobility table.
- Reference and generated Cα coordinate files with frame-to-parent table and hashes; complete frame QC table, no retrospective favorable-only selection.
- Spectrum, relative-mobility plot with missing spans/confidence/TM/pocket annotations, mode-vector view and frame pocket/geometry comparison plot.
- Status explaining whether each frame is Cα analysis only or passed a separate all-atom receptor gate; no automatic docking expansion.

## Review conclusion

A real standard ANM calculation is feasible at low CPU cost for a specified resolved TMC1 monomer. It can satisfy the dynamics-method branch of the supplied stage-1 plan as a new, clearly bounded analysis. It does not itself establish membrane dynamics, channel opening, complete three-species consensus, or chemically valid ensemble docking receptors. The numerical and geometry gates above make a useful result reviewable even when no dynamic docking is run.
