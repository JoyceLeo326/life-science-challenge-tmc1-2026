# Frozen RDKit PAINS and BRENK source catalog

Pinned version: RDKit 2025.09.6 / Release_2025_09_6.
Tag object: 0ece02e9254ef2d5eeade2bd40eb13546522dbd3.

rule_catalog.csv contains all PAINS_A (16), PAINS_B (55), PAINS_C (409) and BRENK (105): 585 exact source entries. rule_id and description are the exact official descriptions and can be mapped directly from GetDescription(). The ten optional reactive_subset rows are a named subset of BRENK, not additional SMARTS or an additional independent catalog. PAINS, BRENK and reactive subset overlap must not be counted as disjoint molecules.

Use RDKit 2025.09.6 built-in FilterCatalogParams.FilterCatalogs.PAINS and BRENK. Preserve the library input identity and rule freeze binding separately. C++ FilterCatalogEntry.cpp constructs queries with mergeHs=true. The source max field is an allowed-match limit, not the matcher upper bound; an alert starts at source max + 1 (or 1 for source max 0). All frozen PAINS/BRENK entries here use source max 0, so min_count=1 and the default matcher max_count=UINT_MAX=4294967295.

The reactive_subset descriptions are acid_halide, aldehyde, isocyanate, ketene, Michael_acceptor_1 through Michael_acceptor_5, and Three-membered_heterocycle. The last is broader than epoxides. A structural hit does not establish actual chemical reactivity, toxicity, PAINS behavior, human TMC1 activity, or drug efficacy. Report these as structural alerts, not measured outcomes.

## Attribution and retained permissions

Official repository source: https://github.com/rdkit/rdkit/tree/Release_2025_09_6/Code/GraphMol/FilterCatalog
Pinned source definitions: https://github.com/rdkit/rdkit/tree/0ece02e9254ef2d5eeade2bd40eb13546522dbd3/Code/GraphMol/FilterCatalog
API context only: https://www.rdkit.org/docs/source/rdkit.Chem.rdfiltercatalog.html

The repository license.txt is BSD 3-Clause and is retained unmodified in official_sources/license.txt. FilterCatalog source files additionally retain their original 2015 Novartis copyright and BSD-style redistribution notice. Preserve both notices and disclaimers in redistribution; no contributor endorsement is implied. SMARTS are copied as shipped by RDKit, not obtained independently from a publication supplement under an invented separate license.

PAINS citation as provided by the official README:
Baell JB, Holloway GA. New Substructure Filters for Removal of Pan Assay Interference Compounds (PAINS) from Screening Libraries and for Their Exclusion in Bioassays. J Med Chem 53 (2010) 2719-2740. DOI: 10.1021/jm901137j.

BRENK citation as provided by the official README/Filters.cpp:
Brenk R et al. Lessons Learnt from Assembling Screening Libraries for Drug Discovery for Neglected Diseases. ChemMedChem 3 (2008) 435-444. DOI: 10.1002/cmdc.200700139.

This freeze performed text/source parsing and SHA checks only, with no RDKit import, molecular matching or 9098-library calculation. RULES_FROZEN.json is the pre-diagnostic gate; diagnostics must record a start time strictly after frozen_at_cst and bind its SHA plus rule_catalog.csv SHA.
