# Evidence Package: Influenza A(H1N1) neuraminidase

- Target ID: `influenza_a_h1n1_na`
- Recommendation: `primary`
- Reason: Best first closed-loop target: clinically validated, clear enzymatic pocket, public co-crystal structures, known inhibitors, and straightforward biochemical assay.

## Evidence Readiness

- Evidence score: `0.9967`
- Readiness: `ready_for_mvp_closed_loop`
- Meaning: source readiness for computational closed-loop screening, not measured efficacy.

| component | score |
|---|---:|
| clinical_controls | 1.0 |
| structure_metadata | 1.0 |
| literature_metadata | 1.0 |
| assay_path | 1.0 |
| open_data_entry_points | 1.0 |
| catalog_prior | 0.9334 |

## Closed-Loop Assets

- Primary structure: `3TI6` / ligand `oseltamivir` / use `recommended_demo_structure`
- Structure metadata: `ok`, resolution `1.69`
- Binding-site source: `co_crystal`
- Reference ligand: `oseltamivir`
- Docking box strategy: Use co-crystallized oseltamivir centroid from 3TI6 or a prepared receptor-ligand complex.
- Positive controls: oseltamivir, zanamivir, peramivir
- Reference controls: laninamivir
- Historical controls: none
- Assay path: Neuraminidase inhibition assay is the most direct biochemical validation route. | Cell-based influenza replication assays are downstream confirmation, not first-line CLI validation.

## Known Drugs / Controls

| drug | mechanism | workflow role | source |
|---|---|---|---|
| oseltamivir | neuraminidase inhibitor | positive_control | [source](https://www.cdc.gov/flu/hcp/antivirals/summary-clinicians.html) |
| zanamivir | neuraminidase inhibitor | positive_control | [source](https://www.cdc.gov/flu/hcp/antivirals/summary-clinicians.html) |
| peramivir | neuraminidase inhibitor | positive_control | [source](https://www.cdc.gov/flu/hcp/antivirals/summary-clinicians.html) |
| laninamivir | neuraminidase inhibitor | reference_control | [source](https://pdb101.rcsb.org/motm/113) |

## Structures

| PDB | status | title | method | resolution | PubMed |
|---|---|---|---|---:|---|
| [3TI6](https://www.rcsb.org/structure/3TI6) | ok | Crystal structure of 2009 pandemic H1N1 neuraminidase complexed with oseltamivir | X-RAY DIFFRACTION | 1.69 | 22028647 |
| [3TI5](https://www.rcsb.org/structure/3TI5) | ok | Crystal structure of 2009 pandemic H1N1 neuraminidase complexed with Zanamivir | X-RAY DIFFRACTION | 1.9 | 22028647 |
| [3TI3](https://www.rcsb.org/structure/3TI3) | ok | Crystal structure of 2009 pandemic H1N1 neuraminidase complexed with laninamivir | X-RAY DIFFRACTION | 1.8 | 22028647 |

## PubMed Query Results

| PMID | year/date | first author | title |
|---|---|---|---|
| [23097226](https://pubmed.ncbi.nlm.nih.gov/23097226/) | 2013 Jan | Tolentino-Lopez L | Outside-binding site mutations modify the active site's shapes in neuraminidase from influenza A H1N1. |
| [19917319](https://pubmed.ncbi.nlm.nih.gov/19917319/) | 2010 Feb | Okomo-Adhiambo M | Host cell selection of influenza neuraminidase variants: implications for drug resistance monitoring in A(H1N1) viruses. |
| [25605596](https://pubmed.ncbi.nlm.nih.gov/25605596/) | 2015 Jan | Gema LR | Targeting a cluster of arginine residues of neuraminidase to avoid oseltamivir resistance in influenza A (H1N1): a theoretical study. |
| [21592407](https://pubmed.ncbi.nlm.nih.gov/21592407/) | 2011 May 19 | Thorlund K | Systematic review of influenza resistance to the neuraminidase inhibitors. |
| [25301400](https://pubmed.ncbi.nlm.nih.gov/25301400/) | 2014 | Zaraket H | Characterization of human Influenza Viruses in Lebanon during 2010-2011 and 2011-2012 post-pandemic seasons. |
| [35304163](https://pubmed.ncbi.nlm.nih.gov/35304163/) | 2022 Apr | Brown SK | Characterization of influenza B viruses with reduced susceptibility to influenza neuraminidase inhibitors. |
| [23575174](https://pubmed.ncbi.nlm.nih.gov/23575174/) | 2013 Sep | Okomo-Adhiambo M | Neuraminidase inhibitor susceptibility surveillance of influenza viruses circulating worldwide during the 2011 Southern Hemisphere season. |
| [26808479](https://pubmed.ncbi.nlm.nih.gov/26808479/) | 2016 Apr | Okomo-Adhiambo M | Standardizing the influenza neuraminidase inhibition assay among United States public health laboratories conducting virological surveillance. |

## Open Data Entry Points

- ChEMBL target search: influenza A neuraminidase
- BindingDB query: influenza neuraminidase oseltamivir
- PubChem compound controls: oseltamivir, zanamivir, peramivir, laninamivir

## Official / Primary Sources

- [CDC Influenza Antiviral Medications: Summary for Clinicians](https://www.cdc.gov/flu/hcp/antivirals/summary-clinicians.html) - Confirms currently used influenza antivirals and drug classes.
- [CDC Antiviral Drug Resistance among Influenza Viruses](https://www.cdc.gov/flu/treatment/antiviralresistance.html) - Supports M2 adamantane resistance caveat and resistance-aware interpretation.
- [WHO Antiviral susceptibility of influenza viruses: neuraminidase inhibitor methods](https://www.who.int/teams/global-influenza-programme/laboratory-network/quality-assurance/antiviral-susceptibility-influenza/neuraminidase-inhibitor) - Supports neuraminidase inhibition assay as a validation route.
- [FDA Influenza Antiviral Drugs and Related Information](https://www.fda.gov/drugs/information-drug-class/influenza-flu-antiviral-drugs-and-related-information) - Regulatory source for influenza antiviral products and labels.

## Evidence Counts

- PDB entries: 3
- PubMed articles: 8
- Positive controls: 3
- Open data queries: 3

## Interpretation Boundary

- This is a target-evidence package for computational screening.
- It does not establish therapeutic efficacy.
- Full-text papers are not redistributed; only metadata, links, and local summaries are stored.
