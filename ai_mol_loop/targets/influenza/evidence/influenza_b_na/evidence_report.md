# Evidence Package: Influenza B neuraminidase

- Target ID: `influenza_b_na`
- Recommendation: `primary_for_b_flu`
- Reason: Good target for B-flu workflows, but first MVP should start with A/H1N1 NA because the target story and demo structure set are simpler.

## Evidence Readiness

- Evidence score: `0.8677`
- Readiness: `ready_for_mvp_closed_loop`
- Meaning: source readiness for computational closed-loop screening, not measured efficacy.

| component | score |
|---|---:|
| clinical_controls | 1.0 |
| structure_metadata | 0.5 |
| literature_metadata | 1.0 |
| assay_path | 1.0 |
| open_data_entry_points | 1.0 |
| catalog_prior | 0.855 |

## Closed-Loop Assets

- Primary structure: `4CPL` / ligand `NA inhibitor or analog` / use `candidate_structure_verify_before_use`
- Structure metadata: `ok`, resolution `2.0`
- Binding-site source: `co_crystal_or_homologous_na_site`
- Reference ligand: `zanamivir or oseltamivir analog`
- Docking box strategy: Use co-crystallized inhibitor centroid or validated NA active-site alignment.
- Positive controls: oseltamivir, zanamivir, peramivir
- Reference controls: none
- Historical controls: none
- Assay path: Neuraminidase inhibition assay with influenza B NA. | Cell-based follow-up only after biochemical enrichment.

## Known Drugs / Controls

| drug | mechanism | workflow role | source |
|---|---|---|---|
| oseltamivir | neuraminidase inhibitor | positive_control | [source](https://www.cdc.gov/flu/hcp/antivirals/summary-clinicians.html) |
| zanamivir | neuraminidase inhibitor | positive_control | [source](https://www.cdc.gov/flu/hcp/antivirals/summary-clinicians.html) |
| peramivir | neuraminidase inhibitor | positive_control | [source](https://www.cdc.gov/flu/hcp/antivirals/summary-clinicians.html) |

## Structures

| PDB | status | title | method | resolution | PubMed |
|---|---|---|---|---:|---|
| [4CPL](https://www.rcsb.org/structure/4CPL) | ok | Structure of the Neuraminidase from the B/Brisbane/60/2008 virus. | X-RAY DIFFRACTION | 2.0 | 24795482 |

## PubMed Query Results

| PMID | year/date | first author | title |
|---|---|---|---|
| [37567255](https://pubmed.ncbi.nlm.nih.gov/37567255/) | 2023 Sep | Ivashchenko AA | Resistance profiles for the investigational neuraminidase inhibitor AV5080 in influenza A and B viruses. |
| [30299356](https://pubmed.ncbi.nlm.nih.gov/30299356/) | 2018 Dec | Lee N | Neuraminidase inhibitor resistance in influenza: a clinical perspective. |
| [11225320](https://pubmed.ncbi.nlm.nih.gov/11225320/) | 2000 Nov | Kaji M | [Neuraminidase inhibitor]. |
| [9651151](https://pubmed.ncbi.nlm.nih.gov/9651151/) | 1998 Jul 2 | Kim CU | Structure-activity relationship studies of novel carbocyclic influenza neuraminidase inhibitors. |
| [29641358](https://pubmed.ncbi.nlm.nih.gov/29641358/) | 2018 | Allen JD | H3N2 influenza viruses in humans: Viral mechanisms, evolution, and evaluation. |
| [37258672](https://pubmed.ncbi.nlm.nih.gov/37258672/) | 2023 Jun | Momont C | A pan-influenza antibody inhibiting neuraminidase via receptor mimicry. |
| [35784305](https://pubmed.ncbi.nlm.nih.gov/35784305/) | 2022 | Bernard MC | Validation of a Harmonized Enzyme-Linked-Lectin-Assay (ELLA-NI) Based Neuraminidase Inhibition Assay Standard Operating Procedure (SOP) for Quantification of N1 Influenza Antibodies and the Use of a Calibrator to Improve the Reproducibility of the ELLA-NI With Reverse Genetics Viral and Recombinant Neuraminidase Antigens: A FLUCOP Collaborative Study. |

## Open Data Entry Points

- ChEMBL target search: influenza B neuraminidase
- BindingDB query: influenza B neuraminidase inhibitor
- PubChem compound controls: oseltamivir, zanamivir, peramivir

## Official / Primary Sources

- [CDC Influenza Antiviral Medications: Summary for Clinicians](https://www.cdc.gov/flu/hcp/antivirals/summary-clinicians.html) - Confirms currently used influenza antivirals and drug classes.
- [CDC Antiviral Drug Resistance among Influenza Viruses](https://www.cdc.gov/flu/treatment/antiviralresistance.html) - Supports M2 adamantane resistance caveat and resistance-aware interpretation.
- [WHO Antiviral susceptibility of influenza viruses: neuraminidase inhibitor methods](https://www.who.int/teams/global-influenza-programme/laboratory-network/quality-assurance/antiviral-susceptibility-influenza/neuraminidase-inhibitor) - Supports neuraminidase inhibition assay as a validation route.
- [FDA Influenza Antiviral Drugs and Related Information](https://www.fda.gov/drugs/information-drug-class/influenza-flu-antiviral-drugs-and-related-information) - Regulatory source for influenza antiviral products and labels.

## Evidence Counts

- PDB entries: 1
- PubMed articles: 7
- Positive controls: 3
- Open data queries: 3

## Interpretation Boundary

- This is a target-evidence package for computational screening.
- It does not establish therapeutic efficacy.
- Full-text papers are not redistributed; only metadata, links, and local summaries are stored.
