# Evidence Package: Influenza B polymerase acidic protein cap-dependent endonuclease

- Target ID: `influenza_b_pa_endonuclease`
- Recommendation: `secondary_for_b_flu`
- Reason: Relevant for B-flu and clinically connected to baloxavir, but more complex than NA for first MVP.

## Evidence Readiness

- Evidence score: `0.8303`
- Readiness: `ready_for_mvp_closed_loop`
- Meaning: source readiness for computational closed-loop screening, not measured efficacy.

| component | score |
|---|---:|
| clinical_controls | 1.0 |
| structure_metadata | 0.5 |
| literature_metadata | 1.0 |
| assay_path | 1.0 |
| open_data_entry_points | 0.6667 |
| catalog_prior | 0.7732 |

## Closed-Loop Assets

- Primary structure: `6FS8` / ligand `baloxavir acid` / use `secondary_b_flu_structure`
- Structure metadata: `ok`, resolution `1.8`
- Binding-site source: `co_crystal`
- Reference ligand: `baloxavir acid`
- Docking box strategy: Use ligand centroid and preserve metal handling.
- Positive controls: baloxavir marboxil, baloxavir acid
- Reference controls: none
- Historical controls: none
- Assay path: Endonuclease or cell-based susceptibility assays.

## Known Drugs / Controls

| drug | mechanism | workflow role | source |
|---|---|---|---|
| baloxavir marboxil | PA cap-dependent endonuclease inhibitor | positive_control | [source](https://www.cdc.gov/flu/hcp/antivirals/summary-clinicians.html) |
| baloxavir acid | PA cap-dependent endonuclease inhibitor | positive_control | [source](https://www.rcsb.org/structure/6FS8) |

## Structures

| PDB | status | title | method | resolution | PubMed |
|---|---|---|---|---:|---|
| [6FS8](https://www.rcsb.org/structure/6FS8) | ok | Influenza B/Memphis/13/03 endonuclease with bound inhibitor, baloxavir acid (BXA) | X-RAY DIFFRACTION | 1.8 | 29941893 |

## PubMed Query Results

| PMID | year/date | first author | title |
|---|---|---|---|
| [32975358](https://pubmed.ncbi.nlm.nih.gov/32975358/) | 2021 May | Abed Y | Fitness of influenza A and B viruses with reduced susceptibility to baloxavir: A mini-review. |
| [37317069](https://pubmed.ncbi.nlm.nih.gov/37317069/) | 2023 Apr 22 | Saim-Mamoun A | Viral Fitness of Baloxavir-Resistant Recombinant Influenza B/Victoria- and B/Yamagata-like Viruses Harboring the I38T PA Change, In Vitro, Ex Vivo and in Guinea Pigs. |
| [36145480](https://pubmed.ncbi.nlm.nih.gov/36145480/) | 2022 Sep 15 | Saim-Mamoun A | Generation and Characterization of Drug-Resistant Influenza B Viruses Selected In Vitro with Baloxavir Acid. |
| [32526195](https://pubmed.ncbi.nlm.nih.gov/32526195/) | 2020 Oct | Ison MG | Early treatment with baloxavir marboxil in high-risk adolescent and adult outpatients with uncomplicated influenza (CAPSTONE-2): a randomised, placebo-controlled, phase 3 trial. |
| [35292289](https://pubmed.ncbi.nlm.nih.gov/35292289/) | 2022 Apr | Govorkova EA | Global update on the susceptibilities of human influenza viruses to neuraminidase inhibitors and the cap-dependent endonuclease inhibitor baloxavir, 2018-2020. |
| [31107922](https://pubmed.ncbi.nlm.nih.gov/31107922/) | 2019 | Fukao K | Baloxavir marboxil, a novel cap-dependent endonuclease inhibitor potently suppresses influenza virus replication and represents therapeutic effects in both immunocompetent and immunocompromised mouse models. |
| [30184455](https://pubmed.ncbi.nlm.nih.gov/30184455/) | 2018 Sep 6 | Hayden FG | Baloxavir Marboxil for Uncomplicated Influenza in Adults and Adolescents. |

## Open Data Entry Points

- ChEMBL target search: influenza B PA endonuclease
- PubChem compound control: baloxavir acid

## Official / Primary Sources

- [CDC Influenza Antiviral Medications: Summary for Clinicians](https://www.cdc.gov/flu/hcp/antivirals/summary-clinicians.html) - Confirms currently used influenza antivirals and drug classes.
- [CDC Antiviral Drug Resistance among Influenza Viruses](https://www.cdc.gov/flu/treatment/antiviralresistance.html) - Supports M2 adamantane resistance caveat and resistance-aware interpretation.
- [WHO Antiviral susceptibility of influenza viruses: neuraminidase inhibitor methods](https://www.who.int/teams/global-influenza-programme/laboratory-network/quality-assurance/antiviral-susceptibility-influenza/neuraminidase-inhibitor) - Supports neuraminidase inhibition assay as a validation route.
- [FDA Influenza Antiviral Drugs and Related Information](https://www.fda.gov/drugs/information-drug-class/influenza-flu-antiviral-drugs-and-related-information) - Regulatory source for influenza antiviral products and labels.

## Evidence Counts

- PDB entries: 1
- PubMed articles: 7
- Positive controls: 2
- Open data queries: 2

## Interpretation Boundary

- This is a target-evidence package for computational screening.
- It does not establish therapeutic efficacy.
- Full-text papers are not redistributed; only metadata, links, and local summaries are stored.
