# Evidence Package: Influenza A polymerase acidic protein cap-dependent endonuclease

- Target ID: `influenza_a_pa_endonuclease`
- Recommendation: `secondary`
- Reason: Clinically validated through baloxavir and attractive as an AI-drug-discovery target, but docking setup is metal-dependent and more complex than NA.

## Evidence Readiness

- Evidence score: `0.8659`
- Readiness: `ready_for_mvp_closed_loop`
- Meaning: source readiness for computational closed-loop screening, not measured efficacy.

| component | score |
|---|---:|
| clinical_controls | 1.0 |
| structure_metadata | 0.5 |
| literature_metadata | 1.0 |
| assay_path | 1.0 |
| open_data_entry_points | 1.0 |
| catalog_prior | 0.8172 |

## Closed-Loop Assets

- Primary structure: `6FS6` / ligand `baloxavir acid` / use `secondary_demo_structure`
- Structure metadata: `ok`, resolution `2.291`
- Binding-site source: `co_crystal`
- Reference ligand: `baloxavir acid`
- Docking box strategy: Use baloxavir acid centroid and preserve catalytic metal treatment during receptor preparation.
- Positive controls: baloxavir marboxil, baloxavir acid
- Reference controls: none
- Historical controls: none
- Assay path: Cap-dependent endonuclease inhibition or cell-based susceptibility assay.

## Known Drugs / Controls

| drug | mechanism | workflow role | source |
|---|---|---|---|
| baloxavir marboxil | PA cap-dependent endonuclease inhibitor | positive_control | [source](https://www.cdc.gov/flu/hcp/antivirals/summary-clinicians.html) |
| baloxavir acid | PA cap-dependent endonuclease inhibitor | positive_control | [source](https://www.rcsb.org/structure/6FS6) |

## Structures

| PDB | status | title | method | resolution | PubMed |
|---|---|---|---|---:|---|
| [6FS6](https://www.rcsb.org/structure/6FS6) | ok | Influenza A/California/04/2009 (pH1N1) endonuclease with bound inhibitor, baloxavir acid (BXA) | X-RAY DIFFRACTION | 2.291 | 29941893 |

## PubMed Query Results

| PMID | year/date | first author | title |
|---|---|---|---|
| [39002800](https://pubmed.ncbi.nlm.nih.gov/39002800/) | 2024 Sep | Chen D | High throughput profiling identified PA-L106R amino acid substitution in A(H1N1)pdm09 influenza virus that confers reduced susceptibility to baloxavir in vitro. |
| [29941893](https://pubmed.ncbi.nlm.nih.gov/29941893/) | 2018 Jun 25 | Omoto S | Characterization of influenza virus variants induced by treatment with the endonuclease inhibitor baloxavir marboxil. |
| [30316915](https://pubmed.ncbi.nlm.nih.gov/30316915/) | 2018 Dec | Noshi T | In vitro characterization of baloxavir acid, a first-in-class cap-dependent endonuclease inhibitor of the influenza virus polymerase PA subunit. |
| [35292289](https://pubmed.ncbi.nlm.nih.gov/35292289/) | 2022 Apr | Govorkova EA | Global update on the susceptibilities of human influenza viruses to neuraminidase inhibitors and the cap-dependent endonuclease inhibitor baloxavir, 2018-2020. |
| [33367751](https://pubmed.ncbi.nlm.nih.gov/33367751/) | 2021 Mar 12 | Ivashchenko AA | Synthesis, inhibitory activity and oral dosing formulation of AV5124, the structural analogue of influenza virus endonuclease inhibitor baloxavir. |
| [32064779](https://pubmed.ncbi.nlm.nih.gov/32064779/) | 2020 Jul | Nakauchi M | Rapid detection of an I38T amino acid substitution in influenza polymerase acidic subunit associated with reduced susceptibility to baloxavir marboxil. |

## Open Data Entry Points

- ChEMBL target search: influenza PA endonuclease
- PubChem compound control: baloxavir acid
- BindingDB query: influenza PA endonuclease baloxavir

## Official / Primary Sources

- [CDC Influenza Antiviral Medications: Summary for Clinicians](https://www.cdc.gov/flu/hcp/antivirals/summary-clinicians.html) - Confirms currently used influenza antivirals and drug classes.
- [CDC Antiviral Drug Resistance among Influenza Viruses](https://www.cdc.gov/flu/treatment/antiviralresistance.html) - Supports M2 adamantane resistance caveat and resistance-aware interpretation.
- [WHO Antiviral susceptibility of influenza viruses: neuraminidase inhibitor methods](https://www.who.int/teams/global-influenza-programme/laboratory-network/quality-assurance/antiviral-susceptibility-influenza/neuraminidase-inhibitor) - Supports neuraminidase inhibition assay as a validation route.
- [FDA Influenza Antiviral Drugs and Related Information](https://www.fda.gov/drugs/information-drug-class/influenza-flu-antiviral-drugs-and-related-information) - Regulatory source for influenza antiviral products and labels.

## Evidence Counts

- PDB entries: 1
- PubMed articles: 6
- Positive controls: 2
- Open data queries: 3

## Interpretation Boundary

- This is a target-evidence package for computational screening.
- It does not establish therapeutic efficacy.
- Full-text papers are not redistributed; only metadata, links, and local summaries are stored.
