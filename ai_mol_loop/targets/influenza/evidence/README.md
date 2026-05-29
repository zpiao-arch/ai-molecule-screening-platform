# Influenza Target Evidence Index

This directory stores target-evidence metadata used by the second-stage target-source workflow.

## Target Packages

| target id | target | evidence score | readiness | primary PDB | PDB entries | PubMed articles | controls | fetch warnings |
|---|---|---:|---|---|---:|---:|---:|---:|
| `influenza_a_h1n1_na` | Influenza A(H1N1) neuraminidase | 0.9967 | ready_for_mvp_closed_loop | 3TI6 | 3 | 8 | 3 | 0 |
| `influenza_b_na` | Influenza B neuraminidase | 0.8677 | ready_for_mvp_closed_loop | 4CPL | 1 | 7 | 3 | 0 |
| `influenza_a_pa_endonuclease` | Influenza A polymerase acidic protein cap-dependent endonuclease | 0.8659 | ready_for_mvp_closed_loop | 6FS6 | 1 | 6 | 2 | 0 |
| `influenza_b_pa_endonuclease` | Influenza B polymerase acidic protein cap-dependent endonuclease | 0.8303 | ready_for_mvp_closed_loop | 6FS8 | 1 | 7 | 2 | 0 |
| `influenza_a_m2` | Influenza A M2 proton channel | 0.5686 | secondary_or_needs_manual_review | 2KQT | 1 | 6 | 0 | 0 |

## Source Types

### official_sources

- [CDC Influenza Antiviral Medications: Summary for Clinicians](https://www.cdc.gov/flu/hcp/antivirals/summary-clinicians.html)
- [CDC Antiviral Drug Resistance among Influenza Viruses](https://www.cdc.gov/flu/treatment/antiviralresistance.html)
- [WHO Antiviral susceptibility of influenza viruses: neuraminidase inhibitor methods](https://www.who.int/teams/global-influenza-programme/laboratory-network/quality-assurance/antiviral-susceptibility-influenza/neuraminidase-inhibitor)
- [FDA Influenza Antiviral Drugs and Related Information](https://www.fda.gov/drugs/information-drug-class/influenza-flu-antiviral-drugs-and-related-information)

### open_data_sources

- [RCSB Protein Data Bank](https://www.rcsb.org/)
- [PubMed / NCBI E-utilities](https://pubmed.ncbi.nlm.nih.gov/)
- [ChEMBL](https://www.ebi.ac.uk/chembl/)
- [BindingDB](https://www.bindingdb.org/)
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
- [UniProt](https://www.uniprot.org/)

## Boundary

The evidence packages store metadata, source links, and local summaries. They do not redistribute copyrighted full text.
