# Evidence Package: Influenza A M2 proton channel

- Target ID: `influenza_a_m2`
- Recommendation: `historical_control_not_primary`
- Reason: Useful as a historical control, but adamantane resistance is widespread and CDC guidance does not recommend amantadine/rimantadine for currently circulating influenza A viruses.

## Evidence Readiness

- Evidence score: `0.5686`
- Readiness: `secondary_or_needs_manual_review`
- Meaning: source readiness for computational closed-loop screening, not measured efficacy.

| component | score |
|---|---:|
| clinical_controls | 0.0 |
| structure_metadata | 0.5 |
| literature_metadata | 1.0 |
| assay_path | 1.0 |
| open_data_entry_points | 0.6667 |
| catalog_prior | 0.5384 |

## Closed-Loop Assets

- Primary structure: `2KQT` / ligand `amantadine` / use `historical_control_structure`
- Structure metadata: `ok`, resolution ``
- Binding-site source: `channel_structure`
- Reference ligand: `amantadine`
- Docking box strategy: Requires membrane/channel-aware preparation; not recommended for first docking MVP.
- Positive controls: none
- Reference controls: none
- Historical controls: amantadine, rimantadine
- Assay path: Channel inhibition and resistance-specific viral assays would be needed.

## Known Drugs / Controls

| drug | mechanism | workflow role | source |
|---|---|---|---|
| amantadine | M2 proton channel blocker | historical_control_not_recommended | [source](https://www.cdc.gov/flu/treatment/antiviralresistance.html) |
| rimantadine | M2 proton channel blocker | historical_control_not_recommended | [source](https://www.cdc.gov/flu/treatment/antiviralresistance.html) |

## Structures

| PDB | status | title | method | resolution | PubMed |
|---|---|---|---|---:|---|
| [2KQT](https://www.rcsb.org/structure/2KQT) | ok | Solid-state NMR structure of the M2 transmembrane peptide of the influenza A virus in DMPC lipid bilayers bound to deuterated amantadine | SOLID-STATE NMR |  | 20130653 |

## PubMed Query Results

| PMID | year/date | first author | title |
|---|---|---|---|
| [24011996](https://pubmed.ncbi.nlm.nih.gov/24011996/) | 2013 Oct | Gu RX | Structural and energetic analysis of drug inhibition of the influenza A M2 proton channel. |
| [37377066](https://pubmed.ncbi.nlm.nih.gov/37377066/) | 2023 Aug 15 | Stampolaki Μ | A Study of the Activity of Adamantyl Amines against Mutant Influenza A M2 Channels Identified a Polycyclic Cage Amine Triple Blocker, Explored by Molecular Dynamics Simulations and Solid-State NMR. |
| [25387967](https://pubmed.ncbi.nlm.nih.gov/25387967/) | 2015 | Gu R | Drug inhibition and proton conduction mechanisms of the influenza a M2 proton channel. |
| [16493107](https://pubmed.ncbi.nlm.nih.gov/16493107/) | 2006 Feb 22 | Weinstock DM | Adamantane resistance in influenza A. |
| [31894969](https://pubmed.ncbi.nlm.nih.gov/31894969/) | 2020 Feb 4 | Thomaston JL | X-ray Crystal Structures of the Influenza M2 Proton Channel Drug-Resistant V27A Mutant Bound to a Spiro-Adamantyl Amine Inhibitor Reveal the Mechanism of Adamantane Resistance. |
| [18258311](https://pubmed.ncbi.nlm.nih.gov/18258311/) | 2008 Aug | Pabbaraju K | Adamantane resistance in circulating human influenza A viruses from Alberta, Canada (1970-2007). |

## Open Data Entry Points

- ChEMBL target search: influenza A M2 proton channel
- PubChem compound controls: amantadine, rimantadine

## Official / Primary Sources

- [CDC Influenza Antiviral Medications: Summary for Clinicians](https://www.cdc.gov/flu/hcp/antivirals/summary-clinicians.html) - Confirms currently used influenza antivirals and drug classes.
- [CDC Antiviral Drug Resistance among Influenza Viruses](https://www.cdc.gov/flu/treatment/antiviralresistance.html) - Supports M2 adamantane resistance caveat and resistance-aware interpretation.
- [WHO Antiviral susceptibility of influenza viruses: neuraminidase inhibitor methods](https://www.who.int/teams/global-influenza-programme/laboratory-network/quality-assurance/antiviral-susceptibility-influenza/neuraminidase-inhibitor) - Supports neuraminidase inhibition assay as a validation route.
- [FDA Influenza Antiviral Drugs and Related Information](https://www.fda.gov/drugs/information-drug-class/influenza-flu-antiviral-drugs-and-related-information) - Regulatory source for influenza antiviral products and labels.

## Evidence Counts

- PDB entries: 1
- PubMed articles: 6
- Positive controls: 0
- Open data queries: 2

## Interpretation Boundary

- This is a target-evidence package for computational screening.
- It does not establish therapeutic efficacy.
- Full-text papers are not redistributed; only metadata, links, and local summaries are stored.
