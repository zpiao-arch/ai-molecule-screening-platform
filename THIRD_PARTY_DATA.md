# Third-Party Data Boundary

This repository is source-available under the root `LICENSE`; it is not an OSI-licensed open-source release.

The public GitHub archive intentionally excludes machine-readable row-level research data from `legacy/` (`.csv`, `.json`, and `.parquet`). The local working copy may contain those files for historical audit purposes, but they are not redistributed by `scripts/build_source_release.py` until each file has a recorded upstream license and redistribution decision.

| Local material | Historical source described by the scripts | Public archive status |
| --- | --- | --- |
| `legacy/multitarget_benchmark/build_candidates_and_panel.py` | DrugCentral 2021 structures and ChEMBL structures | Script retained; generated rows excluded |
| `legacy/multitarget_benchmark/build_bigrun_library.py` | ChEMBL 37 representations plus aligned labels | Script retained; generated rows excluded |
| `scoring/data/` and `examples/` | Small hand-curated CLI examples | Included as example inputs |
| Offline model/data archive | BindingDB, ChEMBL, Uni-Mol, ADMET and receptor assets | Separate local archive; verify upstream terms before sharing |

The builder retains historical Python, shell, Markdown and SVG files because they are source or aggregate documentation, not row-level data packages. A retained script can still contain old absolute-path conventions; it is not a supported production entry point.

Before any public redistribution of excluded data, record the exact upstream version, URL, license/terms, whether redistribution is permitted, and whether the file is a derived dataset. Do not infer permission from the fact that a source can be downloaded.

