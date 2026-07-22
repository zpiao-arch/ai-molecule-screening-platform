# External Assets

`ASSET_MANIFEST.json` lists the exact local assets used by the validated run. The GitHub source archive intentionally excludes them.

The companion offline asset package mirrors the paths expected by the source tree. Extract it at the repository root to make the source package runnable. Do not commit that directory to a normal GitHub repository; use a private artifact store or Git LFS only after checking every upstream license.

Required external binaries:

- `smina`
- `obabel`

Set `DOCK_BIN_DIR` to the directory containing both binaries. The receptor registry uses a package-relative path and will resolve `scoring/receptors/3TI6_protein_only_obabel.pdbqt` when the companion assets are merged.
