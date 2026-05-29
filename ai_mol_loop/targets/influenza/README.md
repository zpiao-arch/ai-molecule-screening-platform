# Influenza Target Catalog

This catalog supports first-stage target selection for the AI molecule design
closed loop.

The first MVP should use:

```text
influenza_a_h1n1_na
```

This is the most practical target because neuraminidase has clinically validated
inhibitors, public co-crystal structures, a clear active site, and a direct
biochemical validation route.

Catalog files:

- `target_catalog.json`: target scoring, rationale, structures, guidance.
- `known_drugs.csv`: verified drugs and workflow roles.
- `pdb_structures.csv`: representative structure entries.
- `assay_plan.md`: computational and experimental validation plan.

Recommended target ranking for first product demo:

1. `influenza_a_h1n1_na`
2. `influenza_b_na`
3. `influenza_a_pa_endonuclease`
4. `influenza_b_pa_endonuclease`
5. `influenza_a_m2` as historical control only

