# Influenza Target Validation Plan

This file defines the first-stage validation logic for influenza targets. It is
intended for computational workflow planning and does not claim therapeutic
efficacy.

## Recommended MVP Target

Use `influenza_a_h1n1_na` as the first complete closed-loop target.

Reasons:

- It has clinically validated neuraminidase inhibitor controls.
- The pocket is an enzyme active site with public co-crystal structures.
- Docking setup is simpler than membrane channels and metal-dependent PA sites.
- The direct biochemical assay is clear: neuraminidase inhibition.

## Computational Validation

For A/H1N1 neuraminidase:

1. Prepare a receptor from a co-crystal structure such as `3TI6`.
2. Use oseltamivir, zanamivir, and peramivir/laninamivir as positive controls.
3. Generate candidates from the target brief.
4. Run RDKit validity and property filters.
5. Dock candidates and positive controls into the same active site.
6. Run pose plausibility checks on top poses.
7. Rank by multi-objective score:
   - docking or rescoring result
   - pose plausibility
   - validity
   - drug-like properties
   - synthetic accessibility proxy
   - novelty and scaffold diversity
8. Benchmark whether known controls rank sensibly and whether generated top hits
   improve over rounds.

## Experimental Validation Path

The primary experimental assay for neuraminidase targets is a neuraminidase
inhibition assay. Cell-based influenza replication assays can be used as
downstream confirmation.

For PA endonuclease targets, use an endonuclease inhibition assay or cell-based
susceptibility workflow. Metal coordination and resistance substitutions must be
handled carefully.

For M2, do not use it as a first MVP target because adamantane resistance is
widespread and the system is a membrane channel rather than a simple soluble
enzyme pocket.

## Caveats

- Docking scores are not measured potency.
- Proxy scores are not biological validation.
- M2 should be presented as historical comparator, not as the recommended
  current influenza target.
- PA endonuclease docking needs metal-aware preparation.
- Any real-world strain or resistance statement must be checked against current
  surveillance and assay data.

