# Stage 24 hotfix6 — knowledge rules-v8

This delta was introduced after auditing a real rules-v7 Stage-24 run.

## Fixed contracts

1. **AlphaMissense materialization.** The target lookup could report exact matches and populate `event_candidates.alpha_*`, while no `K8_ALPHAMISSENSE_SUPPORT` research overlay was emitted. Rules-v8 explicitly materializes the research-only reason after an exact match. AlphaMissense remains CC BY-NC-SA 4.0 and must never enter the open/commercial-safe TSV.

2. **BND/TRA coordinate identity.** Reciprocal Manta BND mates can be on different chromosomes. Rules-v7 collapsed their numeric coordinates with `min/max` while retaining only the first chromosome, creating a fake same-chromosome interval. Rules-v8 preserves a pair of breakpoint loci in the event identity and provenance. GIAB QC tests the real breakpoint chromosome/position pairs, including the mate locus and cross-caller `STRUCTURAL` events.

3. **Patient priority vs call reality.** ClinGen dosage/AnnotSV biological evidence does not prove that a structural call exists in this patient. K5/K6 single-caller events are therefore kept as evidence but demoted to patient priority 5 until independently confirmed. `CROSS_CALLER_CONFIRMED` can retain its evidence-derived priority; GIAB caution still demotes it.

## Real-run audit that motivated v8

Rules-v7 completed the full catalog and remote target lookups, but the audit found:

- 55 AlphaMissense exact matches and zero AlphaMissense research overlays;
- 191 normalized TRA/BND event clusters whose two breakpoint loci require paired-locus handling; 176 of the reciprocal Manta pair clusters observed in the source are cross-chromosomal;
- several single-caller K5/K6 events ranked ahead of technically confirmed small-variant evidence, which would waste Stage-25 web-research slots on a question literature search cannot answer: whether the structural call itself is real.

Stage 25 remains gated until a rules-v8 run is reviewed.
