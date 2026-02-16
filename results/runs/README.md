# results/runs — directory naming (Windows-safe)

Some original run folders had very long names (and in several cases a duplicated nesting),
which can exceed Windows path limits during extraction or copying.

This repo uses short, stable names:

- RunG … baseline BAO harness outputs
- RunH … subspace check
- RunI … Quaia inclusion ladder
- RunJ … Quaia systematics scenarios
- RunK … operator / Banach framing (O2)
- RunLq … operator manifold variant
- RunM … param operator variant
- RunN … two-axis operator variant
- RunO … model selection summary
- RunP … shared operator axis within variants
- RunQ … shared operator across variants (all4)
- RunR … rank-2 (“clumps”) identifiability sweep

Auxiliary bundles:
- MC_LoTSS_CatWISE
- MC_outputs
- Cleaned_bundle
- QuaiaG20_bundle
- QuaiaZbins_bundle

Original zip bundles are kept unchanged in `results/bundles/`.


Windows note: the historical cleaned bundle is provided as `results/runs/Cleaned_bundle.zip` (not extracted) to avoid path-length issues.
