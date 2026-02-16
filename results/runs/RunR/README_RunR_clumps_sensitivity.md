# Protocol D — Run R (clumps / "map distortion" sensitivity)

## Purpose
You suggested that the remaining *directional scatter* might be driven by unrecognised matter clumps (a “Mercator-map” problem),
rather than a clean question of isotropy vs anisotropy.

Run R tests a concrete version of that idea in our operator framework:

- **Run M** (1-axis operator):  E ≈ β w(z) ĝ + t * A(x) ŝ
- **Run N** (2-axis operator):  E ≈ β w(z) ĝ + t * [A1(x) ŝ1 + A2(x) ŝ2]

where t=1 only for Quaia points, so the operator cannot “cheat” by fitting the non-Quaia tracers.

If there is a *genuine second independent drift direction* (rank-2 “clump field” in selection space),
Run N should eventually beat Run M once that component is strong enough.

## Injection experiment
For each (pack, variant) dataset separately:

1) Fit **Run M** and **Run N** to the real data (λ=0).  
   Baseline result reproduces earlier finding: **Run N is disfavoured by ~+24 AICc**.

2) Construct a synthetic “clump” component:
   - Choose a direction ŝ2,true orthogonal to the fitted ŝ (Run M).
   - Choose an amplitude function f(x)=b·x where b is **orthogonal to the fitted Run M coefficient vector a**
     (so the new component is not just a rescaling of the existing operator).
   - Standardise f(x) over Quaia points to mean 0, stdev 1.

3) Inject into Quaia only:
   ΔE_i = λ * f(x_i) * ŝ2,true

4) Refit Run M and Run N at each λ and record:
   ΔAICc(λ) = AICc(RunN) − AICc(RunM)

## Key readout
- If **ΔAICc(λ) > 0**, data prefer **one-axis** (Run M).
- If **ΔAICc(λ) < 0**, data prefer **two-axis** (Run N).

We estimate λ* where ΔAICc crosses 0 (linear interpolation between grid points).

## Results (summary)
See `RunR_thresholds.csv`.

Typical crossing amplitudes:
- boehme / baseline:   λ* ≈ 0.0057
- boehme / 2MRS_clean: λ* ≈ 0.0048
- wagenveld / baseline:   λ* ≈ 0.0057
- wagenveld / 2MRS_clean: λ* ≈ 0.0047

For scale: the median |E| in our datasets is ~0.017 (dipole-amplitude units),
so λ* ~ 0.005 corresponds to ~30% of a “typical” excess-vector norm.

Interpretation (within this toy injection):
- If an unmodelled clump-driven component introduces a *second independent drift direction* in selection space
  at the level of **λ ≳ 0.005**, then **Run N would become AICc-preferred**.
- Because the real data prefer Run M (ΔAICc ~ +24 at λ=0), any such rank-2 component must be **below that rough scale**,
  or structured in a way that projects mostly onto the existing one-axis operator.

## Plots
- `RunR_dAICc_vs_lambda_*` (per dataset and combined)
- `RunR_axis2_recovery_*` shows that as λ increases, Run N recovers ŝ2,true (up to sign).

## Caveats
This is a *sensitivity experiment*, not a physical forward model of LSS:
- the injected f(x) is linear in our feature vector x, chosen orthogonal to the fitted one-axis amplitude;
- λ is in the same units as our dipole vectors.

Still, it answers a useful identifiability question:
**how large a second independent “map distortion / clump” component would have to be before the data demand rank-2 structure.**
