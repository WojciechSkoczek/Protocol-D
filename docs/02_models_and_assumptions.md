# Models and assumptions (Protocol D)

This document lists the **crucial modelling choices** that define Protocol D.

## Coordinate and vector conventions
- Galactic coordinates (l,b) are used for global directions.
- A dipole measurement is treated as a 3D vector:
  - **D_obs_vec** = D_obs * u(l,b)
  - **D_kin_vec** = D_kin * u(l_CMB, b_CMB) (direction fixed to the CMB dipole)

## Error model (3D covariance)
Each observation has:
- σ_mag (amplitude uncertainty)
- σ_dir (directional uncertainty, degrees; may include a conservative floor)

We approximate the 3×3 covariance by Monte Carlo sampling:
1) perturb the direction within a cone of width σ_dir
2) perturb amplitude by N(D_obs, σ_mag)
3) compute empirical covariance of the sampled vectors

This is implemented in `src/protocol_d/core.py` (see also legacy harness in `src/protocol_d_run_input_pack_any_v2_bao.py`).

## Baseline hypothesis test
We compare:

### NULL
E_i = D_obs_i − D_kin_i (no additional coherent structure)

### BAO_P1 (BAO-anchored coherent excess)
E_i ≈ β * w(z_i) * ĝ

where w(z) is a BAO-weight:
w_i = (r_BAO / χ(z_mean_i))^p

Parameters: ĝ (2 angles) + β (scale).  
Model comparison uses AICc (small-sample corrected AIC).

## Operator (Banach framing) model: O2
We treat a subset of points (Quaia) as affected by an additional survey operator term:

E_i ≈ β * w(z_i) * ĝ + t_i * A(bmin, magcut) * ŝ

- t_i = 1 for Quaia points, else 0
- A(bmin, magcut) = a0 + a_bmin*(bmin − 30) + a_mag*(magcut − 20.0)

### Shared-operator joints (Runs P–Q)
We fit multiple datasets jointly:
- each dataset has its own (ĝ, β)
- the operator axis ŝ and coefficients (a0, a_bmin, a_mag) can be **shared** across datasets

### Identifiability (Run R)
We inject a synthetic rank-2 component (“clumps”) and sweep its amplitude to find the threshold where a second axis becomes detectable by ΔAICc.

## Robustness (Huber)
Some runs use a Huber penalty in whitened residual space as a robustness check (see Run K notes).
This is included as an optional switch in `configs/defaults.yaml`.

