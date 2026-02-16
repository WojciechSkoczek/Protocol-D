#!/usr/bin/env python3
"""
Protocol D — run harness using an external input pack (CSV/JSON)

This is a *harness*:
- It compares a NULL model (no coherent excess) to a rank‑1 coherent excess model:
      D_obs_i = D_kin_i + A_i * g_hat
- It uses an approximate directional covariance derived from (sigma_mag, sigma_dir_deg_approx).
- It is NOT a replacement for re‑measuring dipoles from the original catalogues with a consistent estimator.

Usage:
  python protocol_d_run_input_pack_wagenveld2025_tier1.py --input protocol_d_input_pack_wagenveld2025_tier1.csv --variant baseline
  python protocol_d_run_input_pack_wagenveld2025_tier1.py --input protocol_d_input_pack_wagenveld2025_tier1.csv --variant 2MRS_clean

Outputs:
  - protocol_d_run_results_<variant>.txt
"""

from __future__ import annotations
import argparse
import json
import math
import warnings
from dataclasses import dataclass
from typing import Optional, List

import numpy as np

try:
    from scipy.optimize import minimize
except Exception as e:
    raise SystemExit("scipy is required (scipy.optimize.minimize). Install scipy and rerun.") from e

# ---------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------
def galactic_to_unit(l_deg: float, b_deg: float) -> np.ndarray:
    l = math.radians(float(l_deg))
    b = math.radians(float(b_deg))
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)], dtype=float)

def galactic_to_cartesian(l_deg: float, b_deg: float, r: float) -> np.ndarray:
    return float(r) * galactic_to_unit(l_deg, b_deg)

def cartesian_to_galactic(v: np.ndarray) -> tuple[float, float, float]:
    x, y, z = float(v[0]), float(v[1]), float(v[2])
    r = math.sqrt(x*x + y*y + z*z)
    if r == 0:
        return 0.0, 0.0, 0.0
    b = math.degrees(math.asin(max(-1.0, min(1.0, z/r))))
    l = math.degrees(math.atan2(y, x)) % 360.0
    return l, b, r

def random_cone_perturbation(u0: np.ndarray, sigma_deg: float, n: int, rng: np.random.Generator) -> np.ndarray:
    # isotropic Gaussian in angle around u0 (small-angle approx via tangent-plane)
    sigma = math.radians(float(sigma_deg))
    # choose an orthonormal basis (e1,e2) orthogonal to u0
    u0 = u0 / np.linalg.norm(u0)
    tmp = np.array([1.0, 0.0, 0.0], dtype=float)
    if abs(float(np.dot(tmp, u0))) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0], dtype=float)
    e1 = tmp - np.dot(tmp, u0) * u0
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(u0, e1)
    # sample small offsets
    a = rng.normal(0.0, sigma, size=n)
    b = rng.normal(0.0, sigma, size=n)
    u = u0[None, :] + a[:, None]*e1[None, :] + b[:, None]*e2[None, :]
    u /= np.linalg.norm(u, axis=1)[:, None]
    return u

# ---------------------------------------------------------------------
# Dipole observation + covariance
# ---------------------------------------------------------------------
COV_SAMPLES = 20000
COV_FLOOR = 1e-12
RNG = np.random.default_rng(12345)

@dataclass
class DipoleObservation:
    name: str
    variant: str
    survey: str
    l_deg: float
    b_deg: float
    D_obs: float
    sigma_mag: float
    sigma_dir_deg: float
    D_kin: float
    reference: str

    z_mean: float = float("nan")

    # cached
    D_obs_vec: np.ndarray = None
    D_kin_vec: np.ndarray = None
    cov: np.ndarray = None
    cov_inv: np.ndarray = None

    def build(self, rng: np.random.Generator) -> None:
        self.D_obs_vec = galactic_to_cartesian(self.l_deg, self.b_deg, self.D_obs)
        # kinematic expectation is assumed along CMB direction; for this harness, we treat its direction as fixed.
        # We'll use the canonical CMB galactic dipole direction:
        l_cmb, b_cmb = 264.021, 48.253  # standard value used widely (Planck)
        self.D_kin_vec = galactic_to_cartesian(l_cmb, b_cmb, self.D_kin)

        self.cov = self._estimate_covariance(rng=rng)
        self.cov = self.cov + np.eye(3) * COV_FLOOR
        self.cov_inv = np.linalg.inv(self.cov)

    def _estimate_covariance(self, rng: np.random.Generator) -> np.ndarray:
        u0 = galactic_to_unit(self.l_deg, self.b_deg)
        u_samp = random_cone_perturbation(u0, self.sigma_dir_deg, COV_SAMPLES, rng)

        mag = rng.normal(self.D_obs, self.sigma_mag, size=COV_SAMPLES)
        mag = np.clip(mag, 0.0, None)

        v = mag[:, None] * u_samp
        mean = np.mean(v, axis=0)
        dv = v - mean[None, :]
        cov = (dv.T @ dv) / (COV_SAMPLES - 1)
        return cov

# ---------------------------------------------------------------------
# Model fitting
# ---------------------------------------------------------------------
def chi2_from_residual(res: np.ndarray, cov_inv: np.ndarray) -> float:
    return float(res.T @ cov_inv @ res)

def fit_null(observations: List[DipoleObservation]) -> dict:
    chi2 = 0.0
    for obs in observations:
        res = obs.D_obs_vec - obs.D_kin_vec
        chi2 += chi2_from_residual(res, obs.cov_inv)
    n_data = 3 * len(observations)
    k = 0
    aic = chi2 + 2*k
    aicc = aic
    return dict(chi2=chi2, k=k, n=n_data, aic=aic, aicc=aicc)

def fit_rank1(observations: List[DipoleObservation]) -> dict:
    """
    Fit rank‑1 coherent excess:
      D_obs_i = D_kin_i + A_i * g_hat
    Parameters: two angles (l,b) for g_hat + N amplitudes A_i
    """
    N = len(observations)

    mean_excess = np.mean([obs.D_obs_vec - obs.D_kin_vec for obs in observations], axis=0)
    if np.linalg.norm(mean_excess) == 0:
        mean_excess = galactic_to_unit(264.021, 48.253)
    l0, b0, _ = cartesian_to_galactic(mean_excess)
    g0 = galactic_to_unit(l0, b0)
    A0 = np.array([float(np.dot(obs.D_obs_vec - obs.D_kin_vec, g0)) for obs in observations], dtype=float)
    x0 = np.concatenate([[l0, b0], A0])

    def objective(params: np.ndarray) -> float:
        l = float(params[0]) % 360.0
        b = float(np.clip(params[1], -90.0, 90.0))
        g = galactic_to_unit(l, b)
        A = params[2:]
        chi2 = 0.0
        for i, obs in enumerate(observations):
            pred = obs.D_kin_vec + A[i] * g
            res = obs.D_obs_vec - pred
            chi2 += chi2_from_residual(res, obs.cov_inv)
        return chi2

    opt = minimize(objective, x0, method="Nelder-Mead", options={"maxiter": 20000, "xatol": 1e-10, "fatol": 1e-10})
    l_fit = float(opt.x[0]) % 360.0
    b_fit = float(np.clip(opt.x[1], -90.0, 90.0))
    g_fit = galactic_to_unit(l_fit, b_fit)
    A_fit = np.array(opt.x[2:], dtype=float)

    n_data = 3*N
    k = 2 + N
    chi2 = float(opt.fun)
    aic = chi2 + 2*k
    # small-sample correction
    if n_data > k + 1:
        aicc = aic + (2*k*(k+1)) / (n_data - k - 1)
    else:
        aicc = float("nan")

    return dict(chi2=chi2, k=k, n=n_data, aic=aic, aicc=aicc, l=l_fit, b=b_fit, g=g_fit, A=A_fit, success=bool(opt.success))


# ---------------------------------------------------------------------
# BAO-anchored model (D-BAO)
# ---------------------------------------------------------------------
def E_z(z: float, Om0: float) -> float:
    """Dimensionless H(z)/H0 for flat LCDM."""
    return math.sqrt(Om0 * (1.0 + z)**3 + (1.0 - Om0))

def comoving_distance_mpc(z: float, H0: float = 70.0, Om0: float = 0.3, n_steps: int = 4000) -> float:
    """Comoving distance chi(z) in Mpc for flat LCDM."""
    if not np.isfinite(z) or z <= 0:
        return 0.0
    zs = np.linspace(0.0, float(z), int(n_steps) + 1)
    invE = 1.0 / np.array([E_z(float(zz), Om0) for zz in zs], dtype=float)
    integral = np.trapz(invE, zs)
    c_km_s = 299792.458
    return (c_km_s / H0) * float(integral)

def fit_bao_p1(
    observations: List[DipoleObservation],
    p: float = 1.0,
    r_bao_mpc: float = 147.0,
    H0: float = 70.0,
    Om0: float = 0.3,
) -> dict:
    """
    BAO-anchored coherent excess:
      D_obs_i = D_kin_i + beta * w_i * g_hat
    with w_i = (r_BAO / chi(z_mean_i))^p.

    Parameters: (l,b) for g_hat + beta (global scale).
    """
    N = len(observations)

    chis = np.zeros(N, dtype=float)
    w = np.zeros(N, dtype=float)
    for i, o in enumerate(observations):
        chi = comoving_distance_mpc(o.z_mean, H0=H0, Om0=Om0)
        chis[i] = chi
        if chi > 0:
            w[i] = (r_bao_mpc / chi) ** float(p)
        else:
            w[i] = 0.0

    mean_excess = np.mean([obs.D_obs_vec - obs.D_kin_vec for obs in observations], axis=0)
    if np.linalg.norm(mean_excess) == 0:
        mean_excess = galactic_to_unit(264.021, 48.253)
    l0, b0, _ = cartesian_to_galactic(mean_excess)
    g0 = galactic_to_unit(l0, b0)

    proj = np.array([float(np.dot(obs.D_obs_vec - obs.D_kin_vec, g0)) for obs in observations], dtype=float)
    denom = float(np.sum(w*w))
    beta0 = float(np.sum(proj * w) / denom) if denom > 0 else 0.0

    x0 = np.array([l0, b0, beta0], dtype=float)

    def objective(params: np.ndarray) -> float:
        l = float(params[0]) % 360.0
        b = float(np.clip(params[1], -90.0, 90.0))
        g = galactic_to_unit(l, b)
        beta = float(params[2])
        chi2 = 0.0
        for i, obs in enumerate(observations):
            pred = obs.D_kin_vec + (beta * w[i]) * g
            res = obs.D_obs_vec - pred
            chi2 += chi2_from_residual(res, obs.cov_inv)
        return chi2

    opt = minimize(objective, x0, method="Nelder-Mead", options={"maxiter": 20000, "xatol": 1e-10, "fatol": 1e-10})
    l_fit = float(opt.x[0]) % 360.0
    b_fit = float(np.clip(opt.x[1], -90.0, 90.0))
    g_fit = galactic_to_unit(l_fit, b_fit)
    beta_fit = float(opt.x[2])

    n_data = 3*N
    k = 3
    chi2 = float(opt.fun)
    aic = chi2 + 2*k
    if n_data > k + 1:
        aicc = aic + (2*k*(k+1)) / (n_data - k - 1)
    else:
        aicc = float("nan")

    return dict(
        chi2=chi2, k=k, n=n_data, aic=aic, aicc=aicc,
        l=l_fit, b=b_fit, g=g_fit,
        beta=beta_fit, p=float(p), r_bao_mpc=float(r_bao_mpc),
        H0=float(H0), Om0=float(Om0),
        chi=chis, w=w, A_eff=beta_fit*w,
        success=bool(opt.success)
    )


def coherence_angles(observations: List[DipoleObservation], g: np.ndarray) -> np.ndarray:
    g = g / np.linalg.norm(g)
    ang = []
    for obs in observations:
        ex = obs.D_obs_vec - obs.D_kin_vec
        if np.linalg.norm(ex) == 0:
            ang.append(float("nan"))
        else:
            exu = ex / np.linalg.norm(ex)
            ang.append(math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(exu, g)))))))
    return np.array(ang, dtype=float)

# ---------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------
def load_pack(path: str) -> list[dict]:
    if path.lower().endswith(".json"):
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    # CSV
    import csv
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)
    return rows

def to_float(x):
    if x is None:
        return float("nan")
    if isinstance(x, (int, float)):
        return float(x)
    s = str(x).strip()
    if s == "" or s.lower() == "nan":
        return float("nan")
    return float(s)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="CSV or JSON input pack file")
    ap.add_argument("--variant", default="baseline", help="which variant to select (baseline or 2MRS_clean)")
    ap.add_argument("--model", default="rank1", choices=["rank1","bao_p1"], help="model to fit vs NULL")
    ap.add_argument("--bao_p", type=float, default=1.0, help="BAO weight exponent p")
    ap.add_argument("--r_bao_mpc", type=float, default=147.0, help="BAO scale in Mpc")
    ap.add_argument("--H0", type=float, default=70.0, help="H0 in km/s/Mpc")
    ap.add_argument("--Om0", type=float, default=0.3, help="Omega_m")
    args = ap.parse_args()

    rows = load_pack(args.input)
    variant = args.variant

    obs: list[DipoleObservation] = []
    for r in rows:
        if str(r.get("variant", "")).strip() != variant:
            continue
        obs.append(DipoleObservation(
            name=str(r["name"]),
            variant=str(r["variant"]),
            survey=str(r["survey"]),
            l_deg=to_float(r["l_deg"]),
            b_deg=to_float(r["b_deg"]),
            D_obs=to_float(r["D_obs"]),
            sigma_mag=to_float(r["D_obs_err"]),
            sigma_dir_deg=to_float(r.get("sigma_dir_deg_approx", r.get("sigma_dir_deg", 10.0))),
            D_kin=to_float(r["D_kin_from_x_alpha"]),
            reference=str(r["reference"]),
            z_mean=to_float(r.get("z_mean", "nan")),
        ))

    if len(obs) < 2:
        raise SystemExit(f"Not enough tracers found for variant '{variant}'. Found {len(obs)}.")

    for o in obs:
        o.build(rng=RNG)

    # quick table
    print(f"{'Tracer':<12} {'Survey':<32} {'Dobs':>8} {'Dkin':>8} {'ratio':>7} {'(l,b)':>16} {'σmag':>9} {'σdir':>7}")
    print("-"*100)
    for o in obs:
        ratio = o.D_obs / o.D_kin if o.D_kin else float("nan")
        print(f"{o.name:<12} {o.survey[:32]:<32} {o.D_obs:>8.5f} {o.D_kin:>8.5f} {ratio:>7.2f} "
              f"({o.l_deg:>6.1f},{o.b_deg:>5.1f}) {o.sigma_mag:>9.5f} {o.sigma_dir_deg:>7.2f}")
    print()

    null = fit_null(obs)
    if args.model == "rank1":
        fit = fit_rank1(obs)
        model_label = "RANK1"
    else:
        fit = fit_bao_p1(obs, p=float(args.bao_p), r_bao_mpc=float(args.r_bao_mpc), H0=float(args.H0), Om0=float(args.Om0))
        model_label = "BAO_P1"

    dAIC  = fit["aic"]  - null["aic"]
    dAICc = fit["aicc"] - null["aicc"]

    ang = coherence_angles(obs, fit["g"])

    out = []
    out.append("PROTOCOL D — INPUT PACK RUN (HARNESS)")
    out.append("="*78)
    out.append(f"variant: {variant}")
    out.append(f"input: {args.input}")
    out.append("")
    out.append("MODEL COMPARISON")
    out.append(f"  NULL : chi2={null['chi2']:.6f}  k={null['k']}  AIC={null['aic']:.6f}  AICc={null['aicc']:.6f}")
    out.append(f"  {model_label}: chi2={fit['chi2']:.6f}  k={fit['k']}  AIC={fit['aic']:.6f}  AICc={fit['aicc']:.6f}")
    out.append(f"  ΔAIC  ({model_label}-null) = {dAIC:+.6f}")
    out.append(f"  ΔAICc ({model_label}-null) = {dAICc:+.6f}")
    out.append("")
    out.append(f"FITTED GLOBAL EXCESS DIRECTION ({model_label})")
    out.append(f"  g_hat (l,b) = ({fit['l']:.3f}°, {fit['b']:.3f}°)")
    out.append("")
    if args.model == "rank1":
        out.append("PER‑TRACER AMPLITUDES A_i (dimensionless)")
        for i,o in enumerate(obs):
            out.append(f"  {o.name:<12} A = {fit['A'][i]: .8g}")
        out.append("")
    else:
        out.append("BAO‑ANCHORED AMPLITUDES (effective)")
        out.append(f"  beta = {fit['beta']:.8g}   p = {fit['p']:.3f}   r_bao = {fit['r_bao_mpc']:.3f} Mpc")
        out.append(f"  cosmology for chi(z): H0={fit['H0']:.3f}  Om0={fit['Om0']:.3f}")
        out.append("  Per tracer: z_mean, chi[Mpc], w_i, A_eff")
        for i,o in enumerate(obs):
            out.append(f"  {o.name:<12} z={o.z_mean:>6.3f}  chi={fit['chi'][i]:>9.2f}  w={fit['w'][i]:>10.6g}  A_eff={fit['A_eff'][i]:> .8g}")
        out.append("")
    out.append("COHERENCE ANGLES: angle(excess_i, g_hat)  (0° = aligned)")
    for i,o in enumerate(obs):
        out.append(f"  {o.name:<12} angle = {ang[i]:6.2f}°")
    out.append(f"  mean = {float(np.nanmean(ang)):.2f}°   std = {float(np.nanstd(ang)):.2f}°")
    out.append("")
    out.append("CAVEATS")
    out.append("  - Direction uncertainty is an approximation derived from RA/Dec errors.")
    out.append("  - Full rigor requires using the original posterior/covariance from the catalogue dipole fit.")
    out.append("  - This harness assumes the kinematic dipole points exactly in the CMB dipole direction.")
    out.append("")

    out_txt = "\n".join(out)

    from pathlib import Path
    input_stem = Path(args.input).stem
    out_path = f"protocol_d_run_results_{input_stem}_{variant}.txt"
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(out_txt)

    print(out_txt)
    print(f"\n[OK] wrote {out_path}\nSaved: {out_path}")

if __name__ == "__main__":
    main()
