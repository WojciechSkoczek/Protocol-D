from __future__ import annotations

import argparse
import math
import pathlib
import yaml
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Dict, List, Tuple
from scipy.optimize import minimize

import re

from .core import (
    Observation, galactic_to_unit, cartesian_to_galactic, angle_deg,
    bao_weight, robust_whitened_loss
)

# ----------------------------
# O2 operator model
# ----------------------------
def parse_magcut_to_float(magcut) -> float | None:
    if magcut is None or (isinstance(magcut, float) and np.isnan(magcut)):
        return None
    s = str(magcut).strip()
    # accept "G20.0" or numeric
    if re.match(r"^[0-9]+(\.[0-9]+)?$", s):
        return float(s)
    m = re.search(r"([0-9]+\.[0-9]+|[0-9]+)", s)
    return float(m.group(1)) if m else None


def A_of(bmin: float | None, magcut: float | None, a0: float, a_bmin: float, a_mag: float) -> float:
    b = float(bmin) if bmin is not None and np.isfinite(bmin) else 30.0
    m = float(magcut) if magcut is not None and np.isfinite(magcut) else 20.0
    return a0 + a_bmin*(b-30.0) + a_mag*(m-20.0)

def fit_O2_single(observations: List[Observation], cfg: dict) -> dict:
    # params: g_l,g_b,beta, s_l,s_b, a0,a_bmin,a_mag
    rng = np.random.default_rng(123)
    for o in observations:
        o.build(
            rng=rng,
            l_cmb=cfg["cmb"]["l"], b_cmb=cfg["cmb"]["b"],
            n_samples=cfg["cov"]["n_samples"],
            cov_floor=cfg["cov"]["floor"],
        )
    # init from mean excess
    mean_excess = np.mean([o.D_obs_vec - o.D_kin_vec for o in observations], axis=0)
    if np.linalg.norm(mean_excess) == 0:
        mean_excess = galactic_to_unit(cfg["cmb"]["l"], cfg["cmb"]["b"])
    g_l0,g_b0,_ = cartesian_to_galactic(mean_excess)
    # operator axis init: roughly orthogonal to g, but deterministic
    s_l0, s_b0 = (g_l0+90.0)%360.0, max(-70.0, min(70.0, g_b0))
    x0 = np.array([g_l0,g_b0, 0.22, s_l0,s_b0, 0.015, -0.0006, 0.0055], dtype=float)

    huber_delta = cfg["robust"]["huber_delta"] if cfg["robust"]["use_huber"] else None

    def objective(p: np.ndarray) -> float:
        g_l = float(p[0])%360.0
        g_b = float(np.clip(p[1], -90.0, 90.0))
        beta = float(p[2])
        s_l = float(p[3])%360.0
        s_b = float(np.clip(p[4], -90.0, 90.0))
        a0,a_bmin,a_mag = map(float, p[5:8])
        g = galactic_to_unit(g_l,g_b)
        s = galactic_to_unit(s_l,s_b)
        chi2 = 0.0
        for o in observations:
            w = bao_weight(o.z_mean, r_bao_mpc=cfg["cosmo"]["r_bao_mpc"], p=cfg["cosmo"]["bao_p"], H0=cfg["cosmo"]["H0"], Om0=cfg["cosmo"]["Om0"])
            pred = o.D_kin_vec + beta*w*g + float(o.t_quaia)*A_of(o.bmin,o.magcut,a0,a_bmin,a_mag)*s
            res = o.D_obs_vec - pred
            chi2 += robust_whitened_loss(res, o.cov_inv, huber_delta)
        return float(chi2)

    opt = minimize(objective, x0, method="Nelder-Mead", options={"maxiter": 50000, "xatol": 1e-10, "fatol": 1e-10})
    p = opt.x
    g_l = float(p[0])%360.0; g_b=float(np.clip(p[1],-90,90)); beta=float(p[2])
    s_l = float(p[3])%360.0; s_b=float(np.clip(p[4],-90,90))
    a0,a_bmin,a_mag = map(float, p[5:8])
    g = galactic_to_unit(g_l,g_b); s=galactic_to_unit(s_l,s_b)
    k=8; n=3*len(observations)
    chi2=float(opt.fun)
    aic=chi2+2*k
    aicc=aic + (2*k*(k+1))/(n-k-1) if n>k+1 else float("nan")

    return {
        "success": bool(opt.success),
        "chi2": chi2, "aic": aic, "aicc": aicc, "k": k, "n": n,
        "g_l": g_l, "g_b": g_b, "beta": beta,
        "s_l": s_l, "s_b": s_b,
        "a0": a0, "a_bmin": a_bmin, "a_mag": a_mag,
        "g_to_CMB": angle_deg(g, galactic_to_unit(cfg["cmb"]["l"], cfg["cmb"]["b"])),
        "s_to_CMB": angle_deg(s, galactic_to_unit(cfg["cmb"]["l"], cfg["cmb"]["b"])),
    }

def assemble_observations(pack_csv: pathlib.Path, variant: str, quaia_maskscan_csv: pathlib.Path, cfg: dict) -> List[Observation]:
    # From pack: CatWISE, LoTSS, NVSS, Quaia zbins (3)
    rows = pd.read_csv(pack_csv)
    rows = rows[rows["variant"]==variant].copy()

    obs: List[Observation] = []
    for _,r in rows.iterrows():
        name = str(r["name"])
        t = 1.0 if name.lower().startswith("quaia") else 0.0
        # infer magcut numeric for quaia; otherwise None
        mag = None
        if t==1.0:
            # try extract "G<20.0" from name
            m = re.search(r"G<\s*([0-9]+\.?[0-9]*)", name)
            if m:
                mag = float(m.group(1))
            else:
                # fallback: from 'magcut' string if present
                mag = None
        obs.append(Observation(
            name=name,
            l_deg=float(r["l_deg"]), b_deg=float(r["b_deg"]),
            D_obs=float(r["D_obs"]),
            sigma_mag=float(r["D_obs_err"]),
            sigma_dir_deg=float(r.get("sigma_dir_deg_approx", 10.0)),
            D_kin=float(r["D_kin_from_x_alpha"]),
            z_mean=float(r.get("z_mean", np.nan)),
            t_quaia=t,
            bmin=(cfg["operator"]["default_quaia_bmin_for_zbins"] if t==1.0 else None),
            magcut=mag,
        ))

    # Add maskscan points as additional Quaia constraints
    ms = pd.read_csv(quaia_maskscan_csv)
    # magcut column is like "G20.0"
    for _,r in ms.iterrows():
        mag = float(str(r["magcut"]).replace("G",""))
        bmin = float(r["bmin"])
        obs.append(Observation(
            name=f"Quaia_maskscan_{r['magcut']}_b{int(bmin)}",
            l_deg=float(r["l_deg"]), b_deg=float(r["b_deg"]),
            D_obs=float(r["D_obs"]),
            sigma_mag=float(cfg["operator"]["maskscan_error_floors"]["sigma_mag"]),
            sigma_dir_deg=float(cfg["operator"]["maskscan_error_floors"]["sigma_dir_deg"]),
            D_kin=0.0,   # maskscan table is for D_obs only; kinematic term will be handled via beta*w*g (we set D_kin=0 so E=D_obs)
            z_mean=float(r.get("z_mean", np.nan)),
            t_quaia=1.0,
            bmin=bmin,
            magcut=mag,
        ))
    return obs

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()

    cfg = yaml.safe_load(pathlib.Path(args.config).read_text())
    cfg2 = {
        "cosmo": cfg["cosmology"],
        "cmb": cfg["cosmology"]["cmb_dipole_galactic_deg"],
        "cov": {
            "n_samples": cfg["covariance"]["n_samples"],
            "floor": cfg["covariance"]["floor"],
        },
        "robust": cfg["covariance"]["robust_loss"],
        "operator": {
            "default_quaia_bmin_for_zbins": cfg["operator_model_O2"]["default_quaia_bmin_for_zbins"],
            "maskscan_error_floors": cfg["maskscan_error_floors"],
        },
    }


    out_dir = pathlib.Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Default pack mapping: assumes packs are present from cleaned bundle
    pack_map = {
        ("BLM_QG20", "baseline"):   "data/packs/03_BLM_QG20_CW_ZR.csv",
        ("BLM_QG20", "2MRS_clean"): "data/packs/03_BLM_QG20_CW_ZR.csv",
        ("WGV_QG20", "baseline"):   "data/packs/06_WGV_QG20_CW_ZR.csv",
        ("WGV_QG20", "2MRS_clean"): "data/packs/06_WGV_QG20_CW_ZR.csv",
    }

    maskscan = pathlib.Path("data/aux_data/quaia_dipole_maskscan.csv")

    # quick existence check (helps with fresh-zip reproducibility)
    for (_pack, _variant), _path in pack_map.items():
        if not pathlib.Path(_path).exists():
            raise SystemExit(f"Missing pack file: {_path} (check data/packs/)")
    if not maskscan.exists():
        raise SystemExit(f"Missing maskscan file: {maskscan} (check data/aux_data/)")

    rows=[]
    for (pack,variant), path in pack_map.items():
        obs = assemble_observations(pathlib.Path(path), variant, maskscan, cfg2)
        fit = fit_O2_single(obs, cfg2)
        fit.update({"pack": pack, "variant": variant, "N_obs": len(obs)})
        rows.append(fit)

    df = pd.DataFrame(rows)
    df.to_csv(out_dir/"repro_operator_O2_single_fits.csv", index=False)

    # Also write a small markdown note
    md = ["# Operator reproduction (O2)", "", "This run reproduces the *structure* of Runs K–Q using derived inputs.",
          "", "If your AICc numbers do not match the canonical snapshots exactly, check:",
          "- maskscan error floors in configs/defaults.yaml",
          "- whether you want D_kin included for maskscan points (default here: D_kin=0; treat as pure observed vector)"]
    (out_dir/"README_repro_operator.md").write_text("\n".join(md), encoding="utf-8")

if __name__ == "__main__":
    main()
