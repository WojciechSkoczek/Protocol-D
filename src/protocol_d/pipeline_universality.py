#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Protocol D â€” Phase 1b universality runner (stable).
#
# Guarantees:
# - No code executes at import time (everything is inside main()).
# - Robust FITS path resolution: .fit <-> .fits and optional .gz.
# - Robust FITS column alias resolution (RA/Dec, magnitude/flux proxies, cz).
# - Produces artefacts:
#   * run_manifest.csv
#   * model_comparison.csv (H0/H1/H2/C1 using vMF on fitted dipole directions)
#   * transfer_tests.csv (dataset transfer, LODO, RA-holdout; frozen axis + refit kappa on test)
#
# CLI (PowerShell, from repo root):
#   $env:PYTHONPATH = (Resolve-Path .\src).Path
#   python -m protocol_d.pipeline_universality --config .\configs\datasets_phase1b_*.yaml --out .\results\universality_*

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.coordinates import SkyCoord, BarycentricTrueEcliptic
import astropy.units as u
from astropy_healpix import HEALPix


__RUNNER_VERSION__ = "v0.6.7-universality-clean"


# ----------------------------
# YAML
# ----------------------------

def load_yaml(path: str) -> Dict[str, Any]:
    try:
        import yaml  # type: ignore
    except Exception as e:
        raise RuntimeError("PyYAML is required (pip install pyyaml).") from e
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


# ----------------------------
# FITS path resolver
# ----------------------------

def _resolve_fits_path(path: str) -> str:
    """Resolve common FITS filename variants automatically."""
    p = Path(path)
    tried: List[str] = []

    def add(c: Path):
        s = str(c)
        if s not in tried:
            tried.append(s)

    def ok(c: Path) -> bool:
        try:
            return c.exists()
        except Exception:
            return False

    def swap_fit_fits(q: Path) -> Optional[Path]:
        n = q.name.lower()
        if n.endswith(".fits"):
            return q.with_suffix(".fit")
        if n.endswith(".fit"):
            return q.with_suffix(".fits")
        return None

    add(p)
    if ok(p):
        return str(p)

    if p.name.lower().endswith(".gz"):
        un = p.with_name(p.name[:-3])  # strip .gz
        add(un)
        if ok(un):
            return str(un)
        un_sw = swap_fit_fits(un)
        if un_sw is not None:
            add(un_sw)
            if ok(un_sw):
                return str(un_sw)
            sw_gz = p.with_name(un_sw.name + ".gz")
            add(sw_gz)
            if ok(sw_gz):
                return str(sw_gz)
    else:
        sw = swap_fit_fits(p)
        if sw is not None:
            add(sw)
            if ok(sw):
                return str(sw)

        gz = p.with_name(p.name + ".gz")
        add(gz)
        if ok(gz):
            return str(gz)

        if sw is not None:
            sw_gz = sw.with_name(sw.name + ".gz")
            add(sw_gz)
            if ok(sw_gz):
                return str(sw_gz)

    raise FileNotFoundError(f"Missing FITS: {path}. Tried: {tried}")


def load_fits_table(path: str) -> "fits.FITS_rec":
    path = _resolve_fits_path(path)
    with fits.open(path, memmap=True) as hdul:
        for hdu in hdul:
            if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                return hdu.data
    raise RuntimeError(f"No table HDU found in FITS: {path}")


# ----------------------------
# FITS column resolver
# ----------------------------

def _available_cols(tab: "fits.FITS_rec") -> List[str]:
    return list(getattr(tab, "columns").names)


def resolve_col(tab: "fits.FITS_rec", requested: Optional[str], aliases: Sequence[str], *, role: str) -> str:
    names = _available_cols(tab)
    if not names:
        raise KeyError("No columns found in FITS table.")

    if requested and requested in names:
        return requested

    lower_map = {n.lower(): n for n in names}
    if requested and requested.lower() in lower_map:
        return lower_map[requested.lower()]

    for a in aliases:
        if a in names:
            return a
        if a.lower() in lower_map:
            return lower_map[a.lower()]

    head = ", ".join(names[:60])
    more = "" if len(names) <= 60 else f", ... (+{len(names)-60} more)"
    raise KeyError(
        f"Could not resolve {role}. Requested '{requested}'. "
        f"Tried aliases: {list(aliases)[:25]}{'...' if len(aliases)>25 else ''}. "
        f"Available columns: {head}{more}"
    )


RA_ALIASES = ["RA", "RAJ2000", "RA_ICRS", "ALPHA_J2000", "_RAJ2000", "RAdeg", "raj2000"]
DEC_ALIASES = ["DEC", "DEJ2000", "DE_ICRS", "DELTA_J2000", "_DEJ2000", "DEdeg", "dej2000"]
CZ_ALIASES = ["cz", "CZ", "VHELIO", "Vhelio", "cz_helio", "cz_hel", "czHelio"]

MAG_ALIASES = [
    "w2mpro", "W2MPRO", "W2MAG", "w2mag", "W2", "w2",
    "W1_MAG", "W1MAG", "W1", "w1",
    "GMAG", "G_MAG", "gmag", "g_mag",
]
NVSS_FLUX_ALIASES = ["S1.4", "S1_4", "S1400", "S_1_4", "S1p4", "S_1p4"]
RACS_FLUX_ALIASES = ["Ftot", "FTOT", "F_TOTAL", "Fint", "FINT", "Sint", "S_int", "SINT"]


# ----------------------------
# Regions
# ----------------------------

def ra_jackknife_masks(ra_deg: np.ndarray, splits_deg: Sequence[Sequence[float]]) -> Dict[str, np.ndarray]:
    ra = np.asarray(ra_deg, dtype=float) % 360.0
    out: Dict[str, np.ndarray] = {}
    for pair in splits_deg:
        ra_min, ra_max = float(pair[0]) % 360.0, float(pair[1]) % 360.0
        if ra_min < ra_max:
            m = (ra >= ra_min) & (ra < ra_max)
        else:
            m = (ra >= ra_min) | (ra < ra_max)
        out[f"ra_{int(round(ra_min)):03d}_{int(round(ra_max)):03d}"] = m
    return out


def dec_split_masks(dec_deg: np.ndarray) -> Dict[str, np.ndarray]:
    dec = np.asarray(dec_deg, dtype=float)
    return {"dec_north": dec >= 0.0, "dec_south": dec < 0.0}


def ecliptic_lat_split_masks(ra_deg: np.ndarray, dec_deg: np.ndarray, beta_split_deg: float = 30.0) -> Dict[str, np.ndarray]:
    sc = SkyCoord(ra=np.asarray(ra_deg) * u.deg, dec=np.asarray(dec_deg) * u.deg, frame="icrs")
    ecl = sc.transform_to(BarycentricTrueEcliptic())
    beta = ecl.lat.to_value(u.deg)
    ab = np.abs(beta)
    b = float(beta_split_deg)
    return {f"ecl_absbeta_lt_{int(round(b))}": ab < b, f"ecl_absbeta_ge_{int(round(b))}": ab >= b}


# ----------------------------
# Geometry
# ----------------------------

def icrs_to_gal_l_b(ra_deg: np.ndarray, dec_deg: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    sc = SkyCoord(ra=np.asarray(ra_deg) * u.deg, dec=np.asarray(dec_deg) * u.deg, frame="icrs")
    gal = sc.galactic
    l = gal.l.to_value(u.deg) % 360.0
    b = gal.b.to_value(u.deg)
    return l, b


def lonlat_to_unitvec(l_deg: np.ndarray, b_deg: np.ndarray) -> np.ndarray:
    l = np.deg2rad(np.asarray(l_deg, dtype=float))
    b = np.deg2rad(np.asarray(b_deg, dtype=float))
    cb = np.cos(b)
    return np.vstack([cb * np.cos(l), cb * np.sin(l), np.sin(b)]).T


def unitvec_to_lb(v: np.ndarray) -> Tuple[float, float]:
    v = np.asarray(v, dtype=float)
    n = float(np.linalg.norm(v))
    if n == 0.0:
        return float("nan"), float("nan")
    x, y, z = (v / n).tolist()
    l = (math.degrees(math.atan2(y, x)) + 360.0) % 360.0
    b = math.degrees(math.asin(z))
    return l, b


# ----------------------------
# Dipole fit
# ----------------------------

@dataclass
class DipoleFitResult:
    a0: float
    ax: float
    ay: float
    az: float
    d_x: float
    d_y: float
    d_z: float
    dipole_amp: float
    dipole_l_deg: float
    dipole_b_deg: float
    n_pix_used: int


def fit_dipole_lsq_counts(
    l_deg: np.ndarray,
    b_deg: np.ndarray,
    *,
    map_nside: int,
    fit_nside: int,
    footprint_nside: int = 32,
    footprint_min_count: int = 1,
) -> DipoleFitResult:
    hp_map = HEALPix(nside=int(map_nside), order="ring", frame="galactic")
    hp_fit = HEALPix(nside=int(fit_nside), order="ring", frame="galactic")
    hp_foot = HEALPix(nside=int(footprint_nside), order="ring", frame="galactic")

    pix_foot = hp_foot.lonlat_to_healpix(np.asarray(l_deg) * u.deg, np.asarray(b_deg) * u.deg)
    counts_foot = np.bincount(np.asarray(pix_foot, dtype=int), minlength=hp_foot.npix)
    footprint_foot = counts_foot >= int(footprint_min_count)

    pix_all = np.arange(hp_map.npix, dtype=int)
    lon_map, lat_map = hp_map.healpix_to_lonlat(pix_all)
    pix_map_to_foot = hp_foot.lonlat_to_healpix(lon_map, lat_map)
    footprint_map = footprint_foot[np.asarray(pix_map_to_foot, dtype=int)]

    pix_map = hp_map.lonlat_to_healpix(np.asarray(l_deg) * u.deg, np.asarray(b_deg) * u.deg)
    counts_map = np.bincount(np.asarray(pix_map, dtype=int), minlength=hp_map.npix).astype(float)

    pix_use = pix_all[footprint_map]
    y = counts_map[pix_use]

    lon_use, lat_use = hp_map.healpix_to_lonlat(pix_use)
    pix_fit = hp_fit.lonlat_to_healpix(lon_use, lat_use)
    lon_fit, lat_fit = hp_fit.healpix_to_lonlat(pix_fit)
    X = lonlat_to_unitvec(lon_fit.to_value(u.deg), lat_fit.to_value(u.deg))

    A = np.column_stack([np.ones(len(pix_use)), X])
    coeffs, *_ = np.linalg.lstsq(A, y, rcond=None)
    a0, ax, ay, az = coeffs.tolist()

    if a0 == 0.0:
        d = np.array([np.nan, np.nan, np.nan])
        amp = float("nan")
        lhat, bhat = float("nan"), float("nan")
    else:
        d = np.array([ax, ay, az], dtype=float) / float(a0)
        amp = float(np.linalg.norm(d))
        lhat, bhat = unitvec_to_lb(d)

    return DipoleFitResult(
        a0=float(a0), ax=float(ax), ay=float(ay), az=float(az),
        d_x=float(d[0]), d_y=float(d[1]), d_z=float(d[2]),
        dipole_amp=float(amp), dipole_l_deg=float(lhat), dipole_b_deg=float(bhat),
        n_pix_used=int(len(pix_use)),
    )


# ----------------------------
# vMF direction models
# ----------------------------

def vmf_logC(kappa: float) -> float:
    k = float(kappa)
    if k <= 1e-12:
        return -math.log(4.0 * math.pi)
    return math.log(k) - math.log(4.0 * math.pi) - math.log(math.sinh(k))


def A3(kappa: float) -> float:
    k = float(kappa)
    if k <= 1e-8:
        return k / 3.0
    return (1.0 / math.tanh(k)) - (1.0 / k)


def kappa_from_A3(t: float) -> float:
    t = float(t)
    if not math.isfinite(t) or t <= 0.0:
        return 0.0
    if t >= 0.999999:
        return 1e6
    lo, hi = 0.0, 50.0
    while A3(hi) < t and hi < 1e6:
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if A3(mid) < t:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def fit_vmf_axis_and_kappa(unit_dirs: np.ndarray) -> Tuple[np.ndarray, float, float]:
    S = np.sum(unit_dirs, axis=0)
    nS = float(np.linalg.norm(S))
    N = int(unit_dirs.shape[0])
    if N == 0 or nS == 0.0:
        return np.array([np.nan, np.nan, np.nan]), 0.0, 0.0
    v = S / nS
    R = nS / float(N)
    k = kappa_from_A3(R)
    return v, float(k), float(R)


def fit_kappa_given_axis(unit_dirs: np.ndarray, axis: np.ndarray) -> float:
    axis = np.asarray(axis, dtype=float)
    if not np.isfinite(axis).all() or np.linalg.norm(axis) == 0.0:
        return 0.0
    axis = axis / np.linalg.norm(axis)
    t = float(np.mean(np.dot(unit_dirs, axis)))
    t = max(t, 0.0)
    return kappa_from_A3(min(t, 0.999999))


def vmf_loglik(unit_dirs: np.ndarray, axis: np.ndarray, kappa: float) -> float:
    N = int(unit_dirs.shape[0])
    if N == 0:
        return float("nan")
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    return N * vmf_logC(kappa) + float(kappa) * float(np.sum(np.dot(unit_dirs, axis)))


def aicc(logL: float, k: int, n: int) -> float:
    AIC = 2 * int(k) - 2 * float(logL)
    if n <= (k + 1):
        return float("nan")
    return float(AIC + (2 * k * (k + 1)) / (n - k - 1))


def _sha(s: str, n: int = 10) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()[:n]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_yaml(args.config)
    out_root = str(args.out)
    os.makedirs(out_root, exist_ok=True)
    rng = np.random.default_rng(12345)

    defaults = cfg.get("defaults", {})
    fit_nside = int(defaults.get("fit_nside", 64))
    map_nside = int(defaults.get("map_nside", 64))
    global_bmin = float(defaults.get("bmin_deg", 20.0))
    depth_q = defaults.get("depth_bins_quantiles", [0.0, 0.30, 0.70, 1.0])
    z_q = defaults.get("z_bins_quantiles", [0.0, 0.30, 0.70, 1.0])

    region_cfg = defaults.get("region_splits", {})
    ra_splits = region_cfg.get("ra_jackknife_deg", [[0, 90], [90, 180], [180, 270], [270, 360]])
    use_dec_split = bool(region_cfg.get("dec_split", False))
    beta_split = float(region_cfg.get("ecliptic_lat_split_deg", 30.0))

    manifest_rows: List[Dict[str, Any]] = []
    dir_rows: List[Dict[str, Any]] = []

    for ds in cfg.get("datasets", []):
        if not ds.get("enabled", True):
            continue

        ds_id = str(ds["id"])
        ds_label = str(ds.get("label", ds_id))
        input_path = str(ds["input"]["path"])

        tab = load_fits_table(input_path)

        req_ra = ds["columns"].get("ra_deg")
        req_dec = ds["columns"].get("dec_deg")
        ra_col = resolve_col(tab, req_ra, RA_ALIASES + ([str(req_ra)] if req_ra else []), role="RA column")
        dec_col = resolve_col(tab, req_dec, DEC_ALIASES + ([str(req_dec)] if req_dec else []), role="Dec column")

        ra = np.asarray(tab[ra_col], dtype=float)
        dec = np.asarray(tab[dec_col], dtype=float)

        l, b = icrs_to_gal_l_b(ra, dec)

        cuts = ds.get("cuts", {})
        bmin = float(cuts.get("bmin_deg", global_bmin))
        mask_b = np.abs(b) >= bmin
        mask_extra = np.ones_like(mask_b, dtype=bool)

        z_arr: Optional[np.ndarray] = None
        z_cfg = ds["columns"].get("z", None)
        if isinstance(z_cfg, dict) and str(z_cfg.get("kind", "")).lower() == "cz_kms":
            cz_req = z_cfg.get("column")
            cz_col = resolve_col(tab, cz_req, CZ_ALIASES + ([str(cz_req)] if cz_req else []), role="cz column")
            cz = np.asarray(tab[cz_col], dtype=float)
            cz_min = float(cuts.get("cz_kms_min", -np.inf))
            mask_extra &= (cz >= cz_min)
            z_arr = cz / float(z_cfg.get("c_kms", 299792.458))
        elif isinstance(z_cfg, str):
            z_col = resolve_col(tab, z_cfg, [z_cfg], role="z column")
            z_arr = np.asarray(tab[z_col], dtype=float)

        depth_arr: Optional[np.ndarray] = None
        depth_cfg = ds["columns"].get("depth_proxy", None)
        if isinstance(depth_cfg, dict):
            kind = str(depth_cfg.get("kind", "")).lower()
            req = depth_cfg.get("column")
            if kind == "magnitude":
                dcol = resolve_col(tab, req, MAG_ALIASES + ([str(req)] if req else []), role="depth proxy")
                depth_arr = np.asarray(tab[dcol], dtype=float)
                if not bool(depth_cfg.get("brighter_is_higher", True)):
                    depth_arr = -depth_arr
            elif kind == "flux_mjy":
                dcol = resolve_col(tab, req, NVSS_FLUX_ALIASES + RACS_FLUX_ALIASES + ([str(req)] if req else []), role="depth proxy")
                depth_arr = np.asarray(tab[dcol], dtype=float)
                fmin = cuts.get("flux_mjy_min", None)
                if fmin is not None:
                    mask_extra &= (depth_arr >= float(fmin))

        zr = cuts.get("z_range", None)
        if zr is not None and z_arr is not None and isinstance(zr, (list, tuple)) and len(zr) == 2:
            zlo, zhi = float(zr[0]), float(zr[1])
            mask_extra &= (z_arr >= zlo) & (z_arr <= zhi)

        base_mask = mask_b & mask_extra
        if int(base_mask.sum()) < 1000:
            continue

        regions: Dict[str, np.ndarray] = {}
        regions.update(ra_jackknife_masks(ra, ra_splits))
        if use_dec_split:
            regions.update(dec_split_masks(dec))
        if str(ds.get("family", "")).lower().startswith("wise"):
            regions.update(ecliptic_lat_split_masks(ra, dec, beta_split_deg=beta_split))

        binning = ds.get("binning", {})
        use_z_bins = (z_arr is not None) and (binning.get("z_bins", "none") == "quantiles")
        use_depth_bins = (depth_arr is not None) and (binning.get("depth_bins", "none") == "quantiles")

        if use_z_bins:
            xbin = z_arr[base_mask]
            quantiles = z_q
            bin_kind = "z"
        elif use_depth_bins:
            xbin = depth_arr[base_mask]
            quantiles = depth_q
            bin_kind = "depth"
        else:
            xbin = None
            quantiles = [0.0, 1.0]
            bin_kind = "none"

        if xbin is None:
            edges = [float("-inf"), float("inf")]
            bins = [("all", np.ones(int(base_mask.sum()), dtype=bool))]
        else:
            qs = np.quantile(xbin, quantiles)
            edges = [float(v) for v in qs.tolist()]
            bins: List[Tuple[str, np.ndarray]] = []
            for i in range(len(edges) - 1):
                lo, hi = edges[i], edges[i + 1]
                if i < len(edges) - 2:
                    m = (xbin >= lo) & (xbin < hi)
                else:
                    m = (xbin >= lo) & (xbin <= hi)
                bid = f"q{int(100*quantiles[i]):02d}_{int(100*quantiles[i+1]):02d}"
                bins.append((bid, m))

        idx_base = np.where(base_mask)[0]

        for region_id, rmask in regions.items():
            mask_region = base_mask & rmask
            if int(mask_region.sum()) < 2000:
                continue
            idx_region_in_base = np.isin(idx_base, np.where(mask_region)[0])

            for bid, bmask_in_base in bins:
                m_in_base = idx_region_in_base & bmask_in_base
                idx_cell = idx_base[m_in_base]
                if int(idx_cell.size) < 2000:
                    continue

                res = fit_dipole_lsq_counts(
                    l[idx_cell], b[idx_cell],
                    map_nside=map_nside, fit_nside=fit_nside,
                )

                run_key = json.dumps({
                    "dataset": ds_id, "region": region_id, "bin_kind": bin_kind, "bin": bid,
                    "bmin": bmin, "fit_nside": fit_nside, "map_nside": map_nside
                }, sort_keys=True)
                run_id = f"{ds_id}__{region_id}__{bin_kind}_{bid}__{_sha(run_key)}"

                out_dir = os.path.join(out_root, ds_id, run_id)
                os.makedirs(out_dir, exist_ok=True)
                fit_path = os.path.join(out_dir, "fit.csv")

                pd.DataFrame([{
                    "a0": res.a0, "ax": res.ax, "ay": res.ay, "az": res.az,
                    "d_x": res.d_x, "d_y": res.d_y, "d_z": res.d_z,
                    "dipole_amp": res.dipole_amp,
                    "dipole_l_deg": res.dipole_l_deg,
                    "dipole_b_deg": res.dipole_b_deg,
                    "fit_method": "LSQ_cartesian_counts",
                    "fit_nside": fit_nside,
                    "map_nside": map_nside,
                    "n_objects": int(idx_cell.size),
                    "n_pix_used": res.n_pix_used,
                    "bmin_deg": bmin,
                    "dataset_id": ds_id,
                    "region_id": region_id,
                    "bin_kind": bin_kind,
                    "bin_id": bid,
                    "resolved_ra_col": ra_col,
                    "resolved_dec_col": dec_col,
                }]).to_csv(fit_path, index=False)

                manifest_rows.append({
                    "run_id": run_id,
                    "dataset_id": ds_id,
                    "dataset_label": ds_label,
                    "input_path": input_path,
                    "region_id": region_id,
                    "bin_kind": bin_kind,
                    "bin_id": bid,
                    "bin_edges": json.dumps(edges),
                    "bmin_deg": bmin,
                    "fit_nside": fit_nside,
                    "map_nside": map_nside,
                    "n_objects": int(idx_cell.size),
                    "output_dir": out_dir,
                    "fit_result_path": fit_path,
                })

                if np.isfinite(res.dipole_l_deg) and np.isfinite(res.dipole_b_deg):
                    dir_rows.append({
                        "run_id": run_id,
                        "dataset_id": ds_id,
                        "dataset_label": ds_label,
                        "region_id": region_id,
                        "bin_kind": bin_kind,
                        "bin_id": bid,
                        "dipole_l_deg": res.dipole_l_deg,
                        "dipole_b_deg": res.dipole_b_deg,
                        "dipole_amp": res.dipole_amp,
                    })

    pd.DataFrame(manifest_rows).to_csv(os.path.join(out_root, "run_manifest.csv"), index=False)

    dirs_df = pd.DataFrame(dir_rows)
    if len(dirs_df) < 10:
        raise RuntimeError("Too few fitted cells for global comparison.")

    unit_dirs = lonlat_to_unitvec(dirs_df["dipole_l_deg"].to_numpy(), dirs_df["dipole_b_deg"].to_numpy())
    N = int(unit_dirs.shape[0])
    D = int(dirs_df["dataset_id"].nunique())

    logL0 = N * (-math.log(4.0 * math.pi))
    aicc0 = aicc(logL0, k=0, n=N)

    v1, k1, _ = fit_vmf_axis_and_kappa(unit_dirs)
    logL1 = vmf_loglik(unit_dirs, v1, k1)
    aicc1 = aicc(logL1, k=3, n=N)

    k2 = k1
    logL2 = 0.0
    for _dsid, sub in dirs_df.groupby("dataset_id", sort=True):
        u_ds = lonlat_to_unitvec(sub["dipole_l_deg"].to_numpy(), sub["dipole_b_deg"].to_numpy())
        v_ds, _, _ = fit_vmf_axis_and_kappa(u_ds)
        logL2 += vmf_loglik(u_ds, v_ds, k2)
    aicc2 = aicc(logL2, k=(2 * D + 1), n=N)

    vr = rng.normal(size=3)
    vr = vr / np.linalg.norm(vr)
    kr = fit_kappa_given_axis(unit_dirs, vr)
    logLr = vmf_loglik(unit_dirs, vr, kr)
    aiccr = aicc(logLr, k=3, n=N)

    mc = pd.DataFrame([
        {"run_id": "GLOBAL_ALL", "model_id": "H0_isotropic", "k_params": 0, "n_objects": int(N),
         "logL": float(logL0), "AICc": float(aicc0), "axis_l_deg": float("nan"), "axis_b_deg": float("nan"), "kappa": 0.0},
        {"run_id": "GLOBAL_ALL", "model_id": "H1_shared_axis", "k_params": 3, "n_objects": int(N),
         "logL": float(logL1), "AICc": float(aicc1), "axis_l_deg": unitvec_to_lb(v1)[0], "axis_b_deg": unitvec_to_lb(v1)[1], "kappa": float(k1)},
        {"run_id": "GLOBAL_ALL", "model_id": "H2_dataset_axes", "k_params": int(2 * D + 1), "n_objects": int(N),
         "logL": float(logL2), "AICc": float(aicc2), "axis_l_deg": float("nan"), "axis_b_deg": float("nan"), "kappa": float(k2)},
        {"run_id": "GLOBAL_ALL", "model_id": "C1_random_axis", "k_params": 3, "n_objects": int(N),
         "logL": float(logLr), "AICc": float(aiccr), "axis_l_deg": unitvec_to_lb(vr)[0], "axis_b_deg": unitvec_to_lb(vr)[1], "kappa": float(kr)},
    ])
    mc["delta_AICc"] = mc["AICc"] - mc["AICc"].min()
    mc.to_csv(os.path.join(out_root, "model_comparison.csv"), index=False)

    transfers: List[Dict[str, Any]] = []

    for train_ds, train_sub in dirs_df.groupby("dataset_id", sort=True):
        u_tr = lonlat_to_unitvec(train_sub["dipole_l_deg"].to_numpy(), train_sub["dipole_b_deg"].to_numpy())
        v_tr, _, _ = fit_vmf_axis_and_kappa(u_tr)

        for test_ds, test_sub in dirs_df.groupby("dataset_id", sort=True):
            if test_ds == train_ds:
                continue
            u_te = lonlat_to_unitvec(test_sub["dipole_l_deg"].to_numpy(), test_sub["dipole_b_deg"].to_numpy())
            nT = int(u_te.shape[0])
            if nT < 5:
                continue

            k_te = fit_kappa_given_axis(u_te, v_tr)
            logLf = vmf_loglik(u_te, v_tr, k_te)

            logL0T = nT * (-math.log(4.0 * math.pi))
            a_f = aicc(logLf, k=1, n=nT)
            a_0 = aicc(logL0T, k=0, n=nT)
            dA = float(a_f - a_0)

            transfers.append({
                "transfer_id": f"ds_{train_ds}_to_{test_ds}",
                "train_spec": json.dumps({"dataset_id": train_ds}),
                "test_spec": json.dumps({"dataset_id": test_ds}),
                "axis_l_deg": unitvec_to_lb(v_tr)[0],
                "axis_b_deg": unitvec_to_lb(v_tr)[1],
                "kappa_test": float(k_te),
                "AICc_frozen": float(a_f),
                "AICc_baseline": float(a_0),
                "delta_AICc_vs_H0": dA,
                "n_objects_test": nT,
                "passed": bool(dA <= -10.0),
            })

    for held_ds in sorted(dirs_df["dataset_id"].unique().tolist()):
        train_mask = dirs_df["dataset_id"] != held_ds
        test_mask = dirs_df["dataset_id"] == held_ds

        u_tr = lonlat_to_unitvec(dirs_df.loc[train_mask, "dipole_l_deg"].to_numpy(), dirs_df.loc[train_mask, "dipole_b_deg"].to_numpy())
        u_te = lonlat_to_unitvec(dirs_df.loc[test_mask, "dipole_l_deg"].to_numpy(), dirs_df.loc[test_mask, "dipole_b_deg"].to_numpy())
        nT = int(u_te.shape[0])
        if nT < 5:
            continue

        v_tr, _, _ = fit_vmf_axis_and_kappa(u_tr)
        k_te = fit_kappa_given_axis(u_te, v_tr)
        logLf = vmf_loglik(u_te, v_tr, k_te)

        logL0T = nT * (-math.log(4.0 * math.pi))
        a_f = aicc(logLf, k=1, n=nT)
        a_0 = aicc(logL0T, k=0, n=nT)
        dA = float(a_f - a_0)

        transfers.append({
            "transfer_id": f"lodo_holdout_{held_ds}",
            "train_spec": json.dumps({"held_out_dataset": held_ds}),
            "test_spec": json.dumps({"dataset_id": held_ds}),
            "axis_l_deg": unitvec_to_lb(v_tr)[0],
            "axis_b_deg": unitvec_to_lb(v_tr)[1],
            "kappa_test": float(k_te),
            "AICc_frozen": float(a_f),
            "AICc_baseline": float(a_0),
            "delta_AICc_vs_H0": dA,
            "n_objects_test": nT,
            "passed": bool(dA <= -10.0),
        })

    for held in sorted(dirs_df["region_id"].unique().tolist()):
        if not str(held).startswith("ra_"):
            continue

        train_mask = dirs_df["region_id"] != held
        test_mask = dirs_df["region_id"] == held

        u_tr = lonlat_to_unitvec(dirs_df.loc[train_mask, "dipole_l_deg"].to_numpy(), dirs_df.loc[train_mask, "dipole_b_deg"].to_numpy())
        u_te = lonlat_to_unitvec(dirs_df.loc[test_mask, "dipole_l_deg"].to_numpy(), dirs_df.loc[test_mask, "dipole_b_deg"].to_numpy())
        nT = int(u_te.shape[0])
        if nT < 5:
            continue

        v_tr, _, _ = fit_vmf_axis_and_kappa(u_tr)
        k_te = fit_kappa_given_axis(u_te, v_tr)
        logLf = vmf_loglik(u_te, v_tr, k_te)

        logL0T = nT * (-math.log(4.0 * math.pi))
        a_f = aicc(logLf, k=1, n=nT)
        a_0 = aicc(logL0T, k=0, n=nT)
        dA = float(a_f - a_0)

        transfers.append({
            "transfer_id": f"ra_holdout_{held}",
            "train_spec": json.dumps({"held_out_region": held}),
            "test_spec": json.dumps({"region_id": held}),
            "axis_l_deg": unitvec_to_lb(v_tr)[0],
            "axis_b_deg": unitvec_to_lb(v_tr)[1],
            "kappa_test": float(k_te),
            "AICc_frozen": float(a_f),
            "AICc_baseline": float(a_0),
            "delta_AICc_vs_H0": dA,
            "n_objects_test": nT,
            "passed": bool(dA <= -10.0),
        })

    pd.DataFrame(transfers).to_csv(os.path.join(out_root, "transfer_tests.csv"), index=False)

    print(f"Wrote artefacts to: {out_root}")


if __name__ == "__main__":
    main()

