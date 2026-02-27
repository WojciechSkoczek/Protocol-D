#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

# Prefer astropy-healpix (works on Windows). If healpy is available (Linux), that's fine too.
try:
    import healpy as hp  # type: ignore
except Exception:
    try:
        from astropy_healpix import healpy as hp  # healpy-compatible API
    except Exception as e:
        raise SystemExit(
            "Missing HEALPix backend. Install: pip install astropy-healpix\n"
            "(healpy often fails on Windows; astropy-healpix is the supported alternative)."
        ) from e

try:
    from astropy.io import fits
    from astropy.coordinates import SkyCoord
    import astropy.units as u
except Exception as e:
    raise SystemExit("Missing dependency: astropy. Install with: pip install astropy") from e


CMB_L_DEG = 264.021
CMB_B_DEG = 48.253


def angle_deg(uvec: np.ndarray, vvec: np.ndarray) -> float:
    u = uvec / np.linalg.norm(uvec)
    v = vvec / np.linalg.norm(vvec)
    c = float(np.clip(np.dot(u, v), -1.0, 1.0))
    return float(np.degrees(np.arccos(c)))


def gal_to_unit(l_deg: float, b_deg: float) -> np.ndarray:
    l = np.radians(l_deg)
    b = np.radians(b_deg)
    x = np.cos(b) * np.cos(l)
    y = np.cos(b) * np.sin(l)
    z = np.sin(b)
    return np.array([x, y, z], dtype=float)


def fit_dipole_from_counts(nside: int, counts: np.ndarray, mask: np.ndarray) -> tuple[float, float, float]:
    """
    Fit dipole in counts on observed pixels only (mask True).
    Model: N = a0 + ax*x + ay*y + az*z, then dipole vector d = (ax,ay,az)/a0.
    Returns: (amp, l_deg, b_deg) in Galactic.
    """
    pix = np.where(mask)[0]
    if pix.size < 10:
        return float("nan"), float("nan"), float("nan")

    y = counts[pix].astype(float)

    x, yv, z = hp.pix2vec(nside, pix)
    A = np.column_stack([np.ones(pix.size), x, yv, z])

    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    a0, ax, ay, az = coef

    if not np.isfinite(a0) or a0 <= 0:
        return float("nan"), float("nan"), float("nan")

    dvec = np.array([ax, ay, az], dtype=float) / float(a0)
    amp = float(np.linalg.norm(dvec))
    if not np.isfinite(amp) or amp <= 0:
        return float("nan"), float("nan"), float("nan")

    l = float(np.degrees(np.arctan2(dvec[1], dvec[0])) % 360.0)
    b = float(np.degrees(np.arcsin(np.clip(dvec[2] / (amp + 1e-30), -1.0, 1.0))))
    return amp, l, b


def _parse_float_list(s: str) -> list[float]:
    xs = [float(x.strip()) for x in s.split(",") if x.strip()]
    return sorted(set(xs))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fits", required=True, help="Path to LoTSS_DR3_v1.0.srl.fits")
    ap.add_argument("--out", required=True, help="Output CSV path")

    # kept for backward compatibility / logging only (fit uses fit_nside)
    ap.add_argument("--nside", type=int, default=64, help="(Legacy) Work NSIDE for logging only. Fit uses --fit_nside.")
    ap.add_argument("--fit_nside", type=int, default=64, help="HEALPix NSIDE used for ALL counts/maps and dipole fit.")

    ap.add_argument("--flux_cuts_mjy", default="10,20,30,40,60,80", help="Comma-separated Total_flux cuts (mJy)")
    ap.add_argument("--bmins_deg", default="10,20,30", help="Comma-separated Galactic latitude cuts |b|>= (deg)")
    ap.add_argument(
        "--footprint_flux_mjy",
        type=float,
        default=80.0,
        help="Flux threshold used ONLY to build the footprint mask (mJy). "
             "If set higher than the lowest plateau cut, it will be clamped automatically.",
    )
    ap.add_argument(
        "--min_sources",
        type=int,
        default=1000,
        help="Minimum number of sources required to attempt a dipole fit (otherwise NaNs).",
    )
    args = ap.parse_args()

    fits_path = Path(args.fits)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    flux_cuts = _parse_float_list(args.flux_cuts_mjy)
    bmins = _parse_float_list(args.bmins_deg)

    fit_nside = int(args.fit_nside)
    if fit_nside <= 0 or (fit_nside & (fit_nside - 1)) != 0:
        raise SystemExit("[ERROR] fit_nside must be a power of two (e.g. 32, 64, 128).")

    if not flux_cuts:
        raise SystemExit("[ERROR] flux_cuts_mjy is empty.")
    if not bmins:
        raise SystemExit("[ERROR] bmins_deg is empty.")

    # --- Clamp footprint threshold to the lowest plateau cut ---
    min_plateau_cut = float(min(flux_cuts))
    footprint_flux = float(args.footprint_flux_mjy)
    if footprint_flux > min_plateau_cut:
        print(
            f"[WARN] footprint_flux_mjy={footprint_flux:g} mJy is higher than the lowest plateau cut "
            f"({min_plateau_cut:g} mJy). Clamping footprint_flux_mjy -> {min_plateau_cut:g} mJy "
            "(to avoid footprint instability at bright cuts)."
        )
        footprint_flux = min_plateau_cut
    # ----------------------------------------------------------

    # Read required columns only
    with fits.open(fits_path, memmap=True) as hdul:
        data = hdul[1].data
        ra = np.array(data["RA"], dtype=float)
        dec = np.array(data["DEC"], dtype=float)
        tf = np.array(data["Total_flux"], dtype=float)  # mJy
        m = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(tf) & (tf > 0)
        ra, dec, tf = ra[m], dec[m], tf[m]

    # Convert to Galactic
    sc = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
    gal = sc.galactic
    l = gal.l.deg.astype(float)
    b = gal.b.deg.astype(float)

    # HEALPix pixel indices at FIT NSIDE (resolution-invariant estimator)
    theta = np.radians(90.0 - b)
    phi = np.radians(l)
    pix_fit = hp.ang2pix(fit_nside, theta, phi, nest=False)

    cmb_u = gal_to_unit(CMB_L_DEG, CMB_B_DEG)

    rows: list[dict[str, object]] = []
    npix_fit = hp.nside2npix(fit_nside)

    for bmin in bmins:
        # Footprint on FIT NSIDE (stable)
        bright = (np.abs(b) >= bmin) & (tf >= footprint_flux)
        footprint_pix = np.unique(pix_fit[bright])
        mask_fit = np.zeros(npix_fit, dtype=bool)
        mask_fit[footprint_pix] = True
        f_sky = float(mask_fit.mean())

        prev_dir_u: np.ndarray | None = None
        prev_amp: float | None = None

        for fcut in flux_cuts:
            sel = (np.abs(b) >= bmin) & (tf >= fcut)
            sel = sel & mask_fit[pix_fit]  # enforce footprint

            n_src = int(sel.sum())
            if n_src < int(args.min_sources):
                amp, lfit, bfit = float("nan"), float("nan"), float("nan")
                ang_cmb = float("nan")
                d_prev = float("nan")
                d_amp = float("nan")
                dir_u = None
            else:
                counts = np.bincount(pix_fit[sel], minlength=npix_fit).astype(float)
                amp, lfit, bfit = fit_dipole_from_counts(fit_nside, counts, mask_fit)

                if np.isfinite(amp):
                    dir_u = gal_to_unit(lfit, bfit)
                    ang_cmb = angle_deg(dir_u, cmb_u)
                else:
                    dir_u = None
                    ang_cmb = float("nan")

                if dir_u is not None and prev_dir_u is not None:
                    d_prev = angle_deg(dir_u, prev_dir_u)
                else:
                    d_prev = float("nan")

                if np.isfinite(amp) and prev_amp is not None and np.isfinite(prev_amp):
                    d_amp = float(amp - prev_amp)
                else:
                    d_amp = float("nan")

            if dir_u is not None:
                prev_dir_u = dir_u
            if np.isfinite(amp):
                prev_amp = float(amp)

            rows.append(
                dict(
                    nside=int(args.nside),
                    fit_nside=int(fit_nside),
                    bmin_deg=float(bmin),
                    flux_cut_mjy=float(fcut),
                    footprint_flux_mjy=float(footprint_flux),
                    f_sky=float(f_sky),
                    n_src=int(n_src),
                    dipole_amp=float(amp),
                    dipole_l_deg=float(lfit),
                    dipole_b_deg=float(bfit),
                    angle_to_cmb_deg=float(ang_cmb),
                    delta_dir_deg_vs_prev=float(d_prev),
                    delta_amp_vs_prev=float(d_amp),
                )
            )

    cols = [
        "nside",
        "fit_nside",
        "bmin_deg",
        "flux_cut_mjy",
        "footprint_flux_mjy",
        "f_sky",
        "n_src",
        "dipole_amp",
        "dipole_l_deg",
        "dipole_b_deg",
        "angle_to_cmb_deg",
        "delta_dir_deg_vs_prev",
        "delta_amp_vs_prev",
    ]
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})

    print(f"[OK] Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
