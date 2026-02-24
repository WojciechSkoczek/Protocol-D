from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import yaml
from scipy.optimize import minimize

from .core import (
    Observation,
    angle_deg,
    bao_weight,
    galactic_to_unit,
    robust_whitened_loss,
)

_RE_BMIN = re.compile(r"\|b\|\s*[≥>=]\s*([0-9]+)", re.IGNORECASE)
_RE_GCUT = re.compile(r"\bG<\s*([0-9]+(?:\.[0-9]+)?)", re.IGNORECASE)
_RE_W1CUT = re.compile(r"\bW1<\s*([0-9]+(?:\.[0-9]+)?)", re.IGNORECASE)


def _to_float(x: Any) -> float:
    try:
        s = str(x).strip()
        if s == "" or s.lower() == "nan":
            return float("nan")
        return float(s)
    except Exception:
        return float("nan")


def _is_quaia(name: str, survey: str) -> bool:
    s = f"{name} {survey}".lower()
    return "quaia" in s


def _parse_bmin_magcut(cut: str, default_bmin: float = 30.0, default_magcut: float = 20.0) -> Tuple[float, float]:
    bmin = default_bmin
    magcut = default_magcut
    if cut:
        m = _RE_BMIN.search(cut)
        if m:
            bmin = float(m.group(1))
        m = _RE_GCUT.search(cut)
        if m:
            magcut = float(m.group(1))
        else:
            m = _RE_W1CUT.search(cut)
            if m:
                magcut = float(m.group(1))
    return bmin, magcut


def _aicc(loss: float, k: int, n: int) -> float:
    aic = float(loss) + 2.0 * float(k)
    if n > (k + 1):
        return aic + (2.0 * k * (k + 1)) / (n - k - 1)
    return float("nan")


def load_pack_rows(pack_csv: Path, variant: str) -> List[Dict[str, str]]:
    out: List[Dict[str, str]] = []
    with pack_csv.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            if (row.get("variant") or "").strip() == variant:
                out.append(row)
    return out


def build_observations(rows: List[Dict[str, str]], cfg: Dict[str, Any], seed: int) -> List[Observation]:
    rng = np.random.default_rng(seed)

    l_cmb = float(cfg["cosmology"]["cmb_dipole_galactic_deg"]["l"])
    b_cmb = float(cfg["cosmology"]["cmb_dipole_galactic_deg"]["b"])

    n_samples = int(cfg["covariance"]["n_samples"])
    cov_floor = float(cfg["covariance"]["floor"])

    sigma_mag_floor = float(cfg["maskscan_error_floors"]["sigma_mag"])
    sigma_dir_floor = float(cfg["maskscan_error_floors"]["sigma_dir_deg"])

    default_quaia_bmin = float(cfg["operator_model_O2"].get("default_quaia_bmin_for_zbins", 30))

    obs: List[Observation] = []
    for row in rows:
        name = str(row.get("name") or "")
        survey = str(row.get("survey") or "")
        cut = str(row.get("flux_or_mag_cut") or "")

        t = 1.0 if _is_quaia(name, survey) else 0.0
        bmin, magcut = _parse_bmin_magcut(cut, default_bmin=default_quaia_bmin, default_magcut=20.0)

        D_obs = _to_float(row.get("D_obs"))
        D_err = _to_float(row.get("D_obs_err"))
        sig_dir = _to_float(row.get("sigma_dir_deg_approx"))
        D_kin = _to_float(row.get("D_kin_from_x_alpha"))
        z_mean = _to_float(row.get("z_mean"))
        l_deg = _to_float(row.get("l_deg"))
        b_deg = _to_float(row.get("b_deg"))

        sigma_mag = max(float(D_err), sigma_mag_floor) if np.isfinite(D_err) else sigma_mag_floor
        sigma_dir = max(float(sig_dir), sigma_dir_floor) if np.isfinite(sig_dir) else sigma_dir_floor

        o = Observation(
            name=name,
            l_deg=float(l_deg),
            b_deg=float(b_deg),
            D_obs=float(D_obs),
            sigma_mag=float(sigma_mag),
            sigma_dir_deg=float(sigma_dir),
            D_kin=float(D_kin),
            z_mean=float(z_mean),
            t_quaia=float(t),
            bmin=float(bmin),
            magcut=float(magcut),
        )
        o.build(rng=rng, l_cmb=l_cmb, b_cmb=b_cmb, n_samples=n_samples, cov_floor=cov_floor)
        obs.append(o)

    return obs


def fit_bao_p1(obs: List[Observation], cfg: Dict[str, Any]) -> Dict[str, float]:
    H0 = float(cfg["cosmology"]["H0"])
    Om0 = float(cfg["cosmology"]["Om0"])
    r_bao = float(cfg["cosmology"]["r_bao_mpc"])
    p = float(cfg["cosmology"]["bao_p"])

    use_huber = bool(cfg["robust_loss"]["use_huber"])
    huber_delta = float(cfg["robust_loss"]["huber_delta"]) if use_huber else None

    w = np.array([bao_weight(o.z_mean, r_bao_mpc=r_bao, p=p, H0=H0, Om0=Om0) for o in obs], dtype=float)

    mean_ex = np.mean([o.D_obs_vec - o.D_kin_vec for o in obs], axis=0)
    g0 = galactic_to_unit(264.021, 48.253) if np.linalg.norm(mean_ex) == 0 else mean_ex / np.linalg.norm(mean_ex)

    g_l0 = float(np.degrees(np.arctan2(g0[1], g0[0])) % 360.0)
    g_b0 = float(np.degrees(np.arcsin(np.clip(g0[2], -1.0, 1.0))))

    proj = np.array([float(np.dot((o.D_obs_vec - o.D_kin_vec), g0)) for o in obs], dtype=float)
    denom = float(np.sum(w * w))
    beta0 = float(np.sum(proj * w) / denom) if denom > 0 else 0.0

    x0 = np.array([g_l0, g_b0, beta0], dtype=float)

    def obj(x: np.ndarray) -> float:
        g_l = float(x[0]) % 360.0
        g_b = float(np.clip(x[1], -90.0, 90.0))
        beta = float(x[2])
        g = galactic_to_unit(g_l, g_b)
        loss = 0.0
        for i, o in enumerate(obs):
            pred = o.D_kin_vec + (beta * w[i]) * g
            res = o.D_obs_vec - pred
            loss += robust_whitened_loss(res, o.cov_inv, huber_delta)
        return float(loss)

    opt = minimize(obj, x0, method="Nelder-Mead", options={"maxiter": 40000, "xatol": 1e-10, "fatol": 1e-10})
    return {"g_l": float(opt.x[0]) % 360.0, "g_b": float(np.clip(opt.x[1], -90.0, 90.0)), "beta": float(opt.x[2]), "loss": float(opt.fun)}


def fit_o2_single(obs: List[Observation], cfg: Dict[str, Any], bao_init: Dict[str, float]) -> Dict[str, float]:
    H0 = float(cfg["cosmology"]["H0"])
    Om0 = float(cfg["cosmology"]["Om0"])
    r_bao = float(cfg["cosmology"]["r_bao_mpc"])
    p = float(cfg["cosmology"]["bao_p"])

    use_huber = bool(cfg["robust_loss"]["use_huber"])
    huber_delta = float(cfg["robust_loss"]["huber_delta"]) if use_huber else None

    w = np.array([bao_weight(o.z_mean, r_bao_mpc=r_bao, p=p, H0=H0, Om0=Om0) for o in obs], dtype=float)

    g_l0, g_b0, beta0 = float(bao_init["g_l"]), float(bao_init["g_b"]), float(bao_init["beta"])
    g0 = galactic_to_unit(g_l0, g_b0)

    qua = [o for o in obs if o.t_quaia > 0.5]
    if qua:
        mean_q = np.mean([o.D_obs_vec - (o.D_kin_vec + (beta0 * bao_weight(o.z_mean, r_bao, p, H0, Om0)) * g0) for o in qua], axis=0)
        s0 = mean_q / np.linalg.norm(mean_q) if np.linalg.norm(mean_q) > 0 else galactic_to_unit(264.021, 48.253)
    else:
        s0 = galactic_to_unit(264.021, 48.253)

    s_l0 = float(np.degrees(np.arctan2(s0[1], s0[0])) % 360.0)
    s_b0 = float(np.degrees(np.arcsin(np.clip(s0[2], -1.0, 1.0))))

    x0 = np.array([g_l0, g_b0, beta0, s_l0, s_b0, 0.0, 0.0, 0.0], dtype=float)

    def obj(x: np.ndarray) -> float:
        g_l = float(x[0]) % 360.0
        g_b = float(np.clip(x[1], -90.0, 90.0))
        beta = float(x[2])

        s_l = float(x[3]) % 360.0
        s_b = float(np.clip(x[4], -90.0, 90.0))

        a0 = float(x[5])
        a_bmin = float(x[6])
        a_mag = float(x[7])

        g = galactic_to_unit(g_l, g_b)
        s = galactic_to_unit(s_l, s_b)

        loss = 0.0
        for i, o in enumerate(obs):
            A = a0 + a_bmin * ((o.bmin or 30.0) - 30.0) + a_mag * ((o.magcut or 20.0) - 20.0)
            pred = o.D_kin_vec + (beta * w[i]) * g + (o.t_quaia * A) * s
            res = o.D_obs_vec - pred
            loss += robust_whitened_loss(res, o.cov_inv, huber_delta)
        return float(loss)

    opt = minimize(obj, x0, method="Nelder-Mead", options={"maxiter": 80000, "xatol": 1e-10, "fatol": 1e-10})
    return {
        "g_l": float(opt.x[0]) % 360.0,
        "g_b": float(np.clip(opt.x[1], -90.0, 90.0)),
        "beta": float(opt.x[2]),
        "s_l": float(opt.x[3]) % 360.0,
        "s_b": float(np.clip(opt.x[4], -90.0, 90.0)),
        "a0": float(opt.x[5]),
        "a_bmin": float(opt.x[6]),
        "a_mag": float(opt.x[7]),
        "loss": float(opt.fun),
        "success": bool(opt.success),
    }


@dataclass
class Dataset:
    key: str
    variant: str
    obs: List[Observation]


def fit_o2_joint_shared_operator(dsets: List[Dataset], cfg: Dict[str, Any], init_s: Tuple[float, float] = (20.0, 36.0)) -> Dict[str, Any]:
    H0 = float(cfg["cosmology"]["H0"])
    Om0 = float(cfg["cosmology"]["Om0"])
    r_bao = float(cfg["cosmology"]["r_bao_mpc"])
    p = float(cfg["cosmology"]["bao_p"])

    use_huber = bool(cfg["robust_loss"]["use_huber"])
    huber_delta = float(cfg["robust_loss"]["huber_delta"]) if use_huber else None

    ws: List[np.ndarray] = []
    for ds in dsets:
        w = np.array([bao_weight(o.z_mean, r_bao_mpc=r_bao, p=p, H0=H0, Om0=Om0) for o in ds.obs], dtype=float)
        ws.append(w)

    x0_parts: List[float] = []
    for ds in dsets:
        bao = fit_bao_p1(ds.obs, cfg)
        x0_parts.extend([bao["g_l"], bao["g_b"], bao["beta"]])

    s_l0, s_b0 = init_s
    x0_parts.extend([float(s_l0), float(s_b0), 0.0, 0.0, 0.0])
    x0 = np.array(x0_parts, dtype=float)

    n_ds = len(dsets)

    def obj(x: np.ndarray) -> float:
        g_params = x[: 3 * n_ds].reshape((n_ds, 3))
        s_l = float(x[3 * n_ds + 0]) % 360.0
        s_b = float(np.clip(x[3 * n_ds + 1], -90.0, 90.0))
        a0 = float(x[3 * n_ds + 2])
        a_bmin = float(x[3 * n_ds + 3])
        a_mag = float(x[3 * n_ds + 4])

        s = galactic_to_unit(s_l, s_b)

        loss = 0.0
        for j, ds in enumerate(dsets):
            g_l = float(g_params[j, 0]) % 360.0
            g_b = float(np.clip(g_params[j, 1], -90.0, 90.0))
            beta = float(g_params[j, 2])
            g = galactic_to_unit(g_l, g_b)

            wj = ws[j]
            for i, o in enumerate(ds.obs):
                A = a0 + a_bmin * ((o.bmin or 30.0) - 30.0) + a_mag * ((o.magcut or 20.0) - 20.0)
                pred = o.D_kin_vec + (beta * wj[i]) * g + (o.t_quaia * A) * s
                res = o.D_obs_vec - pred
                loss += robust_whitened_loss(res, o.cov_inv, huber_delta)

        return float(loss)

    opt = minimize(obj, x0, method="Nelder-Mead", options={"maxiter": 140000, "xatol": 1e-10, "fatol": 1e-10})

    x = opt.x
    g_params = x[: 3 * n_ds].reshape((n_ds, 3))
    s_l = float(x[3 * n_ds + 0]) % 360.0
    s_b = float(np.clip(x[3 * n_ds + 1], -90.0, 90.0))

    return {
        "g_params": g_params,
        "s_l": s_l,
        "s_b": s_b,
        "a0": float(x[3 * n_ds + 2]),
        "a_bmin": float(x[3 * n_ds + 3]),
        "a_mag": float(x[3 * n_ds + 4]),
        "loss": float(opt.fun),
        "success": bool(opt.success),
    }


def main(argv: List[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="Path to configs/defaults.yaml")
    ap.add_argument("--out_dir", required=True, help="Output directory")
    ap.add_argument("--packs", nargs="*", default=["03_BLM_QG20_CW_ZR", "06_WGV_QG20_CW_ZR"], help="Pack IDs (filename stems in data/packs)")
    args = ap.parse_args(argv)

    cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    l_cmb = float(cfg["cosmology"]["cmb_dipole_galactic_deg"]["l"])
    b_cmb = float(cfg["cosmology"]["cmb_dipole_galactic_deg"]["b"])
    cmb_u = galactic_to_unit(l_cmb, b_cmb)

    def dataset_key(pack_id: str) -> str:
        if pack_id.startswith("03_") or "BLM" in pack_id:
            return "boehme"
        if pack_id.startswith("06_") or "WGV" in pack_id:
            return "wagenveld"
        return pack_id

    single_rows: List[Dict[str, Any]] = []
    dsets_all4: List[Dataset] = []

    for pack_id in args.packs:
        pack_csv = Path("data") / "packs" / f"{pack_id}.csv"
        if not pack_csv.exists():
            raise SystemExit(f"Pack not found: {pack_csv}")

        key = dataset_key(pack_id)

        for variant in ("baseline", "2MRS_clean"):
            rows = load_pack_rows(pack_csv, variant=variant)
            seed = 12345 + (abs(hash((pack_id, variant))) % 100000)
            obs = build_observations(rows, cfg, seed=seed)

            bao = fit_bao_p1(obs, cfg)
            o2 = fit_o2_single(obs, cfg, bao_init=bao)

            n = 3 * len(obs)
            k = 8
            aicc = _aicc(o2["loss"], k=k, n=n)

            g_u = galactic_to_unit(float(o2["g_l"]), float(o2["g_b"]))
            s_u = galactic_to_unit(float(o2["s_l"]), float(o2["s_b"]))

            single_rows.append(
                dict(
                    run="Q",
                    pack=key,
                    variant=variant,
                    model="O2_single",
                    template="",
                    k=float(k),
                    chi2=float(o2["loss"]),
                    aicc=float(aicc),
                    dAICc="",
                    g_l=float(o2["g_l"]),
                    g_b=float(o2["g_b"]),
                    g_to_CMB_deg=float(angle_deg(g_u, cmb_u)),
                    s_l=float(o2["s_l"]),
                    s_b=float(o2["s_b"]),
                    s_to_CMB_deg=float(angle_deg(s_u, cmb_u)),
                )
            )

            dsets_all4.append(Dataset(key=key, variant=variant, obs=obs))

    joint_rows: List[Dict[str, Any]] = []
    for variant in ("baseline", "2MRS_clean"):
        sub = [ds for ds in dsets_all4 if ds.variant == variant and ds.key in ("boehme", "wagenveld")]
        if len(sub) == 2:
            joint = fit_o2_joint_shared_operator(sub, cfg)
            n = 3 * (len(sub[0].obs) + len(sub[1].obs))
            k = 11
            aicc_joint = _aicc(joint["loss"], k=k, n=n)

            sep = sum(float(r["aicc"]) for r in single_rows if r["variant"] == variant and r["pack"] in ("boehme", "wagenveld"))
            dAICc = float(aicc_joint - sep)

            s_u = galactic_to_unit(float(joint["s_l"]), float(joint["s_b"]))

            joint_rows.append(
                dict(
                    run="P",
                    pack="boehme+wagentveld_joint",
                    variant=variant,
                    model="O2_joint_shared_operator",
                    template="",
                    k=float(k),
                    chi2=float(joint["loss"]),
                    aicc=float(aicc_joint),
                    dAICc=float(dAICc),
                    g_l="",
                    g_b="",
                    g_to_CMB_deg="",
                    s_l=float(joint["s_l"]),
                    s_b=float(joint["s_b"]),
                    s_to_CMB_deg=float(angle_deg(s_u, cmb_u)),
                )
            )

    if len(dsets_all4) == 4:
        global_joint = fit_o2_joint_shared_operator(dsets_all4, cfg)
        n = 3 * sum(len(ds.obs) for ds in dsets_all4)
        k = 17
        aicc_global = _aicc(global_joint["loss"], k=k, n=n)

        sep = sum(float(r["aicc"]) for r in single_rows if r["pack"] in ("boehme", "wagenveld"))
        dAICc = float(aicc_global - sep)

        s_u = galactic_to_unit(float(global_joint["s_l"]), float(global_joint["s_b"]))

        joint_rows.append(
            dict(
                run="Q",
                pack="all4_joint",
                variant="baseline+2MRS_clean",
                model="O2_global_joint_shared_operator_all4",
                template="",
                k=float(k),
                chi2=float(global_joint["loss"]),
                aicc=float(aicc_global),
                dAICc=float(dAICc),
                g_l="",
                g_b="",
                g_to_CMB_deg="",
                s_l=float(global_joint["s_l"]),
                s_b=float(global_joint["s_b"]),
                s_to_CMB_deg=float(angle_deg(s_u, cmb_u)),
            )
        )

    out_csv = out_dir / "operator_fits.csv"
    cols = [
        "run", "pack", "variant", "model", "template", "k", "chi2", "aicc", "dAICc",
        "g_l", "g_b", "g_to_CMB_deg", "s_l", "s_b", "s_to_CMB_deg",
    ]

    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in single_rows + joint_rows:
            w.writerow({c: r.get(c, "") for c in cols})

    print(f"[OK] Wrote: {out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
