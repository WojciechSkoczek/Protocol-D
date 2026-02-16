from __future__ import annotations

import math
import numpy as np
from dataclasses import dataclass
from typing import Iterable

# ----------------------------
# Geometry helpers
# ----------------------------
def galactic_to_unit(l_deg: float, b_deg: float) -> np.ndarray:
    l = math.radians(l_deg); b = math.radians(b_deg)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)], dtype=float)

def galactic_to_cartesian(l_deg: float, b_deg: float, r: float) -> np.ndarray:
    return float(r) * galactic_to_unit(l_deg, b_deg)

def cartesian_to_galactic(v: np.ndarray) -> tuple[float,float,float]:
    x,y,z = float(v[0]), float(v[1]), float(v[2])
    r = math.sqrt(x*x+y*y+z*z)
    if r == 0: return 0.0, 0.0, 0.0
    b = math.degrees(math.asin(z/r))
    l = math.degrees(math.atan2(y,x)) % 360.0
    return l,b,r

def angle_deg(u: np.ndarray, v: np.ndarray) -> float:
    uu = u/np.linalg.norm(u); vv=v/np.linalg.norm(v)
    c = float(np.clip(np.dot(uu,vv), -1.0, 1.0))
    return math.degrees(math.acos(c))

def random_cone_perturbation(u0: np.ndarray, sigma_dir_deg: float, n: int, rng: np.random.Generator) -> np.ndarray:
    # sample directions around u0 with small-angle Gaussian in tangent plane
    u0 = u0/np.linalg.norm(u0)
    # build orthonormal basis (e1,e2) perpendicular to u0
    tmp = np.array([1.0,0.0,0.0], dtype=float)
    if abs(np.dot(tmp,u0)) > 0.9:
        tmp = np.array([0.0,1.0,0.0], dtype=float)
    e1 = tmp - np.dot(tmp,u0)*u0
    e1 = e1/np.linalg.norm(e1)
    e2 = np.cross(u0,e1)
    sig = math.radians(max(0.0, float(sigma_dir_deg)))
    a = rng.normal(0.0, sig, size=n)
    b = rng.normal(0.0, sig, size=n)
    u = u0[None,:] + a[:,None]*e1[None,:] + b[:,None]*e2[None,:]
    u = u / np.linalg.norm(u, axis=1)[:,None]
    return u

# ----------------------------
# Cosmology weight (same as harness)
# ----------------------------
def E_z(z: float, Om0: float) -> float:
    return math.sqrt(Om0*(1.0+z)**3 + (1.0-Om0))

def comoving_distance_mpc(z: float, H0: float=70.0, Om0: float=0.3, n_steps: int=4000) -> float:
    if not np.isfinite(z) or z <= 0: return 0.0
    zs = np.linspace(0.0, float(z), int(n_steps)+1)
    invE = 1.0/np.array([E_z(float(zz), Om0) for zz in zs], dtype=float)
    integral = np.trapz(invE, zs)
    c_km_s = 299792.458
    return (c_km_s/H0)*float(integral)

def bao_weight(z_mean: float, r_bao_mpc: float=147.0, p: float=1.0, H0: float=70.0, Om0: float=0.3) -> float:
    chi = comoving_distance_mpc(z_mean, H0=H0, Om0=Om0)
    return (r_bao_mpc/chi)**p if chi>0 else 0.0

# ----------------------------
# Data model
# ----------------------------
@dataclass
class Observation:
    name: str
    l_deg: float
    b_deg: float
    D_obs: float
    sigma_mag: float
    sigma_dir_deg: float
    D_kin: float
    z_mean: float
    t_quaia: float = 0.0
    bmin: float | None = None
    magcut: float | None = None

    D_obs_vec: np.ndarray | None = None
    D_kin_vec: np.ndarray | None = None
    cov_inv: np.ndarray | None = None

    def build(self, rng: np.random.Generator, l_cmb: float, b_cmb: float, n_samples: int, cov_floor: float) -> None:
        self.D_obs_vec = galactic_to_cartesian(self.l_deg, self.b_deg, self.D_obs)
        self.D_kin_vec = galactic_to_cartesian(l_cmb, b_cmb, self.D_kin)
        u0 = galactic_to_unit(self.l_deg, self.b_deg)
        u_samp = random_cone_perturbation(u0, self.sigma_dir_deg, n_samples, rng)
        mag = rng.normal(self.D_obs, self.sigma_mag, size=n_samples)
        mag = np.clip(mag, 0.0, None)
        v = mag[:,None]*u_samp
        mean = np.mean(v, axis=0)
        dv = v-mean[None,:]
        cov = (dv.T@dv)/(n_samples-1)
        cov = cov + np.eye(3)*cov_floor
        self.cov_inv = np.linalg.inv(cov)

def chi2(res: np.ndarray, cov_inv: np.ndarray) -> float:
    return float(res.T @ cov_inv @ res)

def huber(r: float, delta: float) -> float:
    a = abs(r)
    if a <= delta: return 0.5*r*r
    return delta*(a - 0.5*delta)

def robust_whitened_loss(res: np.ndarray, cov_inv: np.ndarray, huber_delta: float | None) -> float:
    # whiten via Cholesky of cov_inv (numerically: use eig)
    if huber_delta is None:
        return chi2(res, cov_inv)
    w, V = np.linalg.eigh(cov_inv)
    w = np.clip(w, 0.0, None)
    Wsqrt = V @ np.diag(np.sqrt(w)) @ V.T
    r = Wsqrt @ res
    return sum(huber(float(x), huber_delta) for x in r)

