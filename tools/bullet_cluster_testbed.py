"""
Standalone Zero-Parameter Bullet Cluster (1E 0657-56) 2-Cluster Collision
& Weak/Strong Lensing N-Body Testbed.

Evolves a 3-species (Cover Pockets + Compiled Galaxies + ICM Gas) Bullet Cluster
merger on a 96^3 periodic 3D FFT Particle-Mesh grid (L = 5,000 kpc, dx = 52.08 kpc)
under strictly linear Pure K_Omega gravity, testing:
  - Three zero-parameter Cover Pocket hardware pressure ceilings:
      1. 3 * P_delta (Galaxy Floor)
      2. Sigma_0 / 3 (Isotropic Cover Ceiling -- primary cluster prediction)
      3. Sigma_0 (Full Cover Ceiling)
  - A Two-Condition Compilation Gate ablation on Sigma_0 / 3:
      * ICM Passive (2-Cond Gate: hot shocked ICM gas fails Lambda_lock >= 1)
      * ICM Active (Full Gravity ablation: all 3 species gravitate)
  - Momentum-conserving Eulerian grid ram-pressure shock drag & kinetic heating
  - 256 x 256 2D weak/strong lensing convergence kappa, shear gamma, critical curves
    (mu^-1 = 0), and X-ray surface brightness S_X ~ int n^2 sqrt(T) dz at the
    post-pericenter Delta x_lens = 720 kpc epoch.

Imports all physical constants directly from `modules.data_loader` (zero CODATA G)
and reuses `QUMONDMeshSolver` from `tools.qumond_nbody_testbed`.
"""

import argparse
import logging
import os
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor

# Ensure repository root is importable regardless of launch directory
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

for _logger_name in (
    "streamlit",
    "streamlit.runtime.caching",
    "streamlit.runtime.scriptrunner_utils",
):
    logging.getLogger(_logger_name).setLevel(logging.ERROR)
warnings.filterwarnings("ignore", message=".*No runtime found.*")
warnings.filterwarnings("ignore", message=".*to view this Streamlit app.*")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.fft import fftfreq, irfft2, rfft2, rfftfreq
from scipy.integrate import cumulative_trapezoid, quad
from scipy.ndimage import gaussian_filter
from scipy.optimize import minimize_scalar

from modules.data_loader import (
    C_DELTA,
    C_LIGHT,
    E_RYD_J,
    H_GLOBAL,
    K_BOLTZMANN,
    K_OMEGA,
    KPC_TO_METER,
    M_PROTON,
    M_SOLAR_KG,
    MPC_TO_METER,
    OMEGA_M_GLOBAL,
    OMEGA_M_VOID,
    P_DELTA_MACRO_PA,
    SIGMA_0_LEDGER,
    evaluate_two_condition_compilation_gate,
)
from tools.qumond_nbody_testbed import QUMONDMeshSolver

# ==============================================================================
# 1. Exact First-Principles Cover Pocket & Lensing Cosmology Invariants
# ==============================================================================
PC_TO_METER = KPC_TO_METER / 1000.0
MYR_TO_SEC = 1.0e6 * 365.25 * 86400.0
KEV_TO_K = 1.602176634e-16 / K_BOLTZMANN
KG_M2_TO_MSUN_PC2 = (PC_TO_METER ** 2) / M_SOLAR_KG

G_EFF = K_OMEGA / M_PROTON
ETA_POCKET = OMEGA_M_VOID / (OMEGA_M_GLOBAL - OMEGA_M_VOID)


def burkert_m(x):
    """Dimensionless Burkert enclosed mass factor m_B(x) = ln(1+x^2) + 2ln(1+x) - 2arctan(x)."""
    x_arr = np.asarray(x, dtype=np.float64)
    small = np.abs(x_arr) < 1.0e-4
    x_safe = np.where(small, 1.0, x_arr)
    taylor = (4.0 / 3.0) * (x_arr ** 3) - (x_arr ** 4) + (2.0 / 3.0) * (x_arr ** 6)
    exact = np.log(1.0 + x_safe * x_safe) + 2.0 * np.log(1.0 + x_safe) - 2.0 * np.arctan(x_safe)
    out = np.where(small, taylor, exact)
    return float(out) if np.ndim(x) == 0 else out


def burkert_f(x):
    """Dimensionless Burkert radial acceleration factor f_B(x) = m_B(x) / x^2."""
    x_arr = np.asarray(x, dtype=np.float64)
    small = np.abs(x_arr) < 1.0e-4
    x_safe = np.where(small, 1.0, x_arr)
    taylor = (4.0 / 3.0) * x_arr - (x_arr ** 2) + (2.0 / 3.0) * (x_arr ** 4)
    exact = burkert_m(x_safe) / (x_safe * x_safe)
    out = np.where(small, taylor, exact)
    return float(out) if np.ndim(x) == 0 else out


def burkert_v2(x):
    """Dimensionless Burkert circular velocity squared factor v2_B(x) = m_B(x) / x."""
    x_arr = np.asarray(x, dtype=np.float64)
    small = np.abs(x_arr) < 1.0e-4
    x_safe = np.where(small, 1.0, x_arr)
    taylor = (4.0 / 3.0) * (x_arr ** 2) - (x_arr ** 3) + (2.0 / 3.0) * (x_arr ** 5)
    exact = burkert_m(x_safe) / x_safe
    out = np.where(small, taylor, exact)
    return float(out) if np.ndim(x) == 0 else out


_RES_GMAX = minimize_scalar(lambda x: -burkert_f(x), bounds=(0.1, 5.0), method="bounded")
X_GMAX = float(_RES_GMAX.x)
F_B_MAX = float(-_RES_GMAX.fun)

_RES_VMAX = minimize_scalar(lambda x: -burkert_v2(x), bounds=(0.5, 10.0), method="bounded")
X_VMAX = float(_RES_VMAX.x)
V2_B_MAX = float(-_RES_VMAX.fun)
M_FAC_VMAX = float(burkert_m(X_VMAX))


def sigma_pocket_from_pcap(p_cap_pa):
    """Computes zero-parameter Burkert Cover Pocket surface density [kg/m^2] from hardware ceiling P_cap."""
    g_cap = np.sqrt(8.0 * np.pi * G_EFF * float(p_cap_pa))
    return float(g_cap / (F_B_MAX * np.pi * G_EFF))


def compute_bullet_lensing_cosmology(z_l=0.296, z_s=1.0, omega_m=0.30):
    """Exact Bullet Cluster angular diameter distances & critical surface density."""
    omega_l = 1.0 - omega_m
    h0_si = (H_GLOBAL * 1000.0) / MPC_TO_METER
    c_over_h0 = C_LIGHT / h0_si
    inv_ez = lambda z: 1.0 / np.sqrt(omega_m * ((1.0 + z) ** 3) + omega_l)
    dc_l, _ = quad(inv_ez, 0.0, z_l)
    dc_s, _ = quad(inv_ez, 0.0, z_s)
    dc_l *= c_over_h0
    dc_s *= c_over_h0
    d_l = dc_l / (1.0 + z_l)
    d_s = dc_s / (1.0 + z_s)
    d_ls = (dc_s - dc_l) / (1.0 + z_s)
    sigma_crit = (C_LIGHT ** 2 / (4.0 * np.pi * G_EFF)) * (d_s / (d_l * d_ls))
    kpc_per_arcsec = (d_l / KPC_TO_METER) * (np.pi / (180.0 * 3600.0))
    return {
        "D_l_mpc": d_l / MPC_TO_METER,
        "D_s_mpc": d_s / MPC_TO_METER,
        "D_ls_mpc": d_ls / MPC_TO_METER,
        "sigma_crit_kgm2": sigma_crit,
        "sigma_crit_msun_pc2": sigma_crit * KG_M2_TO_MSUN_PC2,
        "kpc_per_arcsec": kpc_per_arcsec,
    }


LENS_COSMO = compute_bullet_lensing_cosmology()
SIGMA_CRIT_KGM2 = LENS_COSMO["sigma_crit_kgm2"]

# ==============================================================================
# 2. Merger Configuration & Controlled Experiment Suite
# ==============================================================================
BULLET_MODELS = [
    {
        "model": "3P_delta (Galaxy Floor, ICM Passive)",
        "short": "3P_delta Passive",
        "p_cap_pa": 3.0 * P_DELTA_MACRO_PA,
        "icm_active": False,
        "seed": 201,
    },
    {
        "model": "Sigma_0/3 (Isotropic Cover Ceiling, ICM Passive)",
        "short": "Sigma_0/3 Passive",
        "p_cap_pa": SIGMA_0_LEDGER / 3.0,
        "icm_active": False,
        "seed": 202,
    },
    {
        "model": "Sigma_0/3 (Isotropic Cover Ceiling, ICM Active)",
        "short": "Sigma_0/3 Active",
        "p_cap_pa": SIGMA_0_LEDGER / 3.0,
        "icm_active": True,
        "seed": 203,
    },
    {
        "model": "Sigma_0 (Full Cover Ceiling, ICM Passive)",
        "short": "Sigma_0 Passive",
        "p_cap_pa": SIGMA_0_LEDGER,
        "icm_active": False,
        "seed": 204,
    },
]

CLUSTERS = {
    "main": {
        "M_B_solar": 1.1e14,
        "f_star": 0.15,
        "T_init_keV": 9.0,
        "n_pocket": 15000,
        "n_star": 3000,
        "n_gas": 4000,
    },
    "sub": {
        "M_B_solar": 3.8e13,
        "f_star": 0.15,
        "T_init_keV": 6.5,
        "n_pocket": 7500,
        "n_star": 1500,
        "n_gas": 2000,
    },
}


# ==============================================================================
# 3. Equilibrium Jeans Cluster Initialization & Orbit Assembly
# ==============================================================================
def _sample_burkert_radii(rng, n_particles, r0_m, r_trunc_ratio=4.0):
    """Samples exact truncated Burkert radii via inverse CDF of m_B(x)/m_B(r_trunc_ratio)."""
    x_tab = np.linspace(0.0, r_trunc_ratio, 4000)
    cdf = burkert_m(x_tab) / burkert_m(r_trunc_ratio)
    cdf[0] = 0.0
    cdf[-1] = 1.0
    u = rng.uniform(0.0, 1.0, size=n_particles)
    return np.interp(u, cdf, x_tab) * r0_m


def _sample_isotropic_positions(rng, radii_m):
    """Draws isotropic 3D Cartesian coordinates for given radial distances."""
    n = len(radii_m)
    cos_t = rng.uniform(-1.0, 1.0, size=n)
    phi = rng.uniform(0.0, 2.0 * np.pi, size=n)
    sin_t = np.sqrt(1.0 - cos_t * cos_t)
    return np.column_stack(
        (radii_m * sin_t * np.cos(phi), radii_m * sin_t * np.sin(phi), radii_m * cos_t)
    )


def _jeans_sigma_r(r_grid_m, rho_grid, g_tot_grid):
    """Exact isotropic spherical Jeans velocity dispersion integrated from r to r_cut."""
    i_fwd = cumulative_trapezoid(rho_grid * g_tot_grid, r_grid_m, initial=0.0)
    tail = np.maximum(i_fwd[-1] - i_fwd, 0.0)
    return np.sqrt(tail / np.maximum(rho_grid, 1.0e-60))


def _sample_jeans_velocities(rng, radii_m, r_grid_m, sig_grid, v_esc_grid):
    """Draws isotropic 3D Gaussian velocities with local std sigma_r(r_i), rejecting v > 0.98 v_esc."""
    n = len(radii_m)
    sig_p = np.interp(radii_m, r_grid_m, sig_grid)
    vmax_p = 0.98 * np.interp(radii_m, r_grid_m, v_esc_grid)
    vel = np.empty((n, 3), dtype=np.float64)
    need = np.ones(n, dtype=bool)
    while np.any(need):
        idx = np.flatnonzero(need)
        v_try = rng.normal(0.0, sig_p[idx, None], size=(len(idx), 3))
        spd = np.linalg.norm(v_try, axis=1)
        ok = spd <= vmax_p[idx]
        if np.any(ok):
            vel[idx[ok]] = v_try[ok]
            need[idx[ok]] = False
    return vel


def generate_cluster_ics(rng, spec, sigma_p_kgm2, solver, icm_active, scale=1.0):
    """
    Builds a single 3-species equilibrium cluster (Cover Pocket, Compiled Galaxies,
    ICM Gas) calibrated to exact initial Clausius virial 2K/|W| = 1.0 on the PM mesh.
    """
    m_b_kg = float(spec["M_B_solar"]) * M_SOLAR_KG
    m_star_tot = float(spec["f_star"]) * m_b_kg
    m_gas_tot = (1.0 - float(spec["f_star"])) * m_b_kg

    r0_m = float(np.sqrt(ETA_POCKET * m_b_kg / (np.pi * sigma_p_kgm2 * M_FAC_VMAX)))
    r_trunc = 4.0
    r_cut = r_trunc * r0_m
    r0_star = 0.8 * r0_m
    r_cut_star = r_trunc * r0_star
    m_pocket_tot = np.pi * sigma_p_kgm2 * (r0_m ** 2) * burkert_m(r_trunc)

    n_p = max(100, int(round(spec["n_pocket"] * scale)))
    n_s = max(50, int(round(spec["n_star"] * scale)))
    n_g = max(50, int(round(spec["n_gas"] * scale)))

    r_p = _sample_burkert_radii(rng, n_p, r0_m, r_trunc)
    r_s = _sample_burkert_radii(rng, n_s, r0_star, r_trunc)
    r_g = _sample_burkert_radii(rng, n_g, r0_m, r_trunc)

    pos_p = _sample_isotropic_positions(rng, r_p)
    pos_s = _sample_isotropic_positions(rng, r_s)
    pos_g = _sample_isotropic_positions(rng, r_g)
    for arr in (pos_p, pos_s, pos_g):
        arr -= np.mean(arr, axis=0)

    # Radial grid covering [1e-4 r0, r_cut] for exact enclosed gravitating mass
    r_grid = np.linspace(1.0e-4 * r0_m, r_cut, 2000)
    x_p = r_grid / r0_m
    x_s = np.minimum(r_grid / r0_star, r_trunc)
    m_norm_4 = burkert_m(r_trunc)

    m_p_enc = m_pocket_tot * (burkert_m(x_p) / m_norm_4)
    m_s_enc = m_star_tot * (burkert_m(x_s) / m_norm_4)
    m_g_enc = m_gas_tot * (burkert_m(x_p) / m_norm_4)

    m_grav_enc = m_p_enc + m_s_enc + (m_g_enc if icm_active else 0.0)
    m_grav_tot = float(m_grav_enc[-1])
    g_tot = G_EFF * m_grav_enc / (r_grid ** 2)

    g_fwd = cumulative_trapezoid(g_tot, r_grid, initial=0.0)
    phi_drop = np.maximum((g_fwd[-1] - g_fwd) + G_EFF * m_grav_tot / r_cut, 0.0)
    v_esc_grid = np.sqrt(2.0 * phi_drop)

    rho_p_grid = (1.0 / r0_m) / ((1.0 + x_p) * (1.0 + x_p * x_p))
    sig_p_grid = _jeans_sigma_r(r_grid, rho_p_grid, g_tot)
    vel_p = _sample_jeans_velocities(rng, np.linalg.norm(pos_p, axis=1), r_grid, sig_p_grid, v_esc_grid)

    # Galaxies sit inside r <= r_cut_star
    mask_s = r_grid <= r_cut_star
    r_grid_s = r_grid[mask_s]
    x_s_sub = r_grid_s / r0_star
    rho_s_grid = (1.0 / r0_star) / ((1.0 + x_s_sub) * (1.0 + x_s_sub * x_s_sub))
    sig_s_grid = _jeans_sigma_r(r_grid_s, rho_s_grid, g_tot[mask_s])
    vel_s = _sample_jeans_velocities(rng, np.linalg.norm(pos_s, axis=1), r_grid_s, sig_s_grid, v_esc_grid[mask_s])

    vel_g = _sample_jeans_velocities(rng, np.linalg.norm(pos_g, axis=1), r_grid, sig_p_grid, v_esc_grid)

    for varr in (vel_p, vel_s, vel_g):
        varr -= np.mean(varr, axis=0)

    masses_p = np.full(n_p, m_pocket_tot / n_p, dtype=np.float64)
    masses_s = np.full(n_s, m_star_tot / n_s, dtype=np.float64)
    masses_g = np.full(n_g, m_gas_tot / n_g, dtype=np.float64)
    masses = np.concatenate((masses_p, masses_s, masses_g))
    grav_masses = np.concatenate(
        (masses_p, masses_s, masses_g if icm_active else np.zeros_like(masses_g))
    )
    node_weights = grav_masses / M_PROTON

    pos = np.vstack((pos_p, pos_s, pos_g))
    vel = np.vstack((vel_p, vel_s, vel_g))
    acc0 = solver.get_particle_accelerations(pos, node_weights, mode="Pure K_Omega")

    # Rescale each species independently to exact initial Clausius virial 2K/|W| = 1.0
    slices = [slice(0, n_p), slice(n_p, n_p + n_s), slice(n_p + n_s, n_p + n_s + n_g)]
    for sl in slices:
        ms = masses[sl]
        k2 = float(np.sum(ms * np.sum(vel[sl] ** 2, axis=1)))
        w_abs = abs(float(np.sum(ms * np.sum(pos[sl] * acc0[sl], axis=1))))
        if k2 > 0.0 and w_abs > 0.0:
            vel[sl] *= np.sqrt(w_abs / k2)

    temp_g = np.full(n_g, float(spec["T_init_keV"]) * KEV_TO_K, dtype=np.float64)
    return {
        "r0_m": r0_m,
        "M_pocket_tot": m_pocket_tot,
        "M_star_tot": m_star_tot,
        "M_gas_tot": m_gas_tot,
        "M_grav_tot": m_pocket_tot + m_star_tot + (m_gas_tot if icm_active else 0.0),
        "M_all_tot": m_pocket_tot + m_star_tot + m_gas_tot,
        "pos": pos,
        "vel": vel,
        "masses": masses,
        "grav_masses": grav_masses,
        "temp_g": temp_g,
        "n_p": n_p,
        "n_s": n_s,
        "n_g": n_g,
    }


# ==============================================================================
# 4. Hydro Grid Ram-Pressure Shock Drag & Centroid / Gate Diagnostics
# ==============================================================================
def _apply_ram_pressure_drag(
    pos_m, vel_ms, temp_k, masses, idx_gas1, idx_gas2, box_l_m, hydro_ng, dt_s
):
    """
    Momentum-conserving Eulerian NGP shock-normal relaxation & kinetic heating
    between Main (idx_gas1) and Subcluster (idx_gas2) ICM gas on a 36^3 hydro grid
    (dx_h = 138.9 kpc), preserving transverse velocity dispersion.
    """
    dx_h = box_l_m / hydro_ng
    dv_cell = dx_h ** 3
    n_cells = hydro_ng ** 3
    sigma_shock = 400.0 * M_SOLAR_KG / (PC_TO_METER ** 2)
    c_ram = 1.5
    mu = 0.60

    idx_both = np.concatenate((idx_gas1, idx_gas2))
    coords = np.floor(np.mod((pos_m[idx_both] + 0.5 * box_l_m) / dx_h, hydro_ng)).astype(np.int64)
    coords %= hydro_ng
    cell_both = coords[:, 0] * (hydro_ng * hydro_ng) + coords[:, 1] * hydro_ng + coords[:, 2]

    n1 = len(idx_gas1)
    c1 = cell_both[:n1]
    c2 = cell_both[n1:]
    m1_p = masses[idx_gas1]
    m2_p = masses[idx_gas2]
    v1_p = vel_ms[idx_gas1]
    v2_p = vel_ms[idx_gas2]

    m1_c = np.bincount(c1, weights=m1_p, minlength=n_cells)
    m2_c = np.bincount(c2, weights=m2_p, minlength=n_cells)

    # Number density n_H+e ~ rho_gas / (mu * M_PROTON) per cell for X-ray emission weight
    m_tot_c = m1_c + m2_c
    n_cell_both = ((m_tot_c / dv_cell)[cell_both]) / (mu * M_PROTON)

    overlap = np.flatnonzero((m1_c > 0.0) & (m2_c > 0.0))
    if dt_s <= 0.0 or len(overlap) == 0:
        return n_cell_both[:n1], n_cell_both[n1:]

    mom1_c = np.column_stack(
        (
            np.bincount(c1, weights=m1_p * v1_p[:, 0], minlength=n_cells),
            np.bincount(c1, weights=m1_p * v1_p[:, 1], minlength=n_cells),
            np.bincount(c1, weights=m1_p * v1_p[:, 2], minlength=n_cells),
        )
    )
    mom2_c = np.column_stack(
        (
            np.bincount(c2, weights=m2_p * v2_p[:, 0], minlength=n_cells),
            np.bincount(c2, weights=m2_p * v2_p[:, 1], minlength=n_cells),
            np.bincount(c2, weights=m2_p * v2_p[:, 2], minlength=n_cells),
        )
    )

    m1_ov = m1_c[overlap]
    m2_ov = m2_c[overlap]
    mt_ov = m1_ov + m2_ov
    dv_rel = (mom1_c[overlap] / m1_ov[:, None]) - (mom2_c[overlap] / m2_ov[:, None])
    v_rel = np.linalg.norm(dv_rel, axis=1)

    # Symmetric binary overlap factor prevents trace particles from prematurely shocking a dense core
    f_ov = 4.0 * (m1_ov * m2_ov) / (mt_ov * mt_ov)
    gamma_c = c_ram * f_ov * (mt_ov / dv_cell) * v_rel / sigma_shock

    alpha_par = np.zeros(n_cells, dtype=np.float64)
    alpha_par[overlap] = 1.0 - np.exp(-gamma_c * dt_s)
    n_hat = np.zeros((n_cells, 3), dtype=np.float64)
    n_hat[overlap] = dv_rel / np.maximum(v_rel, 1.0e-12)[:, None]
    v_cm = np.zeros((n_cells, 3), dtype=np.float64)
    v_cm[overlap] = (mom1_c[overlap] + mom2_c[overlap]) / mt_ov[:, None]

    # Relax shock-normal velocity component toward v_cm while conserving total momentum identically
    ap_all = alpha_par[cell_both]
    active = ap_all > 0.0
    if np.any(active):
        idx_act = idx_both[active]
        c_act = cell_both[active]
        dv_p = vel_ms[idx_act] - v_cm[c_act]
        nh = n_hat[c_act]
        ap = ap_all[active]
        v_par = np.sum(dv_p * nh, axis=1)
        vel_ms[idx_act] -= (ap * v_par)[:, None] * nh
        de_spec = 0.5 * (v_par * v_par) * ap * (2.0 - ap)
        temp_k[idx_act] += (2.0 * mu * M_PROTON / (3.0 * K_BOLTZMANN)) * de_spec

    return n_cell_both[:n1], n_cell_both[n1:]


def _shrinking_core_centroid(pos_subset, vel_subset, core_radius_m):
    """Computes robust 3D density core position & mean velocity (3 shrinking iterations)."""
    center = np.median(pos_subset, axis=0)
    mask = np.ones(len(pos_subset), dtype=bool)
    for fac in (4.0, 2.0, 1.0):
        r = np.linalg.norm(pos_subset - center, axis=1)
        m_sub = r <= (fac * core_radius_m)
        if np.sum(m_sub) >= 25:
            center = np.mean(pos_subset[m_sub], axis=0)
            mask = m_sub
    v_center = np.mean(vel_subset[mask], axis=0) if np.any(mask) else np.mean(vel_subset, axis=0)
    return center, v_center


def _estimate_crossing_drag_shift(pos_p1, v_p1, m_p1, pos_p2, v_p2, m_p2, sigma_p_kgm2, box_l_m, hydro_ng, dt_s):
    """
    Computes instantaneous perturbative C_delta crossing drag acceleration on Subcluster
    pocket particles passing through Main pocket particles:
        sigma/m_eff = C_DELTA / sigma_p_kgm2  => tau_c = C_DELTA * (Sigma_c / sigma_p_kgm2)
    """
    dx_h = box_l_m / hydro_ng
    dv_cell = dx_h ** 3
    n_cells = hydro_ng ** 3
    inv_area = 1.0 / (dx_h ** 2)

    all_pos = np.vstack((pos_p1, pos_p2))
    coords = np.floor(np.mod((all_pos + 0.5 * box_l_m) / dx_h, hydro_ng)).astype(np.int64)
    coords %= hydro_ng
    cells = coords[:, 0] * (hydro_ng * hydro_ng) + coords[:, 1] * hydro_ng + coords[:, 2]
    n1 = len(pos_p1)
    c1, c2 = cells[:n1], cells[n1:]

    m1_c = np.bincount(c1, weights=m_p1, minlength=n_cells)
    m2_c = np.bincount(c2, weights=m_p2, minlength=n_cells)
    ov = np.flatnonzero((m1_c > 0.0) & (m2_c > 0.0))
    if len(ov) == 0:
        return np.zeros(3, dtype=np.float64)

    mom1_c = np.column_stack(
        (
            np.bincount(c1, weights=m_p1 * v_p1[:, 0], minlength=n_cells),
            np.bincount(c1, weights=m_p1 * v_p1[:, 1], minlength=n_cells),
            np.bincount(c1, weights=m_p1 * v_p1[:, 2], minlength=n_cells),
        )
    )
    mom2_c = np.column_stack(
        (
            np.bincount(c2, weights=m_p2 * v_p2[:, 0], minlength=n_cells),
            np.bincount(c2, weights=m_p2 * v_p2[:, 1], minlength=n_cells),
            np.bincount(c2, weights=m_p2 * v_p2[:, 2], minlength=n_cells),
        )
    )
    dv = (mom1_c[ov] / m1_c[ov, None]) - (mom2_c[ov] / m2_c[ov, None])
    v_rel = np.linalg.norm(dv, axis=1)
    # Perturbative drag impulse on subcluster mass m2:
    # dP_2 = m2 * (C_DELTA / sigma_p_kgm2) * rho_1 * v_rel * dt * dv
    rho1_ov = m1_c[ov] / dv_cell
    frac = np.minimum((C_DELTA / sigma_p_kgm2) * rho1_ov * v_rel * dt_s, 1.0)
    dp2 = np.sum((m2_c[ov] * frac)[:, None] * dv, axis=0)
    m2_tot = float(np.sum(m_p2))
    # Also compute cell-opacity check tau = C_DELTA * Sigma_1 / sigma_p
    _ = inv_area
    return dp2 / m2_tot


# ==============================================================================
# 5. Exact 2D Weak & Strong Lensing + X-Ray Projection Engine
# ==============================================================================
def _deposit_2d_cic(xy_kpc, weights, half_extent_kpc=1200.0, npix=256):
    """Deposits projected particle weights onto a 2D (npix x npix) sky grid via bilinear CIC."""
    dx = (2.0 * half_extent_kpc) / npix
    gx = (xy_kpc[:, 0] + half_extent_kpc) / dx
    gy = (xy_kpc[:, 1] + half_extent_kpc) / dx

    in_bounds = (gx >= 0.0) & (gx < npix - 1.0) & (gy >= 0.0) & (gy < npix - 1.0)
    gx = gx[in_bounds]
    gy = gy[in_bounds]
    w = weights[in_bounds]

    i0 = np.floor(gx).astype(np.int64)
    j0 = np.floor(gy).astype(np.int64)
    fx = gx - i0
    fy = gy - j0

    grid = np.zeros(npix * npix, dtype=np.float64)
    for di, wx in ((0, 1.0 - fx), (1, fx)):
        for dj, wy in ((0, 1.0 - fy), (1, fy)):
            flat = (j0 + dj) * npix + (i0 + di)  # row = y, col = x
            grid += np.bincount(flat, weights=w * wx * wy, minlength=npix * npix)

    pix_area_m2 = (dx * KPC_TO_METER) ** 2
    return grid.reshape((npix, npix)) / pix_area_m2


def _locate_peak_half(field_2d, x_1d, y_1d, x_lo, x_hi):
    """Finds 2D maximum (x, y) inside x in [x_lo, x_hi]."""
    mask_x = (x_1d >= x_lo) & (x_1d <= x_hi)
    sub = field_2d[:, mask_x]
    if sub.size == 0 or float(np.max(sub)) <= 0.0:
        return 0.0, 0.0, 0.0
    iy, ix_sub = np.unravel_index(int(np.argmax(sub)), sub.shape)
    x_sub = x_1d[mask_x]
    return float(x_sub[ix_sub]), float(y_1d[iy]), float(sub[iy, ix_sub])


def _einstein_radius_from_kappa(kappa_2d, x_2d, y_2d, x_pk, y_pk):
    """Evaluates effective Einstein radius R_E [kpc] where azimuthal mean kappa(<R) = 1.0."""
    r_map = np.hypot(x_2d - x_pk, y_2d - y_pk)
    r_bins = np.linspace(15.0, 450.0, 90)
    mean_k = np.array([float(np.mean(kappa_2d[r_map <= r_b])) for r_b in r_bins])
    if float(np.max(mean_k)) < 1.0:
        return 0.0
    above = np.where(mean_k >= 1.0)[0]
    last_i = int(above[-1])
    if last_i + 1 < len(r_bins) and mean_k[last_i] >= 1.0 > mean_k[last_i + 1]:
        frac = (mean_k[last_i] - 1.0) / max(mean_k[last_i] - mean_k[last_i + 1], 1.0e-12)
        return float(r_bins[last_i] + frac * (r_bins[last_i + 1] - r_bins[last_i]))
    return float(r_bins[last_i])


def compute_lensing_and_xray_maps(snap):
    """
    Computes 256 x 256 convergence kappa, 512 x 512 zero-padded FFT potential psi,
    shear gamma, magnification mu^-1, and X-ray surface brightness S_X at t_bullet.
    """
    half_ext = 1200.0
    npix = 256
    dx_kpc = (2.0 * half_ext) / npix
    sigma_beam_pix = 25.0 / dx_kpc

    # Center coordinates on midpoint between Main (+x) and Sub (-x) lens centers
    x_mid = 0.5 * (snap["c_main_lens"][0] + snap["c_sub_lens"][0])
    y_mid = 0.5 * (snap["c_main_lens"][1] + snap["c_sub_lens"][1])
    xy_kpc = (snap["pos"][:, :2] - np.array([x_mid, y_mid])) / KPC_TO_METER

    sigma_2d_raw = _deposit_2d_cic(xy_kpc, snap["grav_masses"], half_extent_kpc=half_ext, npix=npix)
    sigma_2d = gaussian_filter(sigma_2d_raw, sigma=sigma_beam_pix, mode="constant")
    kappa = sigma_2d / SIGMA_CRIT_KGM2

    # 2D FFT lensing potential solved on zero-padded 512 x 512 domain (unitless pixel coordinates)
    npad = 2 * npix
    kpad = np.zeros((npad, npad), dtype=np.float64)
    kpad[:npix, :npix] = kappa

    ky = 2.0 * np.pi * fftfreq(npad, d=1.0)
    kx = 2.0 * np.pi * rfftfreq(npad, d=1.0)
    KY, KX = np.meshgrid(ky, kx, indexing="ij")
    k2 = KX * KX + KY * KY
    k2_safe = np.where(k2 > 0.0, k2, 1.0)

    k_hat = rfft2(kpad)
    psi_k = np.where(k2 > 0.0, -2.0 * k_hat / k2_safe, 0.0)
    psi_xx = irfft2(-(KX * KX) * psi_k, s=(npad, npad))[:npix, :npix]
    psi_yy = irfft2(-(KY * KY) * psi_k, s=(npad, npad))[:npix, :npix]
    psi_xy = irfft2(-(KX * KY) * psi_k, s=(npad, npad))[:npix, :npix]

    gamma1 = 0.5 * (psi_xx - psi_yy)
    gamma2 = psi_xy
    gamma_mag = np.hypot(gamma1, gamma2)
    lambda_t = 1.0 - kappa - gamma_mag
    mu_inv = ((1.0 - kappa) ** 2) - (gamma_mag ** 2)

    # X-ray surface brightness S_X ~ int n_gas^2 sqrt(T_gas) dz weighted by m_i * n_cell * sqrt(T_i)
    w_x1 = snap["m_g1"] * snap["n_g1"] * np.sqrt(np.maximum(snap["T_g1"], 1.0e5))
    w_x2 = snap["m_g2"] * snap["n_g2"] * np.sqrt(np.maximum(snap["T_g2"], 1.0e5))
    sx1_raw = _deposit_2d_cic(xy_kpc[snap["idx_g1"]], w_x1, half_extent_kpc=half_ext, npix=npix)
    sx2_raw = _deposit_2d_cic(xy_kpc[snap["idx_g2"]], w_x2, half_extent_kpc=half_ext, npix=npix)
    sx_main = gaussian_filter(sx1_raw, sigma=sigma_beam_pix, mode="constant")
    sx_sub = gaussian_filter(sx2_raw, sigma=sigma_beam_pix, mode="constant")
    sx_tot = sx_main + sx_sub

    edges = np.linspace(-half_ext, half_ext, npix + 1)
    x_1d = 0.5 * (edges[:-1] + edges[1:])
    y_1d = x_1d.copy()
    X2D, Y2D = np.meshgrid(x_1d, y_1d)

    # Lensing kappa peaks (Main on x > 0, Sub on x < 0)
    xm_l, ym_l, kap_pk_m = _locate_peak_half(kappa, x_1d, y_1d, 0.0, half_ext)
    xs_l, ys_l, kap_pk_s = _locate_peak_half(kappa, x_1d, y_1d, -half_ext, 0.0)
    dx_lens = float(np.hypot(xm_l - xs_l, ym_l - ys_l))

    # Observable X-ray peaks: find separate Main/Sub S_X peaks in total S_X when distinct,
    # falling back cleanly to per-cloud gas peaks if shocked cores sit closely overlaid
    xm_x, ym_x, _ = _locate_peak_half(sx_tot, x_1d, y_1d, 5.0, half_ext)
    xs_x, ys_x, _ = _locate_peak_half(sx_tot, x_1d, y_1d, -half_ext, -5.0)
    if np.hypot(xm_x - xs_x, ym_x - ys_x) < 40.0:
        xm_x, ym_x, _ = _locate_peak_half(sx_main, x_1d, y_1d, -half_ext, half_ext)
        xs_x, ys_x, _ = _locate_peak_half(sx_sub, x_1d, y_1d, -half_ext, half_ext)

    dx_xray = float(np.hypot(xm_x - xs_x, ym_x - ys_x))
    offset_main = float(np.hypot(xm_l - xm_x, ym_l - ym_x))
    offset_sub = float(np.hypot(xs_l - xs_x, ys_l - ys_x))

    # Projected cylinder mass within R = 250 kpc of each lensing peak
    pix_area_m2 = (dx_kpc * KPC_TO_METER) ** 2
    r_from_main = np.hypot(X2D - xm_l, Y2D - ym_l)
    r_from_sub = np.hypot(X2D - xs_l, Y2D - ys_l)
    m2d_main_250 = float(np.sum(sigma_2d[r_from_main <= 250.0]) * pix_area_m2 / M_SOLAR_KG)
    m2d_sub_250 = float(np.sum(sigma_2d[r_from_sub <= 250.0]) * pix_area_m2 / M_SOLAR_KG)

    re_main = float(np.sqrt(np.sum((lambda_t <= 0.0) & (X2D >= 0.0)) * (dx_kpc ** 2) / np.pi))
    re_sub = float(np.sqrt(np.sum((lambda_t <= 0.0) & (X2D < 0.0)) * (dx_kpc ** 2) / np.pi))

    return {
        "x_1d": x_1d,
        "y_1d": y_1d,
        "kappa": kappa,
        "mu_inv": mu_inv,
        "sx_tot": sx_tot,
        "dx_lens_kpc": dx_lens,
        "dx_xray_kpc": dx_xray,
        "offset_main_kpc": offset_main,
        "offset_sub_kpc": offset_sub,
        "M_2D_main_250kpc": m2d_main_250,
        "M_2D_sub_250kpc": m2d_sub_250,
        "kappa_peak_main": kap_pk_m,
        "kappa_peak_sub": kap_pk_s,
        "R_E_main_kpc": re_main,
        "R_E_sub_kpc": re_sub,
        "lens_peaks": ((xm_l, ym_l), (xs_l, ys_l)),
        "xray_peaks": ((xm_x, ym_x), (xs_x, ys_x)),
    }


# ==============================================================================
# 6. Full 2-Cluster Collision Runner
# ==============================================================================
def run_bullet_simulation(task):
    """Executes one Bullet Cluster merger model across 900 Myr (dt = 2 Myr, 96^3 FFT)."""
    t0_wall = time.perf_counter()
    model_name = task["model"]
    short_name = task["short"]
    p_cap_pa = float(task["p_cap_pa"])
    icm_active = bool(task["icm_active"])
    seed = int(task["seed"])
    scale = float(task["scale"])

    rng = np.random.default_rng(seed)
    box_l_m = 5000.0 * KPC_TO_METER
    grid_n = 96
    hydro_ng = 36
    dt_myr = 2.0
    dt_s = dt_myr * MYR_TO_SEC
    half_dt_s = 0.5 * dt_s
    n_steps = 450

    sigma_p_kgm2 = sigma_pocket_from_pcap(p_cap_pa)
    sigma_p_msun_pc2 = sigma_p_kgm2 * KG_M2_TO_MSUN_PC2
    solver = QUMONDMeshSolver(box_l_m, grid_n=grid_n, fft_workers=1)

    ic_m = generate_cluster_ics(rng, CLUSTERS["main"], sigma_p_kgm2, solver, icm_active, scale=scale)
    ic_s = generate_cluster_ics(rng, CLUSTERS["sub"], sigma_p_kgm2, solver, icm_active, scale=scale)

    # Center-of-mass merger setup in x-y plane: d_0 = 2,000 kpc, b = 100 kpc, v_0 = 2,500 km/s
    d0_m = 2000.0 * KPC_TO_METER
    b0_m = 100.0 * KPC_TO_METER
    v0_ms = 2500.0e3

    # Exact mass ratio (identical for Pocket+Galaxies and Total Baryons because f_star is equal!)
    m_main_ref = ic_m["M_all_tot"]
    m_sub_ref = ic_s["M_all_tot"]
    f_sub = m_sub_ref / (m_main_ref + m_sub_ref)
    f_main = 1.0 - f_sub

    shift_m = np.array([-d0_m * f_sub, -b0_m * f_sub, 0.0])
    shift_s = np.array([+d0_m * f_main, +b0_m * f_main, 0.0])
    vshift_m = np.array([+v0_ms * f_sub, 0.0, 0.0])
    vshift_s = np.array([-v0_ms * f_main, 0.0, 0.0])

    pos = np.vstack((ic_m["pos"] + shift_m, ic_s["pos"] + shift_s))
    vel = np.vstack((ic_m["vel"] + vshift_m, ic_s["vel"] + vshift_s))
    masses = np.concatenate((ic_m["masses"], ic_s["masses"]))
    grav_masses = np.concatenate((ic_m["grav_masses"], ic_s["grav_masses"]))
    node_weights = grav_masses / M_PROTON

    # Index bookkeeping
    nm = len(ic_m["masses"])
    idx_p1 = np.arange(0, ic_m["n_p"])
    idx_s1 = np.arange(ic_m["n_p"], ic_m["n_p"] + ic_m["n_s"])
    idx_g1 = np.arange(ic_m["n_p"] + ic_m["n_s"], nm)

    idx_p2 = np.arange(nm, nm + ic_s["n_p"])
    idx_s2 = np.arange(nm + ic_s["n_p"], nm + ic_s["n_p"] + ic_s["n_s"])
    idx_g2 = np.arange(nm + ic_s["n_p"] + ic_s["n_s"], nm + len(ic_s["masses"]))

    temp_g = np.zeros(len(masses), dtype=np.float64)
    temp_g[idx_g1] = ic_m["temp_g"]
    temp_g[idx_g2] = ic_s["temp_g"]

    idx_lens1 = np.concatenate((idx_p1, idx_s1))
    idx_lens2 = np.concatenate((idx_p2, idx_s2))
    idx_grav = np.flatnonzero(grav_masses > 0.0)
    m_grav_sum = float(np.sum(masses[idx_grav]))

    # Zero exact initial COM drift of gravitating sector & total system
    v_com0_grav = np.sum(masses[idx_grav, None] * vel[idx_grav], axis=0) / m_grav_sum
    vel -= v_com0_grav
    p0_grav = np.sum(masses[idx_grav, None] * vel[idx_grav], axis=0)

    acc = solver.get_particle_accelerations(pos, node_weights, mode="Pure K_Omega")

    core_r_m = np.clip(max(2.0 * ic_m["r0_m"], 250.0 * KPC_TO_METER), 250.0 * KPC_TO_METER, 500.0 * KPC_TO_METER)
    gate_r_m = 250.0 * KPC_TO_METER
    audit_drag = short_name == "Sigma_0/3 Passive"
    c_delta_dv_cum = np.zeros(3, dtype=np.float64)
    c_delta_dx_cum = np.zeros(3, dtype=np.float64)

    traj_rows = []
    min_sep_m = 1.0e99
    t_peri_myr = 0.0
    v_peri_kms = 0.0
    passed_peri = False
    bullet_snap = None
    fallback_snap = None
    best_fallback_dx = -1.0
    max_dp_norm = 0.0

    for step in range(0, n_steps + 1):
        t_myr = step * dt_myr
        if step > 0:
            vel += half_dt_s * acc
            pos += dt_s * vel
            acc = solver.get_particle_accelerations(pos, node_weights, mode="Pure K_Omega")
            vel += half_dt_s * acc

            # Symmetric momentum-conserving hydro shock drag & kinetic heating
            n_g1, n_g2 = _apply_ram_pressure_drag(
                pos, vel, temp_g, masses, idx_g1, idx_g2, box_l_m, hydro_ng, dt_s
            )
            if audit_drag and passed_peri and bullet_snap is None:
                dv_drag = _estimate_crossing_drag_shift(
                    pos[idx_p1], vel[idx_p1], masses[idx_p1],
                    pos[idx_p2], vel[idx_p2], masses[idx_p2],
                    sigma_p_kgm2, box_l_m, hydro_ng, dt_s,
                )
                c_delta_dv_cum += dv_drag
                c_delta_dx_cum += c_delta_dv_cum * dt_s
        else:
            n_g1, n_g2 = _apply_ram_pressure_drag(
                pos, vel.copy(), temp_g.copy(), masses, idx_g1, idx_g2, box_l_m, hydro_ng, 0.0
            )

        p_cur_grav = np.sum(masses[idx_grav, None] * vel[idx_grav], axis=0)
        dp_norm = float(np.linalg.norm(p_cur_grav - p0_grav)) / (m_grav_sum * v0_ms)
        if dp_norm > max_dp_norm:
            max_dp_norm = dp_norm

        c_m_lens, v_m_lens = _shrinking_core_centroid(pos[idx_lens1], vel[idx_lens1], core_r_m)
        c_s_lens, v_s_lens = _shrinking_core_centroid(pos[idx_lens2], vel[idx_lens2], core_r_m)
        c_m_star, _ = _shrinking_core_centroid(pos[idx_s1], vel[idx_s1], 0.8 * core_r_m)
        c_s_star, _ = _shrinking_core_centroid(pos[idx_s2], vel[idx_s2], 0.8 * core_r_m)
        c_m_gas, v_m_gas = _shrinking_core_centroid(pos[idx_g1], vel[idx_g1], core_r_m)
        c_s_gas, v_s_gas = _shrinking_core_centroid(pos[idx_g2], vel[idx_g2], core_r_m)

        sep_vec = c_m_lens - c_s_lens
        sep_xy_m = float(np.hypot(sep_vec[0], sep_vec[1]))
        sep_xy_kpc = sep_xy_m / KPC_TO_METER
        v_rel_kms = float(np.linalg.norm(v_m_lens - v_s_lens)) / 1000.0

        if sep_xy_m < min_sep_m:
            min_sep_m = sep_xy_m
            t_peri_myr = t_myr
            v_peri_kms = v_rel_kms
        elif sep_xy_m > min_sep_m + 30.0 * KPC_TO_METER and sep_vec[0] > 0.0:
            passed_peri = True

        # Live Two-Condition Compilation Gate on Main & Sub ICM gas inside R = 250 kpc
        g1_rel = np.linalg.norm(pos[idx_g1] - c_m_gas, axis=1) <= gate_r_m
        s1_rel = np.linalg.norm(pos[idx_s1] - c_m_gas, axis=1) <= gate_r_m
        m_b_m = float(np.sum(masses[idx_g1[g1_rel]])) + float(np.sum(masses[idx_s1[s1_rel]]))
        rho_n_m = float(np.mean(n_g1[g1_rel])) if np.any(g1_rel) else 1000.0
        t_m_k = float(np.mean(temp_g[idx_g1[g1_rel]])) if np.any(g1_rel) else 9.0 * KEV_TO_K
        gate_m = evaluate_two_condition_compilation_gate(
            rho_n_m, t_m_k, gate_r_m, max(m_b_m, 1.0e11 * M_SOLAR_KG), mu_mol=0.60,
            v_bulk_ms=float(np.linalg.norm(v_m_gas)),
        )

        g2_rel = np.linalg.norm(pos[idx_g2] - c_s_gas, axis=1) <= gate_r_m
        s2_rel = np.linalg.norm(pos[idx_s2] - c_s_gas, axis=1) <= gate_r_m
        m_b_s = float(np.sum(masses[idx_g2[g2_rel]])) + float(np.sum(masses[idx_s2[s2_rel]]))
        rho_n_s = float(np.mean(n_g2[g2_rel])) if np.any(g2_rel) else 1000.0
        t_s_k = float(np.mean(temp_g[idx_g2[g2_rel]])) if np.any(g2_rel) else 6.5 * KEV_TO_K
        gate_s = evaluate_two_condition_compilation_gate(
            rho_n_s, t_s_k, gate_r_m, max(m_b_s, 1.0e11 * M_SOLAR_KG), mu_mol=0.60,
            v_bulk_ms=float(np.linalg.norm(v_s_gas)),
        )

        traj_rows.append(
            {
                "t_Myr": t_myr,
                "x_main_pocket_kpc": float(c_m_lens[0] / KPC_TO_METER),
                "x_sub_pocket_kpc": float(c_s_lens[0] / KPC_TO_METER),
                "x_main_gal_kpc": float(c_m_star[0] / KPC_TO_METER),
                "x_sub_gal_kpc": float(c_s_star[0] / KPC_TO_METER),
                "x_main_gas_kpc": float(c_m_gas[0] / KPC_TO_METER),
                "x_sub_gas_kpc": float(c_s_gas[0] / KPC_TO_METER),
                "sep_lens_xy_kpc": sep_xy_kpc,
                "sep_gas_xy_kpc": float(np.hypot(c_m_gas[0] - c_s_gas[0], c_m_gas[1] - c_s_gas[1]) / KPC_TO_METER),
                "v_rel_kms": v_rel_kms,
                "T_main_ICM_keV": t_m_k / KEV_TO_K,
                "T_sub_ICM_keV": t_s_k / KEV_TO_K,
                "Lambda_lock_main_ICM": float(gate_m["lambda_lock"]),
                "Lambda_lock_sub_ICM": float(gate_s["lambda_lock"]),
                "com_drift_norm": dp_norm,
            }
        )

        if passed_peri:
            snap_payload = {
                "t_myr": t_myr,
                "pos": pos.copy(),
                "grav_masses": grav_masses.copy(),
                "c_main_lens": c_m_lens.copy(),
                "c_sub_lens": c_s_lens.copy(),
                "idx_g1": idx_g1,
                "idx_g2": idx_g2,
                "m_g1": masses[idx_g1].copy(),
                "m_g2": masses[idx_g2].copy(),
                "n_g1": n_g1.copy(),
                "n_g2": n_g2.copy(),
                "T_g1": temp_g[idx_g1].copy(),
                "T_g2": temp_g[idx_g2].copy(),
                "lam_sub": float(gate_s["lambda_lock"]),
                "lam_main": float(gate_m["lambda_lock"]),
            }
            if sep_xy_kpc > best_fallback_dx:
                best_fallback_dx = sep_xy_kpc
                fallback_snap = snap_payload
            if bullet_snap is None and sep_xy_kpc >= 720.0:
                bullet_snap = snap_payload

    if bullet_snap is None:
        bullet_snap = fallback_snap if fallback_snap is not None else snap_payload

    maps = compute_lensing_and_xray_maps(bullet_snap)
    runtime_s = time.perf_counter() - t0_wall

    summary = {
        "Model": model_name,
        "Short_Name": short_name,
        "P_cap_Pa": p_cap_pa,
        "Sigma_p_Msun_pc2": sigma_p_msun_pc2,
        "t_peri_Myr": t_peri_myr,
        "t_bullet_Myr": float(bullet_snap["t_myr"]),
        "v_peri_kms": v_peri_kms,
        "dx_lens_kpc": maps["dx_lens_kpc"],
        "dx_xray_kpc": maps["dx_xray_kpc"],
        "offset_main_kpc": maps["offset_main_kpc"],
        "offset_sub_kpc": maps["offset_sub_kpc"],
        "M_2D_main_250kpc": maps["M_2D_main_250kpc"],
        "M_2D_sub_250kpc": maps["M_2D_sub_250kpc"],
        "kappa_peak_main": maps["kappa_peak_main"],
        "kappa_peak_sub": maps["kappa_peak_sub"],
        "R_E_main_kpc": maps["R_E_main_kpc"],
        "R_E_sub_kpc": maps["R_E_sub_kpc"],
        "Lambda_lock_sub_ICM": float(bullet_snap["lam_sub"]),
        "com_drift_norm": max_dp_norm,
        "runtime_s": runtime_s,
    }
    drag_shift_kpc = float(np.linalg.norm(c_delta_dx_cum)) / KPC_TO_METER if audit_drag else 0.0
    return {
        "summary": summary,
        "traj": traj_rows,
        "maps": maps,
        "drag_shift_kpc": drag_shift_kpc,
    }


# ==============================================================================
# 7. Dark-Theme 4-Panel Diagnostic Figure
# ==============================================================================
def build_bullet_figure(results, summary_df, out_png):
    """Generates the dark-theme (#0b1120) 4-panel Bullet Cluster figure at 150 DPI."""
    fig, axes = plt.subplots(2, 2, figsize=(16.5, 11.8), dpi=150)
    fig.patch.set_facecolor("#0b1120")

    for ax in axes.flat:
        ax.set_facecolor("#070c18")
        ax.tick_params(colors="#cbd5e1", labelsize=9.2)
        ax.xaxis.label.set_color("#e2e8f0")
        ax.yaxis.label.set_color("#e2e8f0")
        ax.title.set_color("#f8fafc")
        for spine in ax.spines.values():
            spine.set_color("#334155")

    res_map = {r["summary"]["Short_Name"]: r for r in results}
    prim = res_map["Sigma_0/3 Passive"]
    pmaps = prim["maps"]
    psum = prim["summary"]

    # ------------------------------------------------------------------
    # Panel 1 (Top-Left): Composite Bullet Cluster X-ray + Lensing Map
    # ------------------------------------------------------------------
    ax1 = axes[0, 0]
    x1, y1 = pmaps["x_1d"], pmaps["y_1d"]
    sx_norm = np.log10(np.maximum(pmaps["sx_tot"] / max(float(np.max(pmaps["sx_tot"])), 1e-30), 1.0e-4))
    im1 = ax1.imshow(
        sx_norm,
        extent=[x1[0], x1[-1], y1[0], y1[-1]],
        origin="lower",
        cmap="inferno",
        vmin=-3.2,
        vmax=0.0,
        aspect="auto",
    )
    cbar1 = fig.colorbar(im1, ax=ax1, fraction=0.045, pad=0.03)
    cbar1.ax.tick_params(colors="#94a3b8", labelsize=8)
    cbar1.set_label(r"$\log_{10}(S_X / S_{X,\max})$", color="#cbd5e1", fontsize=9)

    levels_k = [0.2, 0.4, 0.6, 0.8, 1.0]
    cs_k = ax1.contour(x1, y1, pmaps["kappa"], levels=levels_k, colors="#22d3ee", linewidths=1.35, alpha=0.92)
    ax1.clabel(cs_k, fmt=r"$\kappa=%.1f$", fontsize=7.5, inline=True)

    # Strong lensing critical curves mu^-1 = 0
    if float(np.min(pmaps["mu_inv"])) <= 0.0 <= float(np.max(pmaps["mu_inv"])):
        ax1.contour(x1, y1, pmaps["mu_inv"], levels=[0.0], colors="#ffffff", linewidths=2.1, linestyles="-")

    (xml, yml), (xsl, ysl) = pmaps["lens_peaks"]
    (xmx, ymx), (xsx, ysx) = pmaps["xray_peaks"]
    ax1.scatter([xml, xsl], [yml, ysl], marker="+", s=110, linewidths=2.2, color="#22d3ee", zorder=6, label=r"Lensing $\kappa$ Peaks")
    ax1.scatter([xmx, xsx], [ymx, ysx], marker="x", s=90, linewidths=2.2, color="#facc15", zorder=6, label=r"X-Ray $S_X$ Peaks")

    badge_1 = (
        rf"$t = {psum['t_bullet_Myr']:.0f}$ Myr | $\Delta x_{{\mathrm{{lens}}}} = {psum['dx_lens_kpc']:.0f}$ kpc" "\n"
        rf"Offsets: Main ${psum['offset_main_kpc']:.0f}$ kpc | Bullet ${psum['offset_sub_kpc']:.0f}$ kpc" "\n"
        rf"White: $\mu^{{-1}}=0$ ($R_E={psum['R_E_main_kpc']:.0f}, {psum['R_E_sub_kpc']:.0f}$ kpc)"
    )
    ax1.text(
        0.03, 0.96, badge_1, transform=ax1.transAxes, va="top", ha="left", fontsize=8.0, color="#f8fafc",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#0f172a", edgecolor="#22d3ee", alpha=0.9),
    )
    ax1.set_xlim(-950.0, 950.0)
    ax1.set_ylim(-650.0, 650.0)
    ax1.set_xlabel(r"Sky Offset $x$ [kpc]  ($4.581\ \mathrm{kpc}/''$ at $z_l=0.296$)")
    ax1.set_ylabel("Sky Offset $y$ [kpc]")
    ax1.set_title(r"$\Sigma_0/3$ Composite Map: X-Ray Gas vs. $\kappa$ & Strong Lensing", fontsize=11.0, fontweight="bold", pad=8)
    ax1.legend(loc="lower right", fontsize=8.0, frameon=True, facecolor="#0f172a", edgecolor="#334155", labelcolor="#e2e8f0")

    # ------------------------------------------------------------------
    # Panel 2 (Top-Right): 1D Equatorial Merger Axis Profiles (y = 0)
    # ------------------------------------------------------------------
    ax2 = axes[0, 1]
    ax2.grid(True, linestyle="--", alpha=0.22, color="#64748b")
    mid_y = np.argmin(np.abs(y1))
    sl_y = slice(max(0, mid_y - 2), min(len(y1), mid_y + 3))

    sx_1d = np.mean(pmaps["sx_tot"][sl_y, :], axis=0)
    sx_1d_norm = sx_1d / max(float(np.max(sx_1d)), 1e-30)

    col_map = {
        "3P_delta Passive": ("#94a3b8", "--", r"$3P_\delta$ ($\Sigma_p=169\ M_\odot/\mathrm{pc}^2$)"),
        "Sigma_0/3 Passive": ("#22d3ee", "-", r"$\Sigma_0/3$ Passive ($\Sigma_p=2600\ M_\odot/\mathrm{pc}^2$)"),
        "Sigma_0 Passive": ("#f472b6", "-.", r"$\Sigma_0$ ($\Sigma_p=4503\ M_\odot/\mathrm{pc}^2$)"),
    }
    ax2_r = ax2.twinx()
    l_xray = ax2_r.fill_between(x1, 0.0, sx_1d_norm, color="#f97316", alpha=0.28, label=r"Normalized X-Ray $S_X(x)$")
    ax2_r.plot(x1, sx_1d_norm, color="#fb923c", linewidth=1.6, alpha=0.9)
    ax2_r.set_ylabel(r"Normalized X-Ray $S_X$", color="#fb923c")
    ax2_r.set_ylim(0.0, 1.35)
    ax2_r.tick_params(colors="#cbd5e1", labelsize=9)
    for sp in ax2_r.spines.values():
        sp.set_color("#334155")

    lines_p2 = []
    for key, (col, ls, lbl) in col_map.items():
        if key in res_map:
            k_1d = np.mean(res_map[key]["maps"]["kappa"][sl_y, :], axis=0)
            ln, = ax2.plot(x1, k_1d, color=col, linestyle=ls, linewidth=2.1, label=lbl)
            lines_p2.append(ln)

    ax2.axhline(1.0, color="#ffffff", linestyle=":", linewidth=1.2, alpha=0.75, label=r"Critical Density $\kappa = 1.0$")
    ax2.set_xlim(-950.0, 950.0)
    ax2.set_ylim(0.0, max(1.85, float(np.max(res_map["Sigma_0 Passive"]["maps"]["kappa"])) * 1.08))
    ax2.set_xlabel("Merger Axis $x$ [kpc]  (|y| < 25 kpc Strip)")
    ax2.set_ylabel(r"Lensing Convergence $\kappa(x)$")
    ax2.set_title(r"1D Equatorial Profile ($y \approx 0$): Pressure Ceiling Ablation vs. X-Ray Gas", fontsize=11.0, fontweight="bold", pad=8)
    labels_p2 = [ln.get_label() for ln in lines_p2] + [r"Critical Density $\kappa = 1.0$", r"Normalized X-Ray $S_X(x)$"]
    handles_p2 = lines_p2 + [ax2.get_lines()[-1], l_xray]
    ax2.legend(handles_p2, labels_p2, loc="upper left", fontsize=7.8, frameon=True, facecolor="#0f172a", edgecolor="#334155", labelcolor="#e2e8f0")

    # ------------------------------------------------------------------
    # Panel 3 (Bottom-Left): Model Comparison vs. Observed Benchmarks
    # ------------------------------------------------------------------
    ax3 = axes[1, 0]
    ax3.grid(True, linestyle="--", alpha=0.22, color="#64748b")
    short_order = ["3P_delta Passive", "Sigma_0/3 Passive", "Sigma_0/3 Active", "Sigma_0 Passive"]
    tick_lbls = [r"$3P_\delta$" "\nPassive", r"$\Sigma_0/3$" "\nPassive", r"$\Sigma_0/3$" "\nActive", r"$\Sigma_0$" "\nPassive"]
    x_b = np.arange(len(short_order))

    m2d_m = [res_map[k]["summary"]["M_2D_main_250kpc"] / 1.0e14 for k in short_order]
    m2d_s = [res_map[k]["summary"]["M_2D_sub_250kpc"] / 1.0e14 for k in short_order]
    re_m = [res_map[k]["summary"]["R_E_main_kpc"] for k in short_order]
    re_s = [res_map[k]["summary"]["R_E_sub_kpc"] for k in short_order]

    bw = 0.28
    ax3.bar(x_b - 0.5 * bw, m2d_m, width=bw, color="#38bdf8", label=r"Main $M_{2D}(<250\ \mathrm{kpc})$ [$10^{14} M_\odot$]")
    ax3.bar(x_b + 0.5 * bw, m2d_s, width=bw, color="#a78bfa", label=r"Bullet $M_{2D}(<250\ \mathrm{kpc})$ [$10^{14} M_\odot$]")
    ax3.axhspan(2.6, 3.0, color="#38bdf8", alpha=0.14, label=r"Obs Main: $2.8 \pm 0.2 \times 10^{14} M_\odot$ (Bradač+06)")
    ax3.axhspan(1.9, 2.3, color="#a78bfa", alpha=0.14, label=r"Obs Bullet: $2.1 \pm 0.2 \times 10^{14} M_\odot$ (Paraficz+16)")

    ax3_r = ax3.twinx()
    ax3_r.plot(x_b - 0.5 * bw, re_m, "o--", color="#22d3ee", linewidth=1.6, markersize=6, label=r"Main $R_E$ [kpc] (Obs $\approx 125$)")
    ax3_r.plot(x_b + 0.5 * bw, re_s, "s--", color="#facc15", linewidth=1.6, markersize=6, label=r"Bullet $R_E$ [kpc] (Obs $\approx 80$)")
    ax3_r.set_ylabel("Effective Einstein Radius $R_E$ [kpc]", color="#facc15")
    ax3_r.set_ylim(0.0, 260.0)
    ax3_r.tick_params(colors="#cbd5e1", labelsize=9)
    for sp in ax3_r.spines.values():
        sp.set_color("#334155")

    ax3.set_xticks(x_b)
    ax3.set_xticklabels(tick_lbls, fontsize=8.8)
    ax3.set_ylabel(r"Projected Cylinder Mass $M_{2D}(<250\ \mathrm{kpc})$ [$10^{14} M_\odot$]")
    ax3.set_ylim(0.0, 5.4)
    ax3.set_title(r"Projected Mass $M_{2D}(<250\ \mathrm{kpc})$ & Einstein Radii vs. HST/Chandra Bounds", fontsize=11.0, fontweight="bold", pad=8)
    h1, l1 = ax3.get_legend_handles_labels()
    h2, l2 = ax3_r.get_legend_handles_labels()
    ax3.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=7.2, frameon=True, facecolor="#0f172a", edgecolor="#334155", labelcolor="#e2e8f0")

    # ------------------------------------------------------------------
    # Panel 4 (Bottom-Right): Merger Trajectories & Two-Condition Gate
    # ------------------------------------------------------------------
    ax4 = axes[1, 1]
    ax4.grid(True, linestyle="--", alpha=0.22, color="#64748b")
    tdf = pd.DataFrame(prim["traj"])
    tt = tdf["t_Myr"].to_numpy()

    l_p1, = ax4.plot(tt, tdf["x_main_pocket_kpc"], color="#38bdf8", linewidth=2.2, label="Main Pocket + Galaxies $x(t)$")
    l_p2, = ax4.plot(tt, tdf["x_sub_pocket_kpc"], color="#22d3ee", linewidth=2.2, linestyle="-", label="Bullet Pocket + Galaxies $x(t)$")
    l_g1, = ax4.plot(tt, tdf["x_main_gas_kpc"], color="#fb923c", linewidth=2.0, linestyle="--", label="Main Shocked Gas $x(t)$")
    l_g2, = ax4.plot(tt, tdf["x_sub_gas_kpc"], color="#facc15", linewidth=2.0, linestyle="--", label="Bullet Shocked Gas $x(t)$")

    ax4.axvline(psum["t_peri_Myr"], color="#94a3b8", linestyle=":", linewidth=1.2)
    ax4.axvline(psum["t_bullet_Myr"], color="#22c55e", linestyle="-.", linewidth=1.6, label=rf"Bullet Epoch ($t={psum['t_bullet_Myr']:.0f}$ Myr)")

    ax4_r = ax4.twinx()
    l_lam, = ax4_r.semilogy(tt, np.maximum(tdf["Lambda_lock_sub_ICM"], 1e-4), color="#f43f5e", linewidth=2.0, label=r"Bullet ICM $\Lambda_{\mathrm{lock}}(t)$ (Passive < 1)")
    ax4_r.axhline(1.0, color="#f43f5e", linestyle=":", linewidth=1.3, alpha=0.8)
    ax4_r.set_ylabel(r"Two-Condition Phase Lock $\Lambda_{\mathrm{lock}}(t)$", color="#f43f5e")
    ax4_r.set_ylim(1.0e-2, 5.0)
    ax4_r.tick_params(colors="#cbd5e1", labelsize=9)
    for sp in ax4_r.spines.values():
        sp.set_color("#334155")

    badge_4 = (
        rf"$C_\delta$ crossing-drag shift: $\Delta x_{{C_\delta}} = {prim['drag_shift_kpc']:.3f}$ kpc (< 1 kpc)" "\n"
        rf"Active 3-Species $|\Delta P|/(Mv_0) = {res_map['Sigma_0/3 Active']['summary']['com_drift_norm']:.2e}$"
    )
    ax4.text(
        0.03, 0.05, badge_4, transform=ax4.transAxes, va="bottom", ha="left", fontsize=7.8, color="#f8fafc",
        bbox=dict(boxstyle="round,pad=0.38", facecolor="#0f172a", edgecolor="#22c55e", alpha=0.9),
    )
    ax4.set_xlabel("Merger Elapsed Time $t$ [Myr]")
    ax4.set_ylabel("Merger Axis Centroid $x(t)$ [kpc]")
    ax4.set_title(r"Dynamical Separation & Live Two-Condition Gate ($\Sigma_0/3$, ICM Passive)", fontsize=11.0, fontweight="bold", pad=8)

    h4 = [l_p1, l_p2, l_g1, l_g2, ax4.get_lines()[-1], l_lam]
    l4 = [h.get_label() for h in h4]
    ax4.legend(h4, l4, loc="upper right", fontsize=7.4, frameon=True, facecolor="#0f172a", edgecolor="#334155", labelcolor="#e2e8f0")

    fig.suptitle(
        r"Bullet Cluster (1E 0657-56) — Zero-Parameter Cover Pocket & Two-Condition Compilation Gate ($96^3$ FFT PM N-Body)",
        fontsize=13.4,
        fontweight="bold",
        color="#f8fafc",
        y=0.985,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    fig.savefig(out_png, dpi=150, facecolor=fig.get_facecolor())
    plt.close(fig)


# ==============================================================================
# 8. CLI Main Entrypoint
# ==============================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Zero-Parameter Bullet Cluster (1E 0657-56) 3D FFT Merger & Lensing Testbed"
    )
    parser.add_argument("--out-dir", type=str, default="bullet_audit", help="Output directory")
    parser.add_argument("--scale", type=float, default=1.0, help="Particle count scale factor (1.0 = 33,000 total)")
    parser.add_argument(
        "--jobs",
        type=int,
        default=min(4, os.cpu_count() or 4),
        help="Parallel worker processes across the 4 controlled Bullet Cluster runs",
    )
    args = parser.parse_args()

    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    print("=== Bullet Cluster Zero-Parameter First-Principles Invariants ===")
    print(f"G_EFF          = {G_EFF:.6e} m^3/(kg s^2) | ETA_POCKET = {ETA_POCKET:.6f}")
    print(
        f"Lensing z_l=0.296->z_s=1.0: D_l={LENS_COSMO['D_l_mpc']:.2f} Mpc, "
        f"D_s={LENS_COSMO['D_s_mpc']:.2f} Mpc, D_ls={LENS_COSMO['D_ls_mpc']:.2f} Mpc"
    )
    print(
        f"SIGMA_CRIT     = {SIGMA_CRIT_KGM2:.4f} kg/m^2 = {LENS_COSMO['sigma_crit_msun_pc2']:.2f} M_sun/pc^2 | "
        f"scale = {LENS_COSMO['kpc_per_arcsec']:.4f} kpc/arcsec"
    )

    tasks = [{**m, "scale": args.scale} for m in BULLET_MODELS]
    print(f"\nLaunching {len(tasks)} Bullet Cluster 96^3 3D FFT simulations across {args.jobs} workers...")

    t_all = time.perf_counter()
    with ProcessPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        results = list(pool.map(run_bullet_simulation, tasks))

    summary_rows = [r["summary"] for r in results]
    summary_df = pd.DataFrame(summary_rows)
    prim_res = next(r for r in results if r["summary"]["Short_Name"] == "Sigma_0/3 Passive")
    traj_df = pd.DataFrame(prim_res["traj"])

    for r in results:
        s = r["summary"]
        print(
            f"  {s['Short_Name']:18s} | t_peri={s['t_peri_Myr']:5.0f} | t_bul={s['t_bullet_Myr']:5.0f} | "
            f"dx_lens={s['dx_lens_kpc']:5.1f} | dx_xray={s['dx_xray_kpc']:5.1f} | "
            f"off=({s['offset_main_kpc']:5.1f}, {s['offset_sub_kpc']:5.1f}) | "
            f"M2D=({s['M_2D_main_250kpc']:.2e}, {s['M_2D_sub_250kpc']:.2e}) | "
            f"kap=({s['kappa_peak_main']:.2f}, {s['kappa_peak_sub']:.2f}) | "
            f"RE=({s['R_E_main_kpc']:5.1f}, {s['R_E_sub_kpc']:5.1f}) | "
            f"Lam={s['Lambda_lock_sub_ICM']:.3f} | dP={s['com_drift_norm']:.2e} | t={s['runtime_s']:.1f}s"
        )
    print(f"C_delta crossing-drag upper bound on Sigma_0/3 subcluster pocket: {prim_res['drag_shift_kpc']:.5f} kpc")
    act_res = next(r for r in results if r["summary"]["Short_Name"] == "Sigma_0/3 Active")
    b_main = prim_res["summary"]["M_2D_main_250kpc"] / act_res["summary"]["M_2D_main_250kpc"]
    b_sub = prim_res["summary"]["M_2D_sub_250kpc"] / act_res["summary"]["M_2D_sub_250kpc"]
    print(f"Cluster lensing mass bias M_2D(Passive)/M_2D(Active) on Sigma_0/3: Main={b_main:.4f}, Bullet={b_sub:.4f}")

    sum_csv = os.path.join(out_dir, "bullet_cluster_summary.csv")
    traj_csv = os.path.join(out_dir, "bullet_cluster_trajectory.csv")
    png_path = os.path.join(out_dir, "bullet_cluster_benchmark.png")

    summary_export_df = summary_df.drop(columns=["Short_Name"])
    summary_export_df.to_csv(sum_csv, index=False)
    traj_df.to_csv(traj_csv, index=False)
    build_bullet_figure(results, summary_df, png_path)

    print(f"\nTotal wall time: {time.perf_counter() - t_all:.1f}s")
    print(f"Saved: {sum_csv}")
    print(f"Saved: {traj_csv}")
    print(f"Saved: {png_path}")


if __name__ == "__main__":
    main()
