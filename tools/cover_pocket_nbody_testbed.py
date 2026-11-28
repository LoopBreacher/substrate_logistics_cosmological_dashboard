"""
Standalone Zero-Parameter Topological Cover Pocket N-Body Testbed.

Tests the Topological Cover Pocket architecture (strictly linear Pure K_Omega
gravity + cored Burkert Cover Pockets anchored by the 3D Topological Drag
Pressure Floor 3 P_delta) across three built-in verification experiments:
  - Experiment A: 6-Satellite live-star Plummer benchmark (Palomar 14, NGC 2419,
    Segue 1, Draco dSph, Ursa Minor dSph, Sculptor dSph) with exact isotropic
    spherical Jeans initialization in Pure K_Omega and analytic 3 P_delta pockets.
  - Experiment B: Live two-component self-gravitating Draco dSph proof
    (10,000 live stars + 20,000 live collisionless Burkert pocket nodes) evolved
    for 50 crossing times under strictly linear Pure K_Omega 3D FFT PM gravity.
  - Experiment C: Zero-parameter mass-locked Burkert + Freeman thin-disk
    rotation curve cross-check on DDO 154, NGC 3198, Milky Way, and UGC 2885.

Imports all physical constants directly from `modules.data_loader` (zero CODATA G),
reuses `QUMONDMeshSolver` from `tools.qumond_nbody_testbed`, and compares against
`calculate_rotation_profile` from `modules.substrate_rotation_curves` without
modifying either file.
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
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize_scalar

from modules.data_loader import (
    A_OMEGA,
    K_OMEGA,
    KPC_TO_METER,
    M_PROTON,
    M_SOLAR_KG,
    OMEGA_M_GLOBAL,
    OMEGA_M_VOID,
    P_DELTA_MACRO_PA,
)
from modules.substrate_rotation_curves import calculate_rotation_profile
from tools.qumond_nbody_testbed import QUMONDMeshSolver

# ==============================================================================
# 1. Exact First-Principles Cover Pocket Algebraic Invariants (Zero CODATA G)
# ==============================================================================
PC_TO_METER = KPC_TO_METER / 1000.0
G_EFF = K_OMEGA / M_PROTON
DELTA_OMEGA_M = OMEGA_M_GLOBAL - OMEGA_M_VOID
ETA_POCKET = OMEGA_M_VOID / DELTA_OMEGA_M
G_3PDELTA = float(np.sqrt(8.0 * np.pi * G_EFF * 3.0 * P_DELTA_MACRO_PA))


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

SIGMA_POCKET_KGM2 = G_3PDELTA / (F_B_MAX * np.pi * G_EFF)
SIGMA_POCKET_MSUN_PC2 = SIGMA_POCKET_KGM2 * (PC_TO_METER ** 2) / M_SOLAR_KG
A_BTFR = np.pi * G_EFF * SIGMA_POCKET_KGM2 * ETA_POCKET * ((V2_B_MAX ** 2) / M_FAC_VMAX)

# ==============================================================================
# 2. Benchmark Catalogs
# ==============================================================================
SATELLITE_CATALOG = [
    {
        "target": "Palomar 14",
        "kind": "GC",
        "M_solar": 1.0e4,
        "r_h_pc": 25.0,
        "obs_sigma_kms": 0.38,
        "obs_err_kms": 0.12,
        "ref": "Jordi+2009",
        "r0_modes": [("Pure K_Omega", 0.0)],
    },
    {
        "target": "NGC 2419",
        "kind": "GC",
        "M_solar": 9.0e5,
        "r_h_pc": 20.0,
        "obs_sigma_kms": 4.14,
        "obs_err_kms": 0.48,
        "ref": "Baumgardt+2009",
        "r0_modes": [("Pure K_Omega", 0.0)],
    },
    {
        "target": "Segue 1",
        "kind": "UFD",
        "M_solar": 3.4e2,
        "r_h_pc": 29.0,
        "obs_sigma_kms": 3.70,
        "obs_err_kms": 1.40,
        "ref": "Simon+2011",
        "r0_modes": [
            ("Pure K_Omega", 0.0),
            ("Cover Pocket (r0=r_h)", 29.0),
            ("Cover Pocket (r0=300pc)", 300.0),
        ],
    },
    {
        "target": "Draco dSph",
        "kind": "dSph",
        "M_solar": 3.0e5,
        "r_h_pc": 220.0,
        "obs_sigma_kms": 9.10,
        "obs_err_kms": 1.20,
        "ref": "Walker+2007",
        "r0_modes": [
            ("Pure K_Omega", 0.0),
            ("Cover Pocket (r0=r_h)", 220.0),
            ("Cover Pocket (r0=300pc)", 300.0),
        ],
    },
    {
        "target": "Ursa Minor dSph",
        "kind": "dSph",
        "M_solar": 2.9e5,
        "r_h_pc": 181.0,
        "obs_sigma_kms": 9.50,
        "obs_err_kms": 1.20,
        "ref": "Walker+2009",
        "r0_modes": [
            ("Pure K_Omega", 0.0),
            ("Cover Pocket (r0=r_h)", 181.0),
            ("Cover Pocket (r0=300pc)", 300.0),
        ],
    },
    {
        "target": "Sculptor dSph",
        "kind": "dSph",
        "M_solar": 2.3e6,
        "r_h_pc": 280.0,
        "obs_sigma_kms": 9.20,
        "obs_err_kms": 1.10,
        "ref": "Walker+2009",
        "r0_modes": [
            ("Pure K_Omega", 0.0),
            ("Cover Pocket (r0=r_h)", 280.0),
            ("Cover Pocket (r0=300pc)", 300.0),
        ],
    },
]

DISK_CATALOG = [
    {
        "name": "DDO 154",
        "M_B_solar": 3.8e8,
        "R_d_kpc": 1.20,
        "r_focus_kpc": 5.0,
        "r_max_kpc": 10.0,
        "v_obs_kms": 47.0,
    },
    {
        "name": "NGC 3198",
        "M_B_solar": 3.2e10,
        "R_d_kpc": 2.72,
        "r_focus_kpc": 20.0,
        "r_max_kpc": 35.0,
        "v_obs_kms": 150.0,
    },
    {
        "name": "Milky Way",
        "M_B_solar": 8.8e10,
        "R_d_kpc": 2.30,
        "r_focus_kpc": 8.12,
        "r_max_kpc": 30.0,
        "v_obs_kms": 229.0,
    },
    {
        "name": "UGC 2885",
        "M_B_solar": 3.5e11,
        "R_d_kpc": 11.50,
        "r_focus_kpc": 50.0,
        "r_max_kpc": 70.0,
        "v_obs_kms": 298.0,
    },
]


# ==============================================================================
# 3. Spherical Helpers, Burkert Field & Exact Isotropic Jeans Samplers
# ==============================================================================
def burkert_acceleration_vec(pos_m, r0_m):
    """Analytic radial inward acceleration from a cored 3 P_delta Burkert Cover Pocket."""
    if r0_m <= 0.0:
        return np.zeros_like(pos_m)
    r = np.linalg.norm(pos_m, axis=1)
    r_safe = np.where(r > 0.0, r, 1.0)
    g_mag = np.pi * G_EFF * SIGMA_POCKET_KGM2 * burkert_f(r_safe / r0_m)
    scale = np.where(r > 0.0, -g_mag / r_safe, 0.0)
    return pos_m * scale[:, None]


def _sample_isotropic_positions(rng, radii_m):
    """Draws isotropic 3D Cartesian coordinates for given radial distances."""
    n = len(radii_m)
    cos_t = rng.uniform(-1.0, 1.0, size=n)
    phi = rng.uniform(0.0, 2.0 * np.pi, size=n)
    sin_t = np.sqrt(1.0 - cos_t * cos_t)
    return np.column_stack(
        (radii_m * sin_t * np.cos(phi), radii_m * sin_t * np.sin(phi), radii_m * cos_t)
    )


def _sample_plummer_radii(rng, n_particles, a_p_m, r_trunc_ratio):
    """Samples exact truncated Plummer radii via inverse CDF inside r_cut = r_trunc_ratio * a_p."""
    u_max = ((r_trunc_ratio ** 2) / (1.0 + r_trunc_ratio ** 2)) ** 1.5
    u = rng.uniform(1.0e-12, u_max, size=n_particles)
    return a_p_m / np.sqrt(u ** (-2.0 / 3.0) - 1.0)


def _sample_burkert_radii(rng, n_particles, r0_m, r_trunc_ratio=4.0):
    """Samples exact truncated Burkert radii via inverse CDF of m_B(x)/m_B(r_trunc_ratio)."""
    x_tab = np.linspace(0.0, r_trunc_ratio, 4000)
    cdf = burkert_m(x_tab) / burkert_m(r_trunc_ratio)
    cdf[0] = 0.0
    cdf[-1] = 1.0
    u = rng.uniform(0.0, 1.0, size=n_particles)
    return np.interp(u, cdf, x_tab) * r0_m


def _jeans_sigma_r(r_grid_m, rho_grid, g_tot_grid):
    """
    Computes exact isotropic spherical Jeans velocity dispersion sigma_r(r) by integrating
    strictly from r to r_cut (= r_grid_m[-1]) so boundary pressure vanishes cleanly.
    """
    integrand = rho_grid * g_tot_grid
    i_fwd = cumulative_trapezoid(integrand, r_grid_m, initial=0.0)
    tail = np.maximum(i_fwd[-1] - i_fwd, 0.0)
    return np.sqrt(tail / np.maximum(rho_grid, 1.0e-60))


def _sample_jeans_velocities(rng, radii_m, r_grid_m, sigma_r_grid, v_esc_grid):
    """
    Draws isotropic 3D Gaussian velocities with local std sigma_r(r_i) in each axis,
    rejecting unbound speeds v > 0.98 v_esc(r_i).
    """
    n = len(radii_m)
    sig_p = np.interp(radii_m, r_grid_m, sigma_r_grid)
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


def _compute_metrics(pos, vel, acc, masses, origin_pos=None):
    """Evaluates mean 1D LOS velocity dispersion, half-mass radius, and Clausius virial 2K/|W|."""
    total_m = float(np.sum(masses))
    v_com = np.sum(masses[:, None] * vel, axis=0) / total_m
    v_rel = vel - v_com
    # Mass-weighted 1D LOS velocity dispersion averaged across 3 orthogonal projections
    sig_axes = np.sqrt(np.sum(masses[:, None] * (v_rel * v_rel), axis=0) / total_m)
    sigma_los = float(np.mean(sig_axes))

    if origin_pos is None:
        origin_pos = np.sum(masses[:, None] * pos, axis=0) / total_m
    p_rel = pos - origin_pos
    r_vals = np.linalg.norm(p_rel, axis=1)
    order = np.argsort(r_vals)
    cum_m = np.cumsum(masses[order])
    idx_half = int(np.searchsorted(cum_m, 0.5 * total_m))
    idx_half = min(max(idx_half, 0), len(r_vals) - 1)
    r_half = float(r_vals[order[idx_half]])

    k2 = float(np.sum(masses * np.sum(v_rel * v_rel, axis=1)))
    w_clausius = float(np.sum(masses * np.sum(p_rel * acc, axis=1)))
    virial = k2 / abs(w_clausius) if abs(w_clausius) > 0.0 else 1.0
    return sigma_los, r_half, virial, v_com


def generate_jeans_plummer_ics(
    n_stars, a_p_m, total_star_mass_kg, solver, r0_m=0.0, r_trunc_ratio=4.5, seed=42
):
    """
    Initializes stars from a truncated Plummer density profile using the exact
    2,000-point spherical Jeans integral in g_star(r) + g_Burkert(r), then rescales
    velocities so the initial Clausius virial ratio on the 64^3 PM mesh is 1.0000.
    """
    rng = np.random.default_rng(seed)
    radii = _sample_plummer_radii(rng, n_stars, a_p_m, r_trunc_ratio)
    pos = _sample_isotropic_positions(rng, radii)
    pos -= np.mean(pos, axis=0)
    radii_centered = np.linalg.norm(pos, axis=1)

    r_cut = r_trunc_ratio * a_p_m
    r_grid = np.linspace(1.0e-4 * a_p_m, r_cut, 2000)
    x_s = r_grid / a_p_m
    f_trunc = ((r_trunc_ratio ** 2) / (1.0 + r_trunc_ratio ** 2)) ** 1.5
    rho_star = (
        (3.0 * total_star_mass_kg / (4.0 * np.pi * (a_p_m ** 3) * f_trunc))
        * ((1.0 + x_s * x_s) ** (-2.5))
    )
    m_star_enc = total_star_mass_kg * (((x_s * x_s) / (1.0 + x_s * x_s)) ** 1.5) / f_trunc
    g_tot = G_EFF * m_star_enc / (r_grid ** 2)

    phi_ext_tail = G_EFF * total_star_mass_kg / r_cut
    if r0_m > 0.0:
        g_burk = np.pi * G_EFF * SIGMA_POCKET_KGM2 * burkert_f(r_grid / r0_m)
        g_tot = g_tot + g_burk
        # Exterior Burkert potential drop from r_cut out to infinity (100 r0 numerical + analytic tail)
        r_ext_hi = max(r_cut * 1.01, 100.0 * r0_m)
        r_ext = np.geomspace(r_cut, r_ext_hi, 1000)
        g_burk_ext = np.pi * G_EFF * SIGMA_POCKET_KGM2 * burkert_f(r_ext / r0_m)
        x_end = r_ext_hi / r0_m
        phi_burk_tail = (
            np.pi * G_EFF * SIGMA_POCKET_KGM2 * r0_m * (3.0 * np.log(x_end) + 3.0 - np.pi) / x_end
        )
        phi_ext_tail += float(np.trapezoid(g_burk_ext, r_ext)) + phi_burk_tail

    sigma_r_grid = _jeans_sigma_r(r_grid, rho_star, g_tot)
    g_fwd = cumulative_trapezoid(g_tot, r_grid, initial=0.0)
    phi_drop = np.maximum((g_fwd[-1] - g_fwd) + phi_ext_tail, 0.0)
    v_esc_grid = np.sqrt(2.0 * phi_drop)

    vel = _sample_jeans_velocities(rng, radii_centered, r_grid, sigma_r_grid, v_esc_grid)
    vel -= np.mean(vel, axis=0)

    m_star = total_star_mass_kg / n_stars
    masses = np.full(n_stars, m_star, dtype=np.float64)
    node_weights = masses / M_PROTON

    acc0 = solver.get_particle_accelerations(pos, node_weights, mode="Pure K_Omega")
    if r0_m > 0.0:
        acc0 = acc0 + burkert_acceleration_vec(pos, r0_m)

    origin = np.zeros(3) if r0_m > 0.0 else None
    p_ref = pos if origin is None else (pos - origin)
    if origin is None:
        p_ref = pos - np.mean(pos, axis=0)
    k2 = float(np.sum(masses * np.sum(vel * vel, axis=1)))
    w_abs = abs(float(np.sum(masses * np.sum(p_ref * acc0, axis=1))))
    if k2 > 0.0 and w_abs > 0.0:
        vel *= np.sqrt(w_abs / k2)

    return pos, vel, acc0, masses, node_weights


# ==============================================================================
# 4. Experiment A: 6-Satellite Live-Star Plummer Benchmark
# ==============================================================================
def run_experiment_a_case(task):
    """Worker function executing a single Experiment A satellite run over 50 crossings."""
    t_start = time.perf_counter()
    target = task["target"]
    mode_label = task["mode"]
    r0_pc = float(task["r0_pc"])
    n_stars = int(task["n_stars"])
    grid_n = int(task["grid_n"])
    crossings = int(task["crossings"])
    seed = int(task["seed"])

    m_total_kg = float(task["M_solar"]) * M_SOLAR_KG
    a_p_m = float(task["r_h_pc"]) * PC_TO_METER
    r0_m = r0_pc * PC_TO_METER if r0_pc > 0.0 else 0.0
    box_l = 12.0 * a_p_m

    solver = QUMONDMeshSolver(box_l, grid_n=grid_n, fft_workers=1)
    pos, vel, acc, masses, node_weights = generate_jeans_plummer_ics(
        n_stars, a_p_m, m_total_kg, solver, r0_m=r0_m, r_trunc_ratio=4.5, seed=seed
    )

    origin0 = np.zeros(3) if r0_m > 0.0 else None
    sig0, r_half_0, vir0, v_com0 = _compute_metrics(pos, vel, acc, masses, origin_pos=origin0)
    sigma_3d_0 = np.sqrt(3.0) * sig0
    t_cross = 2.0 * r_half_0 / max(sigma_3d_0, 1.0)
    steps_per_cross = 25
    dt = t_cross / steps_per_cross
    half_dt = 0.5 * dt
    total_steps = crossings * steps_per_cross

    ts_rows = [
        {
            "Target": target,
            "Mode": mode_label,
            "t_cross": 0,
            "sigma_LOS_kms": sig0 / 1000.0,
            "r_half_star_pc": r_half_0 / PC_TO_METER,
            "r_half_star_over_init": 1.0,
            "r_half_pocket_pc": r0_pc if r0_pc > 0.0 else 0.0,
            "r_half_pocket_over_init": 1.0,
            "rho_pocket_core_msun_pc3": 0.0,
            "rho_pocket_core_over_init": 1.0,
            "virial_2K_over_W": vir0,
            "com_drift_kms": 0.0,
        }
    ]

    for step in range(1, total_steps + 1):
        vel += half_dt * acc
        pos += dt * vel
        acc = solver.get_particle_accelerations(pos, node_weights, mode="Pure K_Omega")
        if r0_m > 0.0:
            acc = acc + burkert_acceleration_vec(pos, r0_m)
        vel += half_dt * acc

        if step % steps_per_cross == 0:
            c_idx = step // steps_per_cross
            sig, r_half, vir, v_com = _compute_metrics(pos, vel, acc, masses, origin_pos=origin0)
            com_drift = float(np.linalg.norm(v_com - v_com0)) / 1000.0
            ts_rows.append(
                {
                    "Target": target,
                    "Mode": mode_label,
                    "t_cross": c_idx,
                    "sigma_LOS_kms": sig / 1000.0,
                    "r_half_star_pc": r_half / PC_TO_METER,
                    "r_half_star_over_init": r_half / r_half_0,
                    "r_half_pocket_pc": r0_pc if r0_pc > 0.0 else 0.0,
                    "r_half_pocket_over_init": 1.0,
                    "rho_pocket_core_msun_pc3": 0.0,
                    "rho_pocket_core_over_init": 1.0,
                    "virial_2K_over_W": vir,
                    "com_drift_kms": com_drift,
                }
            )

    runtime_s = time.perf_counter() - t_start
    ts_df = pd.DataFrame(ts_rows)
    late = ts_df[ts_df["t_cross"] >= min(20, crossings)]
    sig_mean = float(late["sigma_LOS_kms"].mean())
    sig_std = float(late["sigma_LOS_kms"].std(ddof=0))
    obs_sig = float(task["obs_sigma_kms"])
    obs_err = float(task["obs_err_kms"])
    max_drift = float(ts_df["com_drift_kms"].max())

    summary = {
        "Target": target,
        "Mode": mode_label,
        "r0_pc": r0_pc,
        "N_stars": n_stars,
        "N_pocket": 0,
        "sigma_LOS_mean_kms": sig_mean,
        "sigma_LOS_std_kms": sig_std,
        "observed_sigma_kms": obs_sig,
        "observed_err_kms": obs_err,
        "z_score": (sig_mean - obs_sig) / obs_err,
        "r_half_final_over_init": float(ts_df["r_half_star_over_init"].iloc[-1]),
        "virial_2K_over_W_mean": float(late["virial_2K_over_W"].mean()),
        "com_drift_kms": max_drift,
        "com_drift_norm": max_drift / max(sig_mean, 1.0e-12),
        "runtime_s": runtime_s,
    }
    return summary, ts_rows


# ==============================================================================
# 5. Experiment B: Live Two-Component Self-Gravitating Draco dSph Proof
# ==============================================================================
def _pocket_radial_profile(pos_pocket_m, com_m, m_pocket_kg, r0_m, r_edges_m):
    """Computes spherically averaged volumetric pocket density profile rho(r) [M_sun / pc^3]."""
    r_p = np.linalg.norm(pos_pocket_m - com_m, axis=1)
    counts, _ = np.histogram(r_p, bins=r_edges_m)
    vol_m3 = (4.0 / 3.0) * np.pi * (r_edges_m[1:] ** 3 - r_edges_m[:-1] ** 3)
    rho_kg_m3 = (counts * m_pocket_kg) / np.maximum(vol_m3, 1.0e-30)
    rho_msun_pc3 = rho_kg_m3 * (PC_TO_METER ** 3) / M_SOLAR_KG
    r_mid_pc = np.sqrt(r_edges_m[:-1] * r_edges_m[1:]) / PC_TO_METER
    return r_mid_pc, rho_msun_pc3


def run_experiment_b_draco_live(n_stars=10000, n_pocket=20000, grid_n=64, crossings=50, seed=1042):
    """
    Live 2-component Draco dSph N-body test (10,000 stars + 20,000 collisionless
    Burkert pocket nodes) under strictly linear Pure K_Omega 3D FFT PM gravity with
    zero artificial COM recentering.
    """
    t_start = time.perf_counter()
    rng = np.random.default_rng(seed)

    r0_pc = 220.0
    r0_m = r0_pc * PC_TO_METER
    r_trunc_ratio = 4.0
    r_cut = r_trunc_ratio * r0_m
    box_l = 12.0 * r0_m

    m_star_tot_kg = 3.0e5 * M_SOLAR_KG
    m_pocket_tot_kg = np.pi * SIGMA_POCKET_KGM2 * (r0_m ** 2) * burkert_m(r_trunc_ratio)
    m_star = m_star_tot_kg / n_stars
    m_pocket = m_pocket_tot_kg / n_pocket

    # Sample both spatial distributions inside r_cut = 4.0 * r_0
    r_stars = _sample_plummer_radii(rng, n_stars, r0_m, r_trunc_ratio)
    r_pocket = _sample_burkert_radii(rng, n_pocket, r0_m, r_trunc_ratio)
    pos_s = _sample_isotropic_positions(rng, r_stars)
    pos_p = _sample_isotropic_positions(rng, r_pocket)
    pos_s -= np.mean(pos_s, axis=0)
    pos_p -= np.mean(pos_p, axis=0)

    r_s_cen = np.linalg.norm(pos_s, axis=1)
    r_p_cen = np.linalg.norm(pos_p, axis=1)

    # Exact self-consistent spherical Jeans integrals on [1e-4 r0, r_cut] in g_tot = g_* + g_pocket
    r_grid = np.linspace(1.0e-4 * r0_m, r_cut, 2000)
    x = r_grid / r0_m
    f_trunc = ((r_trunc_ratio ** 2) / (1.0 + r_trunc_ratio ** 2)) ** 1.5

    rho_s = (3.0 * m_star_tot_kg / (4.0 * np.pi * (r0_m ** 3) * f_trunc)) * ((1.0 + x * x) ** (-2.5))
    m_s_enc = m_star_tot_kg * (((x * x) / (1.0 + x * x)) ** 1.5) / f_trunc
    rho_p = (SIGMA_POCKET_KGM2 / r0_m) / ((1.0 + x) * (1.0 + x * x))
    m_p_enc = np.pi * SIGMA_POCKET_KGM2 * (r0_m ** 2) * burkert_m(x)

    g_tot = G_EFF * (m_s_enc + m_p_enc) / (r_grid ** 2)
    sig_r_s = _jeans_sigma_r(r_grid, rho_s, g_tot)
    sig_r_p = _jeans_sigma_r(r_grid, rho_p, g_tot)

    g_fwd = cumulative_trapezoid(g_tot, r_grid, initial=0.0)
    phi_drop = np.maximum((g_fwd[-1] - g_fwd) + G_EFF * (m_star_tot_kg + m_pocket_tot_kg) / r_cut, 0.0)
    v_esc = np.sqrt(2.0 * phi_drop)

    vel_s = _sample_jeans_velocities(rng, r_s_cen, r_grid, sig_r_s, v_esc)
    vel_p = _sample_jeans_velocities(rng, r_p_cen, r_grid, sig_r_p, v_esc)
    vel_s -= np.mean(vel_s, axis=0)
    vel_p -= np.mean(vel_p, axis=0)

    pos = np.vstack((pos_s, pos_p))
    vel = np.vstack((vel_s, vel_p))
    masses_s = np.full(n_stars, m_star, dtype=np.float64)
    masses_p = np.full(n_pocket, m_pocket, dtype=np.float64)
    masses = np.concatenate((masses_s, masses_p))
    node_weights = masses / M_PROTON
    total_m = float(np.sum(masses))

    solver = QUMONDMeshSolver(box_l, grid_n=grid_n, fft_workers=1)
    acc = solver.get_particle_accelerations(pos, node_weights, mode="Pure K_Omega")

    # Calibrate each species independently to exact initial 2K_species / |W_species| = 1.0000
    k2_s = float(np.sum(masses_s * np.sum(vel[:n_stars] ** 2, axis=1)))
    w_s_abs = abs(float(np.sum(masses_s * np.sum(pos[:n_stars] * acc[:n_stars], axis=1))))
    if k2_s > 0.0 and w_s_abs > 0.0:
        vel[:n_stars] *= np.sqrt(w_s_abs / k2_s)

    k2_p = float(np.sum(masses_p * np.sum(vel[n_stars:] ** 2, axis=1)))
    w_p_abs = abs(float(np.sum(masses_p * np.sum(pos[n_stars:] * acc[n_stars:], axis=1))))
    if k2_p > 0.0 and w_p_abs > 0.0:
        vel[n_stars:] *= np.sqrt(w_p_abs / k2_p)

    # Zero exact initial total momentum once at t=0 (NO recentering after t=0)
    p_tot_0 = np.sum(masses[:, None] * vel, axis=0)
    vel -= (p_tot_0 / total_m)
    p_tot_0 = np.sum(masses[:, None] * vel, axis=0)

    def measure_state(p_arr, v_arr, a_arr):
        com_tot = np.sum(masses[:, None] * p_arr, axis=0) / total_m
        p_tot = np.sum(masses[:, None] * v_arr, axis=0)
        dp_norm_v = float(np.linalg.norm(p_tot - p_tot_0)) / total_m

        sig_s, r_h_s, _, _ = _compute_metrics(
            p_arr[:n_stars], v_arr[:n_stars], a_arr[:n_stars], masses_s, origin_pos=com_tot
        )
        _, r_h_p, _, _ = _compute_metrics(
            p_arr[n_stars:], v_arr[n_stars:], a_arr[n_stars:], masses_p, origin_pos=com_tot
        )
        _, _, vir_tot, _ = _compute_metrics(p_arr, v_arr, a_arr, masses, origin_pos=com_tot)

        r_p_dist = np.linalg.norm(p_arr[n_stars:] - com_tot, axis=1)
        core_r = 0.5 * r0_m
        core_mass = float(np.sum(r_p_dist < core_r)) * m_pocket
        rho_core_kg_m3 = core_mass / ((4.0 / 3.0) * np.pi * (core_r ** 3))
        rho_core_msun_pc3 = rho_core_kg_m3 * (PC_TO_METER ** 3) / M_SOLAR_KG
        return com_tot, dp_norm_v, sig_s, r_h_s, r_h_p, rho_core_msun_pc3, vir_tot

    com0, _, sig_s0, r_hs0, r_hp0, rho_core0, vir0 = measure_state(pos, vel, acc)
    r_edges_m = np.geomspace(0.20 * r0_m, r_cut, 19)
    r_prof_pc, rho_prof_0 = _pocket_radial_profile(pos[n_stars:], com0, m_pocket, r0_m, r_edges_m)

    sigma_3d_s0 = np.sqrt(3.0) * sig_s0
    t_cross = 2.0 * r_hs0 / max(sigma_3d_s0, 1.0)
    steps_per_cross = 25
    dt = t_cross / steps_per_cross
    half_dt = 0.5 * dt
    total_steps = crossings * steps_per_cross

    ts_rows = [
        {
            "Target": "Draco dSph",
            "Mode": "Live 2-Component Pocket",
            "t_cross": 0,
            "sigma_LOS_kms": sig_s0 / 1000.0,
            "r_half_star_pc": r_hs0 / PC_TO_METER,
            "r_half_star_over_init": 1.0,
            "r_half_pocket_pc": r_hp0 / PC_TO_METER,
            "r_half_pocket_over_init": 1.0,
            "rho_pocket_core_msun_pc3": rho_core0,
            "rho_pocket_core_over_init": 1.0,
            "virial_2K_over_W": vir0,
            "com_drift_kms": 0.0,
        }
    ]

    com_final = com0
    for step in range(1, total_steps + 1):
        vel += half_dt * acc
        pos += dt * vel
        acc = solver.get_particle_accelerations(pos, node_weights, mode="Pure K_Omega")
        vel += half_dt * acc

        if step % steps_per_cross == 0:
            c_idx = step // steps_per_cross
            com_t, dp_v, sig_s, r_hs, r_hp, rho_core, vir_t = measure_state(pos, vel, acc)
            com_final = com_t
            ts_rows.append(
                {
                    "Target": "Draco dSph",
                    "Mode": "Live 2-Component Pocket",
                    "t_cross": c_idx,
                    "sigma_LOS_kms": sig_s / 1000.0,
                    "r_half_star_pc": r_hs / PC_TO_METER,
                    "r_half_star_over_init": r_hs / r_hs0,
                    "r_half_pocket_pc": r_hp / PC_TO_METER,
                    "r_half_pocket_over_init": r_hp / r_hp0,
                    "rho_pocket_core_msun_pc3": rho_core,
                    "rho_pocket_core_over_init": rho_core / rho_core0,
                    "virial_2K_over_W": vir_t,
                    "com_drift_kms": dp_v / 1000.0,
                }
            )

    _, rho_prof_final = _pocket_radial_profile(pos[n_stars:], com_final, m_pocket, r0_m, r_edges_m)
    x_prof = (r_prof_pc * PC_TO_METER) / r0_m
    rho_analytic = (
        (SIGMA_POCKET_KGM2 / r0_m) / ((1.0 + x_prof) * (1.0 + x_prof * x_prof))
    ) * (PC_TO_METER ** 3) / M_SOLAR_KG

    runtime_s = time.perf_counter() - t_start
    ts_df = pd.DataFrame(ts_rows)
    late = ts_df[ts_df["t_cross"] >= min(20, crossings)]
    sig_mean = float(late["sigma_LOS_kms"].mean())
    sig_std = float(late["sigma_LOS_kms"].std(ddof=0))
    obs_sig, obs_err = 9.10, 1.20
    max_drift = float(ts_df["com_drift_kms"].max())

    summary = {
        "Target": "Draco dSph",
        "Mode": "Live 2-Component Pocket",
        "r0_pc": r0_pc,
        "N_stars": n_stars,
        "N_pocket": n_pocket,
        "sigma_LOS_mean_kms": sig_mean,
        "sigma_LOS_std_kms": sig_std,
        "observed_sigma_kms": obs_sig,
        "observed_err_kms": obs_err,
        "z_score": (sig_mean - obs_sig) / obs_err,
        "r_half_final_over_init": float(ts_df["r_half_star_over_init"].iloc[-1]),
        "virial_2K_over_W_mean": float(late["virial_2K_over_W"].mean()),
        "com_drift_kms": max_drift,
        "com_drift_norm": max_drift / max(sig_mean, 1.0e-12),
        "runtime_s": runtime_s,
    }
    pocket_stability = {
        "r_half_pocket_final_over_init": float(ts_df["r_half_pocket_over_init"].iloc[-1]),
        "rho_pocket_core_init_msun_pc3": float(rho_core0),
        "rho_pocket_core_final_over_init": float(ts_df["rho_pocket_core_over_init"].iloc[-1]),
        "r_prof_pc": r_prof_pc,
        "rho_prof_0": rho_prof_0,
        "rho_prof_final": rho_prof_final,
        "rho_analytic": rho_analytic,
    }
    return summary, ts_rows, pocket_stability


# ==============================================================================
# 6. Experiment C: Zero-Parameter Spiral & Dwarf Disk Rotation Curve Check
# ==============================================================================
def run_experiment_c_disks():
    """
    Evaluates the zero-parameter mass-locked Burkert Cover Pocket + Freeman thin
    exponential disk rotation curve against the legacy a_Omega underflow formula
    and observed flat velocities across 4 benchmark disk galaxies.
    """
    rows = []
    curves = {}
    for disk in DISK_CATALOG:
        name = disk["name"]
        m_b_solar = float(disk["M_B_solar"])
        m_b_kg = m_b_solar * M_SOLAR_KG
        r_d_kpc = float(disk["R_d_kpc"])
        r_focus = float(disk["r_focus_kpc"])
        r_max = float(disk["r_max_kpc"])
        v_obs = float(disk["v_obs_kms"])

        r0_m = float(np.sqrt(ETA_POCKET * m_b_kg / (np.pi * SIGMA_POCKET_KGM2 * M_FAC_VMAX)))
        r0_kpc = r0_m / KPC_TO_METER
        m_pocket_tot_solar = (
            np.pi * SIGMA_POCKET_KGM2 * (r0_m ** 2) * M_FAC_VMAX
        ) / M_SOLAR_KG

        r_grid = np.linspace(0.05, r_max, 240)
        v_old_aomega, v_kepler, g_freeman, _, _, v_flat_old = calculate_rotation_profile(
            m_b_solar,
            r_grid,
            a_omega_active=A_OMEGA,
            R_d_kpc=r_d_kpc,
            disk_geometry="Freeman Thin Disk",
        )
        r_m = r_grid * KPC_TO_METER
        g_burkert = np.pi * G_EFF * SIGMA_POCKET_KGM2 * burkert_f(r_m / r0_m)
        v_pocket = np.sqrt((g_freeman + g_burkert) * r_m) / 1000.0

        # Also evaluate exact point at r_focus
        v_old_foc, v_kep_foc, g_free_foc, _, _, _ = calculate_rotation_profile(
            m_b_solar,
            np.array([r_focus, r_max], dtype=np.float64),
            a_omega_active=A_OMEGA,
            R_d_kpc=r_d_kpc,
            disk_geometry="Freeman Thin Disk",
        )
        r_pt_m = np.array([r_focus, r_max], dtype=np.float64) * KPC_TO_METER
        g_burk_pt = np.pi * G_EFF * SIGMA_POCKET_KGM2 * burkert_f(r_pt_m / r0_m)
        v_pocket_pt = np.sqrt((g_free_foc + g_burk_pt) * r_pt_m) / 1000.0

        rows.append(
            {
                "Galaxy": name,
                "M_B_solar": m_b_solar,
                "R_d_kpc": r_d_kpc,
                "r0_locked_kpc": r0_kpc,
                "M_pocket_at_vmax_over_M_B": m_pocket_tot_solar / m_b_solar,
                "r_focus_kpc": r_focus,
                "r_max_kpc": r_max,
                "v_obs_kms": v_obs,
                "v_pocket_at_focus_kms": float(v_pocket_pt[0]),
                "v_pocket_at_rmax_kms": float(v_pocket_pt[1]),
                "v_pocket_max_kms": float(np.max(v_pocket)),
                "v_old_aOmega_at_focus_kms": float(v_old_foc[0]),
                "v_old_aOmega_flat_kms": float(v_flat_old),
                "v_kepler_at_focus_kms": float(v_kep_foc[0]),
                "pocket_rel_err_at_focus": (float(v_pocket_pt[0]) - v_obs) / v_obs,
                "a_BTFR_ms2": A_BTFR,
                "A_OMEGA_ms2": A_OMEGA,
                "a_BTFR_over_A_OMEGA": A_BTFR / A_OMEGA,
            }
        )
        curves[name] = {
            "r_kpc": r_grid,
            "v_pocket": v_pocket,
            "v_old": v_old_aomega,
            "v_kepler": v_kepler,
            "v_obs": v_obs,
            "r0_kpc": r0_kpc,
        }
    return pd.DataFrame(rows), curves


# ==============================================================================
# 7. Dark-Theme 4-Panel Diagnostic Figure
# ==============================================================================
def build_diagnostic_figure(summary_df, ts_df, pocket_stab, disk_curves, out_png, out_dir):
    """Generates the dark-theme (#0b1120) 4-panel Cover Pocket diagnostic figure at 150 DPI."""
    fig, axes = plt.subplots(2, 2, figsize=(16.5, 11.5), dpi=150)
    fig.patch.set_facecolor("#0b1120")

    for ax in axes.flat:
        ax.set_facecolor("#070c18")
        ax.tick_params(colors="#cbd5e1", labelsize=9.5)
        ax.xaxis.label.set_color("#e2e8f0")
        ax.yaxis.label.set_color("#e2e8f0")
        ax.title.set_color("#f8fafc")
        for spine in ax.spines.values():
            spine.set_color("#334155")
        ax.grid(True, linestyle="--", alpha=0.22, color="#64748b")

    # Optional load of sibling QUMOND + MW EFE results for Panel 1
    qumond_map = {}
    qumond_csv = os.path.join(REPO_ROOT, "qumond_audit", "qumond_satellite_summary.csv")
    if os.path.exists(qumond_csv):
        try:
            qdf = pd.read_csv(qumond_csv)
            efe_rows = qdf[qdf["Mode"] == "QUMOND + MW EFE"]
            for _, r in efe_rows.iterrows():
                qumond_map[str(r["Target"])] = float(r["sigma_LOS_mean_kms"])
        except Exception:
            qumond_map = {}

    # ------------------------------------------------------------------
    # Panel 1 (Top-Left): 6-Satellite sigma_LOS Comparison
    # ------------------------------------------------------------------
    ax1 = axes[0, 0]
    targets = [s["target"] for s in SATELLITE_CATALOG]
    short_names = ["Pal 14\n(GC)", "NGC 2419\n(GC)", "Segue 1\n(UFD)", "Draco\n(dSph)", "UMi\n(dSph)", "Sculptor\n(dSph)"]
    x_idx = np.arange(len(targets))

    def _lookup_sig(t_name, m_name):
        sub = summary_df[(summary_df["Target"] == t_name) & (summary_df["Mode"] == m_name)]
        return float(sub["sigma_LOS_mean_kms"].iloc[0]) if not sub.empty else np.nan

    pure_vals = [_lookup_sig(t, "Pure K_Omega") for t in targets]
    qum_vals = [qumond_map.get(t, np.nan) for t in targets]
    rh_vals = [_lookup_sig(t, "Cover Pocket (r0=r_h)") for t in targets]
    r300_vals = [_lookup_sig(t, "Cover Pocket (r0=300pc)") for t in targets]
    obs_vals = [s["obs_sigma_kms"] for s in SATELLITE_CATALOG]
    obs_errs = [s["obs_err_kms"] for s in SATELLITE_CATALOG]

    w = 0.20
    ax1.bar(x_idx - 1.5 * w, pure_vals, width=w, color="#64748b", label=r"Pure $\mathcal{K}_\Omega$ (No Pocket)")
    if any(np.isfinite(qum_vals)):
        ax1.bar(x_idx - 0.5 * w, np.nan_to_num(qum_vals, nan=0.0), width=w, color="#f97316", alpha=0.85,
                label="QUMOND + MW EFE (prev.)")
    ax1.bar(x_idx + 0.5 * w, np.nan_to_num(rh_vals, nan=0.0), width=w, color="#38bdf8",
            label=r"Cover Pocket ($r_0 = r_h$)")
    ax1.bar(x_idx + 1.5 * w, np.nan_to_num(r300_vals, nan=0.0), width=w, color="#a78bfa",
            label=r"Cover Pocket ($r_0 = 300\ \mathrm{pc}$)")

    ax1.errorbar(
        x_idx,
        obs_vals,
        yerr=obs_errs,
        fmt="o",
        color="#facc15",
        ecolor="#fde047",
        elinewidth=1.8,
        capsize=4.5,
        markersize=6.5,
        zorder=6,
        label=r"Observed $\sigma_{\mathrm{LOS}} \pm 1\sigma$",
    )
    ax1.set_xticks(x_idx)
    ax1.set_xticklabels(short_names, fontsize=8.8)
    ax1.set_ylabel(r"$\sigma_{\mathrm{LOS}}$ [km/s]")
    ax1.set_ylim(0.0, 14.2)
    ax1.set_title(
        r"Exp A: Zero-Parameter $3P_\delta$ Cover Pocket vs. Pure $\mathcal{K}_\Omega$ Across 6 Satellites",
        fontsize=11.0,
        fontweight="bold",
        pad=10,
    )
    ax1.legend(loc="upper left", fontsize=8.0, frameon=True, facecolor="#0f172a", edgecolor="#334155", labelcolor="#e2e8f0")

    # ------------------------------------------------------------------
    # Panel 2 (Top-Right): Experiment B Draco Live 2-Component Evolution
    # ------------------------------------------------------------------
    ax2 = axes[0, 1]
    draco_live = ts_df[(ts_df["Target"] == "Draco dSph") & (ts_df["Mode"] == "Live 2-Component Pocket")]
    t_arr = draco_live["t_cross"].to_numpy()
    sig_t = draco_live["sigma_LOS_kms"].to_numpy()
    rhs_ratio = draco_live["r_half_star_over_init"].to_numpy()
    rhp_ratio = draco_live["r_half_pocket_over_init"].to_numpy()

    ax2.axhspan(9.10 - 1.20, 9.10 + 1.20, color="#22c55e", alpha=0.14, label=r"Draco Obs: $9.10 \pm 1.20$ km/s")
    ax2.axhline(9.10, color="#22c55e", linestyle=":", linewidth=1.3)
    l1 = ax2.plot(t_arr, sig_t, color="#38bdf8", linewidth=2.2, label=r"Live Stars $\sigma_{\mathrm{LOS}}(t)$ [km/s]")
    ax2.set_xlabel(r"Stellar Crossing Times ($t / t_{\mathrm{cross}}$)")
    ax2.set_ylabel(r"Stellar $\sigma_{\mathrm{LOS}}$ [km/s]", color="#38bdf8")
    ax2.set_ylim(5.5, 13.5)

    ax2_r = ax2.twinx()
    ax2_r.tick_params(colors="#cbd5e1", labelsize=9.0)
    ax2_r.yaxis.label.set_color("#f472b6")
    l2 = ax2_r.plot(t_arr, rhs_ratio, color="#f472b6", linewidth=1.8, linestyle="-", label=r"Stars $r_{1/2}(t)/r_{1/2}(0)$")
    l3 = ax2_r.plot(t_arr, rhp_ratio, color="#facc15", linewidth=1.8, linestyle="--", label=r"Pocket $r_{1/2}(t)/r_{1/2}(0)$")
    ax2_r.axhline(1.0, color="#94a3b8", linestyle=":", linewidth=1.0, alpha=0.6)
    ax2_r.set_ylabel(r"Half-Mass Radius Ratio $r_{1/2}(t) / r_{1/2}(0)$")
    ax2_r.set_ylim(0.75, 1.25)
    for sp in ax2_r.spines.values():
        sp.set_color("#334155")

    lines_2 = l1 + l2 + l3
    labels_2 = [ln.get_label() for ln in lines_2]
    ax2.legend(lines_2, labels_2, loc="upper right", fontsize=8.0, frameon=True, facecolor="#0f172a", edgecolor="#334155", labelcolor="#e2e8f0")
    ax2.set_title(
        r"Exp B: Draco Live 2-Component ($10\mathrm{k}$ Stars + $20\mathrm{k}$ Pocket Nodes)",
        fontsize=11.0,
        fontweight="bold",
        pad=10,
    )

    # ------------------------------------------------------------------
    # Panel 3 (Bottom-Left): Experiment B Live Pocket Radial Profile + Momentum Badge
    # ------------------------------------------------------------------
    ax3 = axes[1, 0]
    r_pc = pocket_stab["r_prof_pc"]
    ax3.loglog(r_pc, pocket_stab["rho_analytic"], color="#94a3b8", linestyle="--", linewidth=1.8,
               label=rf"Analytic Burkert ($\Sigma_p = {SIGMA_POCKET_MSUN_PC2:.1f}\ M_\odot/\mathrm{{pc}}^2$)")
    ax3.loglog(r_pc, pocket_stab["rho_prof_0"], "o-", color="#38bdf8", linewidth=2.0, markersize=4.5,
               label=r"Live Pocket $t = 0$")
    ax3.loglog(r_pc, pocket_stab["rho_prof_final"], "s-", color="#22c55e", linewidth=2.0, markersize=4.5,
               label=rf"Live Pocket $t = {int(t_arr[-1])}\,t_{{\mathrm{{cross}}}}$")

    # NFW r^-1 inner cusp reference slope anchored at outer point for visual contrast
    r_ref_nfw = r_pc[r_pc <= 260.0]
    if len(r_ref_nfw) >= 2:
        rho_anchor = pocket_stab["rho_analytic"][len(r_ref_nfw) - 1]
        rho_nfw_slope = rho_anchor * (r_ref_nfw[-1] / r_ref_nfw)
        ax3.loglog(r_ref_nfw, rho_nfw_slope, ":", color="#f97316", linewidth=1.5, label=r"NFW $r^{-1}$ Cusp Guide")

    draco_sum = summary_df[(summary_df["Target"] == "Draco dSph") & (summary_df["Mode"] == "Live 2-Component Pocket")].iloc[0]
    dp_norm = float(draco_sum["com_drift_norm"])
    rho_ratio = float(pocket_stab["rho_pocket_core_final_over_init"])
    rhp_fin = float(pocket_stab["r_half_pocket_final_over_init"])
    badge_txt = (
        rf"$|\Delta P_{{\mathrm{{tot}}}}| / (M_{{\mathrm{{tot}}}}\sigma_\star) = {dp_norm:.2e}$ (IEEE-754)" "\n"
        rf"Core $\rho(<0.5\,r_0)_{{50}} / \rho_0 = {rho_ratio:.3f}$" "\n"
        rf"Pocket $r_{{1/2,50}} / r_{{1/2,0}} = {rhp_fin:.3f}$"
    )
    ax3.text(
        0.04,
        0.07,
        badge_txt,
        transform=ax3.transAxes,
        fontsize=8.5,
        color="#f8fafc",
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#0f172a", edgecolor="#22c55e", alpha=0.92),
    )
    ax3.set_xlabel(r"Radius $r$ [pc]")
    ax3.set_ylabel(r"Pocket Density $\rho_{\mathrm{pocket}}(r)$ [$M_\odot/\mathrm{pc}^3$]")
    ax3.set_title(
        r"Exp B: Cored Burkert Stability & Exact Momentum Conservation",
        fontsize=11.0,
        fontweight="bold",
        pad=10,
    )
    ax3.legend(loc="upper right", fontsize=8.0, frameon=True, facecolor="#0f172a", edgecolor="#334155", labelcolor="#e2e8f0")

    # ------------------------------------------------------------------
    # Panel 4 (Bottom-Right): Experiment C Disk Rotation Curves (DDO 154 & NGC 3198)
    # ------------------------------------------------------------------
    ax4 = axes[1, 1]
    ddo = disk_curves["DDO 154"]
    ngc = disk_curves["NGC 3198"]

    ax4.plot(ngc["r_kpc"], ngc["v_pocket"], color="#38bdf8", linewidth=2.3,
             label=rf"NGC 3198 Cover Pocket ($r_0={ngc['r0_kpc']:.2f}$ kpc)")
    ax4.plot(ngc["r_kpc"], ngc["v_old"], color="#38bdf8", linewidth=1.6, linestyle="--", alpha=0.75,
             label=r"NGC 3198 Old $a_\Omega$ Underflow")
    ax4.plot(ngc["r_kpc"], ngc["v_kepler"], color="#38bdf8", linewidth=1.2, linestyle=":", alpha=0.55,
             label="NGC 3198 Freeman Baryons Only")
    ax4.axhline(ngc["v_obs"], color="#bae6fd", linestyle="-.", linewidth=1.1, alpha=0.7,
                label=rf"NGC 3198 Obs $v_{{\mathrm{{flat}}}} = {ngc['v_obs']:.0f}$ km/s")

    ax4.plot(ddo["r_kpc"], ddo["v_pocket"], color="#f472b6", linewidth=2.3,
             label=rf"DDO 154 Cover Pocket ($r_0={ddo['r0_kpc']:.2f}$ kpc)")
    ax4.plot(ddo["r_kpc"], ddo["v_old"], color="#f472b6", linewidth=1.6, linestyle="--", alpha=0.75,
             label=r"DDO 154 Old $a_\Omega$ Underflow")
    ax4.axhline(ddo["v_obs"], color="#fbcfe8", linestyle="-.", linewidth=1.1, alpha=0.7,
                label=rf"DDO 154 Obs $v_{{\mathrm{{flat}}}} = {ddo['v_obs']:.0f}$ km/s")

    btfr_txt = (
        rf"$a_{{\mathrm{{BTFR}}}} = \pi G_{{\mathrm{{eff}}}} \Sigma_p \eta_{{\mathrm{{pocket}}}} (v^4_{{\mathrm{{max}}}}/m_{{\mathrm{{max}}}}) = {A_BTFR:.3e}\ \mathrm{{m/s}}^2$" "\n"
        rf"$A_\Omega = {A_OMEGA:.3e}\ \mathrm{{m/s}}^2 \quad (a_{{\mathrm{{BTFR}}}}/A_\Omega = {A_BTFR/A_OMEGA:.4f})$"
    )
    ax4.text(
        0.04,
        0.47,
        btfr_txt,
        transform=ax4.transAxes,
        fontsize=7.9,
        color="#f8fafc",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#0f172a", edgecolor="#38bdf8", alpha=0.90),
    )
    ax4.set_xlabel("Galactocentric Radius $r$ [kpc]")
    ax4.set_ylabel("Circular Velocity $v_c(r)$ [km/s]")
    ax4.set_xlim(0.0, 35.0)
    ax4.set_ylim(0.0, 200.0)
    ax4.set_title(
        r"Exp C: Mass-Locked Cover Pocket vs. Old $a_\Omega$ Underflow (DDO 154 & NGC 3198)",
        fontsize=11.0,
        fontweight="bold",
        pad=10,
    )
    ax4.legend(loc="lower right", fontsize=7.6, frameon=True, facecolor="#0f172a", edgecolor="#334155", labelcolor="#e2e8f0")

    fig.suptitle(
        r"Topological Cover Pocket Architecture — Zero-Parameter $3P_\delta$ Burkert + Strictly Linear Pure $\mathcal{K}_\Omega$ N-Body Verification",
        fontsize=13.5,
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
        description="Zero-Parameter Topological Cover Pocket 3D FFT Particle-Mesh N-Body Testbed"
    )
    parser.add_argument("--particles", type=int, default=10000, help="Number of live star particles")
    parser.add_argument(
        "--pocket-particles",
        type=int,
        default=20000,
        help="Number of live Burkert Cover Pocket particles in Experiment B",
    )
    parser.add_argument("--grid", type=int, default=64, help="3D PM grid resolution per axis")
    parser.add_argument("--crossings", type=int, default=50, help="Crossing times per N-body run")
    parser.add_argument("--out-dir", type=str, default="cover_pocket_audit", help="Output directory")
    parser.add_argument(
        "--jobs",
        type=int,
        default=min(6, os.cpu_count() or 4),
        help="Parallel worker processes for Experiment A + B",
    )
    args = parser.parse_args()

    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    print("=== Zero-Parameter Cover Pocket First-Principles Invariants ===")
    print(f"G_EFF                 = {G_EFF:.6e} m^3/(kg s^2)")
    print(f"DELTA_OMEGA_M         = {DELTA_OMEGA_M:.6f}")
    print(f"ETA_POCKET            = {ETA_POCKET:.6f}")
    print(f"G_3PDELTA             = {G_3PDELTA:.6e} m/s^2")
    print(f"Burkert peak grad     : x_gmax={X_GMAX:.4f}, f_B_max={F_B_MAX:.6f}")
    print(f"Burkert peak v_c^2    : x_vmax={X_VMAX:.4f}, v2_B_max={V2_B_MAX:.5f}, m_fac={M_FAC_VMAX:.4f}")
    print(f"SIGMA_POCKET          = {SIGMA_POCKET_KGM2:.6f} kg/m^2 = {SIGMA_POCKET_MSUN_PC2:.2f} M_sun/pc^2")
    print(f"Emergent A_BTFR       = {A_BTFR:.6e} m/s^2  (vs A_OMEGA = {A_OMEGA:.6e}, ratio={A_BTFR/A_OMEGA:.4f})")

    # Build Experiment A tasks
    exp_a_tasks = []
    seed_counter = 100
    for sat in SATELLITE_CATALOG:
        for mode_label, r0_pc in sat["r0_modes"]:
            exp_a_tasks.append(
                {
                    "target": sat["target"],
                    "M_solar": sat["M_solar"],
                    "r_h_pc": sat["r_h_pc"],
                    "obs_sigma_kms": sat["obs_sigma_kms"],
                    "obs_err_kms": sat["obs_err_kms"],
                    "mode": mode_label,
                    "r0_pc": r0_pc,
                    "n_stars": args.particles,
                    "grid_n": args.grid,
                    "crossings": args.crossings,
                    "seed": seed_counter,
                }
            )
            seed_counter += 1

    print(
        f"\nLaunching {len(exp_a_tasks)} Experiment A runs (N_stars={args.particles}) + "
        f"1 Experiment B Live 2-Component Draco run (N={args.particles}+{args.pocket_particles}) "
        f"across {args.jobs} workers..."
    )

    t0_all = time.perf_counter()
    summary_rows = []
    ts_all_rows = []
    with ProcessPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        fut_b = pool.submit(
            run_experiment_b_draco_live,
            args.particles,
            args.pocket_particles,
            args.grid,
            args.crossings,
            1042,
        )
        fut_a_list = [pool.submit(run_experiment_a_case, t) for t in exp_a_tasks]

        for fut in fut_a_list:
            s_row, ts_rows = fut.result()
            summary_rows.append(s_row)
            ts_all_rows.extend(ts_rows)
            print(
                f"  [Exp A] {s_row['Target']:16s} | {s_row['Mode']:23s} | "
                f"sig={s_row['sigma_LOS_mean_kms']:6.2f} km/s (z={s_row['z_score']:+5.2f}) | "
                f"rh50/rh0={s_row['r_half_final_over_init']:.3f} | "
                f"2K/W={s_row['virial_2K_over_W_mean']:.3f} | "
                f"t={s_row['runtime_s']:.1f}s"
            )

        b_sum, b_ts, pocket_stab = fut_b.result()
        summary_rows.append(b_sum)
        ts_all_rows.extend(b_ts)
        print(
            f"  [Exp B] {b_sum['Target']:16s} | {b_sum['Mode']:23s} | "
            f"sig={b_sum['sigma_LOS_mean_kms']:6.2f} km/s (z={b_sum['z_score']:+5.2f}) | "
            f"dP/Mv={b_sum['com_drift_norm']:.2e} | "
            f"rh_p50/rh_p0={pocket_stab['r_half_pocket_final_over_init']:.3f} | "
            f"rho_c50/rho_c0={pocket_stab['rho_pocket_core_final_over_init']:.3f} | "
            f"t={b_sum['runtime_s']:.1f}s"
        )

    # Experiment C: Disk Rotation Curve Cross-Check
    disk_df, disk_curves = run_experiment_c_disks()

    summary_df = pd.DataFrame(summary_rows)
    ts_df = pd.DataFrame(ts_all_rows)

    sum_csv = os.path.join(out_dir, "cover_pocket_satellite_summary.csv")
    ts_csv = os.path.join(out_dir, "cover_pocket_time_series.csv")
    disk_csv = os.path.join(out_dir, "cover_pocket_disk_check.csv")
    png_path = os.path.join(out_dir, "cover_pocket_nbody_benchmark.png")

    summary_df.to_csv(sum_csv, index=False)
    ts_df.to_csv(ts_csv, index=False)
    disk_df.to_csv(disk_csv, index=False)
    build_diagnostic_figure(summary_df, ts_df, pocket_stab, disk_curves, png_path, out_dir)

    print(f"\nTotal benchmark wall time: {time.perf_counter() - t0_all:.1f}s")
    print(f"Saved: {sum_csv}")
    print(f"Saved: {ts_csv}")
    print(f"Saved: {disk_csv}")
    print(f"Saved: {png_path}")


if __name__ == "__main__":
    main()
