"""
Standalone Momentum-Conserving 3D FFT QUMOND Particle-Mesh N-Body Testbed.

Evolves equilibrium 10,000-star Plummer spheres for 50+ crossing times using exact
first-principles Substrate constants from `modules.data_loader` (zero CODATA G) across
three dynamical regimes:
  1. Pure K_Omega (linear 7-point FFT Poisson solver)
  2. Isolated QUMOND (unscreened a_Omega via 4-FFT discrete Milgrom QUMOND)
  3. Total-Field QUMOND + MW EFE (free-falling frame with external host field)

Also executes a 500-star unequal-mass asymmetric momentum-conservation proof contrasting
3D FFT QUMOND against naive pairwise g_eff self-acceleration.
"""

import os
import sys
import time
import argparse
import warnings
import logging
from concurrent.futures import ProcessPoolExecutor

# Ensure repository root is importable regardless of launch directory
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

for logger_name in ("streamlit", "streamlit.runtime.caching", "streamlit.runtime.scriptrunner_utils"):
    logging.getLogger(logger_name).setLevel(logging.ERROR)
warnings.filterwarnings("ignore", message=".*No runtime found.*")
warnings.filterwarnings("ignore", message=".*to view this Streamlit app.*")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.fft import rfftn, irfftn, fftfreq, rfftfreq

from modules.data_loader import (
    K_OMEGA,
    M_PROTON,
    M_SOLAR_KG,
    A_OMEGA,
    KPC_TO_METER,
)

PC_TO_METER = KPC_TO_METER / 1000.0

# Built-in Satellite Benchmark Targets (Stellar Mass, Plummer r_h, Galactocentric D_GC, Obs sigma_LOS)
SATELLITE_TARGETS = {
    "palomar14": {
        "name": "Palomar 14",
        "M_solar": 1.0e4,
        "r_h_pc": 25.0,
        "D_GC_kpc": 71.0,
        "obs_sigma_kms": 0.38,
        "obs_err_kms": 0.12,
        "ref": "Jordi+2009",
    },
    "ngc2419": {
        "name": "NGC 2419",
        "M_solar": 9.0e5,
        "r_h_pc": 20.0,
        "D_GC_kpc": 87.5,
        "obs_sigma_kms": 4.14,
        "obs_err_kms": 0.48,
        "ref": "Baumgardt+2009",
    },
    "draco": {
        "name": "Draco dSph",
        "M_solar": 3.0e5,
        "r_h_pc": 220.0,
        "D_GC_kpc": 76.0,
        "obs_sigma_kms": 9.10,
        "obs_err_kms": 1.20,
        "ref": "Walker+2007",
    },
}

MODES = ("Pure K_Omega", "Isolated QUMOND", "QUMOND + MW EFE")


class QUMONDMeshSolver:
    """
    Zero-CODATA-G, momentum-conserving 3D Particle-Mesh (CIC + 7-point FFT) QUMOND solver.
    Uses identical trilinear CIC weights for node deposition and force interpolation,
    with a pre-allocated periodic halo buffer for zero-allocation central differences.
    """

    def __init__(self, box_length_m, grid_n=64, fft_workers=1):
        self.L = float(box_length_m)
        self.ng = int(grid_n)
        self.dx = self.L / self.ng
        self.inv2dx = 0.5 / self.dx
        self.inv_dx3 = 1.0 / (self.dx ** 3)
        self.fft_workers = fft_workers

        ng = self.ng
        self.ng3 = ng * ng * ng
        self.buf = np.empty((ng + 2, ng + 2, ng + 2), dtype=np.float64)

        kx = 2.0 * np.pi * fftfreq(ng, d=self.dx)
        ky = 2.0 * np.pi * fftfreq(ng, d=self.dx)
        kz = 2.0 * np.pi * rfftfreq(ng, d=self.dx)
        KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")

        lap7 = (2.0 / (self.dx * self.dx)) * (
            3.0 - np.cos(KX * self.dx) - np.cos(KY * self.dx) - np.cos(KZ * self.dx)
        )
        self.G7 = np.zeros_like(lap7)
        mask = lap7 > 0.0
        self.G7[mask] = -1.0 / lap7[mask]

    def _fill_halo(self, arr):
        """Copies 3D periodic boundary faces into pre-allocated (ng+2)^3 halo buffer."""
        b = self.buf
        b[1:-1, 1:-1, 1:-1] = arr
        b[0, 1:-1, 1:-1] = arr[-1, :, :]
        b[-1, 1:-1, 1:-1] = arr[0, :, :]
        b[1:-1, 0, 1:-1] = arr[:, -1, :]
        b[1:-1, -1, 1:-1] = arr[:, 0, :]
        b[1:-1, 1:-1, 0] = arr[:, :, -1]
        b[1:-1, 1:-1, -1] = arr[:, :, 0]
        return b

    def _grad_cd(self, phi):
        """Computes (-dPhi/dx, -dPhi/dy, -dPhi/dz) using central differences."""
        b = self._fill_halo(phi)
        s = self.inv2dx
        gx = -(b[2:, 1:-1, 1:-1] - b[:-2, 1:-1, 1:-1]) * s
        gy = -(b[1:-1, 2:, 1:-1] - b[1:-1, :-2, 1:-1]) * s
        gz = -(b[1:-1, 1:-1, 2:] - b[1:-1, 1:-1, :-2]) * s
        return gx, gy, gz

    def _div_cd(self, fx, fy, fz):
        """Computes central-difference divergence dFx/dx + dFy/dy + dFz/dz."""
        s = self.inv2dx
        b = self._fill_halo(fx)
        div_f = (b[2:, 1:-1, 1:-1] - b[:-2, 1:-1, 1:-1]) * s
        b = self._fill_halo(fy)
        div_f += (b[1:-1, 2:, 1:-1] - b[1:-1, :-2, 1:-1]) * s
        b = self._fill_halo(fz)
        div_f += (b[1:-1, 1:-1, 2:] - b[1:-1, 1:-1, :-2]) * s
        return div_f

    def cic_indices_and_weights(self, pos):
        """
        Computes trilinear CIC flat cell indices and weights once per position array,
        wrapping coordinates modulo L inside [0, ng) so particles can evolve without wrapping pos.
        """
        u = np.mod((pos + 0.5 * self.L) / self.dx, self.ng)
        i0 = np.floor(u).astype(np.int64)
        f = u - i0
        i0 %= self.ng
        i1 = (i0 + 1) % self.ng

        wx0, wx1 = 1.0 - f[:, 0], f[:, 0]
        wy0, wy1 = 1.0 - f[:, 1], f[:, 1]
        wz0, wz1 = 1.0 - f[:, 2], f[:, 2]

        ng = self.ng
        ng2 = ng * ng
        ix0, ix1 = i0[:, 0] * ng2, i1[:, 0] * ng2
        iy0, iy1 = i0[:, 1] * ng, i1[:, 1] * ng
        iz0, iz1 = i0[:, 2], i1[:, 2]

        w000 = wx0 * wy0 * wz0
        w001 = wx0 * wy0 * wz1
        w010 = wx0 * wy1 * wz0
        w011 = wx0 * wy1 * wz1
        w100 = wx1 * wy0 * wz0
        w101 = wx1 * wy0 * wz1
        w110 = wx1 * wy1 * wz0
        w111 = wx1 * wy1 * wz1

        c000 = ix0 + iy0 + iz0
        c001 = ix0 + iy0 + iz1
        c010 = ix0 + iy1 + iz0
        c011 = ix0 + iy1 + iz1
        c100 = ix1 + iy0 + iz0
        c101 = ix1 + iy0 + iz1
        c110 = ix1 + iy1 + iz0
        c111 = ix1 + iy1 + iz1

        return (c000, c001, c010, c011, c100, c101, c110, c111), (
            w000, w001, w010, w011, w100, w101, w110, w111
        )

    def deposit_cic(self, node_weights, cic_data):
        """Deposits node counts N_i onto the 3D density grid rho_N [nodes/m^3]."""
        cells, weights = cic_data
        n3 = self.ng3
        grid = np.bincount(cells[0], weights=node_weights * weights[0], minlength=n3)
        for k in range(1, 8):
            grid += np.bincount(cells[k], weights=node_weights * weights[k], minlength=n3)
        return (grid * self.inv_dx3).reshape(self.ng, self.ng, self.ng)

    def interpolate_cic(self, ax, ay, az, cic_data):
        """Interpolates 3D acceleration grid back onto particle coordinates with identical CIC weights."""
        cells, weights = cic_data
        ax_f = ax.ravel()
        ay_f = ay.ravel()
        az_f = az.ravel()

        p_ax = ax_f[cells[0]] * weights[0]
        p_ay = ay_f[cells[0]] * weights[0]
        p_az = az_f[cells[0]] * weights[0]
        for k in range(1, 8):
            w = weights[k]
            p_ax += ax_f[cells[k]] * w
            p_ay += ay_f[cells[k]] * w
            p_az += az_f[cells[k]] * w
        return np.column_stack((p_ax, p_ay, p_az))

    def solve_grid_accelerations(self, rho_N, mode="Isolated QUMOND", g_ext_N_x=0.0):
        """Solves for the 3D acceleration field on the grid under Pure K_Omega or QUMOND."""
        ng_shape = (self.ng, self.ng, self.ng)
        phi_N_k = (4.0 * np.pi * K_OMEGA) * rfftn(rho_N, workers=self.fft_workers) * self.G7
        phi_N = irfftn(phi_N_k, s=ng_shape, workers=self.fft_workers)

        gx, gy, gz = self._grad_cd(phi_N)
        if mode == "Pure K_Omega":
            return gx, gy, gz

        if mode == "QUMOND + MW EFE" and g_ext_N_x != 0.0:
            gx_tot = gx + g_ext_N_x
        else:
            gx_tot = gx

        g_tot_mag = np.sqrt(gx_tot * gx_tot + gy * gy + gz * gz)
        nu = np.sqrt(1.0 + 1.0 / np.maximum(g_tot_mag / A_OMEGA, 1.0e-12))

        if mode == "QUMOND + MW EFE" and g_ext_N_x != 0.0:
            nu_ext = np.sqrt(1.0 + 1.0 / max(abs(g_ext_N_x) / A_OMEGA, 1.0e-12))
            fx = nu * gx_tot - nu_ext * g_ext_N_x
        else:
            fx = nu * gx_tot
        fy = nu * gy
        fz = nu * gz

        div_f = self._div_cd(fx, fy, fz)
        phi_k = (-rfftn(div_f, workers=self.fft_workers)) * self.G7
        phi = irfftn(phi_k, s=ng_shape, workers=self.fft_workers)
        return self._grad_cd(phi)

    def get_particle_accelerations(self, pos, node_weights, mode="Isolated QUMOND", g_ext_N_x=0.0):
        """Full CIC deposit -> FFT QUMOND solve -> CIC interpolation pass."""
        cic_data = self.cic_indices_and_weights(pos)
        rho_N = self.deposit_cic(node_weights, cic_data)
        ax, ay, az = self.solve_grid_accelerations(rho_N, mode=mode, g_ext_N_x=g_ext_N_x)
        return self.interpolate_cic(ax, ay, az, cic_data)


def generate_plummer_ics(n_stars, scale_a_m, total_mass_kg, solver, mode, g_ext_N_x=0.0, seed=42):
    """
    Generates an isotropic Plummer sphere (truncated at 4.5 * a_p inside L = 12 * a_p)
    using exact Aarseth, Henon & Wielen (1974) velocity rejection sampling, then scales
    velocities by the exact initial virial ratio 2K / |W| = 1.0 from the active field solver.
    """
    rng = np.random.default_rng(seed)

    # Radial inverse CDF truncated cleanly at r_max = 4.5 * a_p
    r_trunc_ratio = 4.5
    u_max = (r_trunc_ratio ** 2 / (1.0 + r_trunc_ratio ** 2)) ** 1.5
    u = rng.uniform(0.0, u_max, size=n_stars)
    radii = scale_a_m / np.sqrt(u ** (-2.0 / 3.0) - 1.0)

    cos_t = rng.uniform(-1.0, 1.0, size=n_stars)
    phi = rng.uniform(0.0, 2.0 * np.pi, size=n_stars)
    sin_t = np.sqrt(1.0 - cos_t * cos_t)
    pos = np.column_stack((radii * sin_t * np.cos(phi), radii * sin_t * np.sin(phi), radii * cos_t))

    # Aarseth, Henon & Wielen (1974) rejection sampling for q in [0, 1], g(q) = q^2 * (1 - q^2)^3.5
    q_samples = np.empty(n_stars, dtype=np.float64)
    filled = 0
    while filled < n_stars:
        batch_n = max(1024, (n_stars - filled) * 12)
        q_try = rng.uniform(0.0, 1.0, size=batch_n)
        y_try = rng.uniform(0.0, 0.1, size=batch_n)
        acc = q_try[y_try < (q_try ** 2) * ((1.0 - q_try ** 2) ** 3.5)]
        take = min(len(acc), n_stars - filled)
        q_samples[filled : filled + take] = acc[:take]
        filled += take

    v_mag = q_samples * ((1.0 + (radii / scale_a_m) ** 2) ** (-0.25))
    cos_tv = rng.uniform(-1.0, 1.0, size=n_stars)
    phiv = rng.uniform(0.0, 2.0 * np.pi, size=n_stars)
    sin_tv = np.sqrt(1.0 - cos_tv * cos_tv)
    vel = v_mag[:, None] * np.column_stack((sin_tv * np.cos(phiv), sin_tv * np.sin(phiv), cos_tv))

    # Zero exact initial COM position & velocity
    pos -= np.mean(pos, axis=0)
    vel -= np.mean(vel, axis=0)

    m_star = total_mass_kg / n_stars
    node_weights = np.full(n_stars, m_star / M_PROTON, dtype=np.float64)

    # Virial equilibrium scaling against the active QUMOND / K_Omega field
    acc0 = solver.get_particle_accelerations(pos, node_weights, mode=mode, g_ext_N_x=g_ext_N_x)
    W0 = float(np.sum(m_star * np.sum(pos * acc0, axis=1)))
    K0 = 0.5 * m_star * float(np.sum(vel * vel))
    scale_v = np.sqrt(abs(W0) / max(2.0 * K0, 1.0e-30))
    vel *= scale_v
    vel -= np.mean(vel, axis=0)

    return pos, vel, acc0, m_star, node_weights, abs(W0)


def _compute_diagnostics(pos, vel, acc, m_star, total_mass_kg, a_p, r_half_init=None):
    """Computes COM drift (uncorrected), r_half, bound sigma_LOS [km/s], and 2K/|W|."""
    com_pos = np.mean(pos, axis=0)
    com_vel = np.mean(vel, axis=0)
    com_drift_kms = float(np.linalg.norm(com_vel)) / 1000.0

    pos_rel = pos - com_pos
    vel_rel = vel - com_vel
    r_rel = np.linalg.norm(pos_rel, axis=1)
    r_half = float(np.median(r_rel))

    bound = r_rel <= 5.5 * a_p
    vel_b = vel_rel[bound] if np.sum(bound) >= 20 else vel_rel
    sigma_3d_sq = float(np.mean(np.sum(vel_b * vel_b, axis=1)))
    sigma_los_kms = np.sqrt(sigma_3d_sq / 3.0) / 1000.0

    W = float(np.sum(m_star * np.sum(pos_rel * acc, axis=1)))
    K = 0.5 * m_star * float(np.sum(vel_rel * vel_rel))
    virial = (2.0 * K) / max(abs(W), 1.0e-30)

    r_half_ratio = 1.0 if r_half_init is None else (r_half / r_half_init)
    return com_pos, com_vel, pos_rel, vel_rel, sigma_los_kms, r_half, r_half_ratio, virial, com_drift_kms


def run_single_simulation(task_args):
    """Runs one (target, mode) N-body integration over 50+ crossing times."""
    target_key, mode, n_stars, grid_n, crossings, mw_mass_solar, fft_workers = task_args
    cfg = SATELLITE_TARGETS[target_key]

    t_start = time.perf_counter()
    a_p = cfg["r_h_pc"] * PC_TO_METER
    box_L = 12.0 * a_p
    total_mass_kg = cfg["M_solar"] * M_SOLAR_KG

    mu_mw = (mw_mass_solar * M_SOLAR_KG / M_PROTON) * K_OMEGA
    g_ext_N_x = mu_mw / ((cfg["D_GC_kpc"] * KPC_TO_METER) ** 2) if mode == "QUMOND + MW EFE" else 0.0

    solver = QUMONDMeshSolver(box_L, grid_n=grid_n, fft_workers=fft_workers)
    pos, vel, acc, m_star, node_weights, abs_W0 = generate_plummer_ics(
        n_stars, a_p, total_mass_kg, solver, mode, g_ext_N_x=g_ext_N_x, seed=42
    )

    sigma_3d_0 = np.sqrt(abs_W0 / total_mass_kg)
    t_cross = 2.0 * a_p / sigma_3d_0
    steps_per_cross = 25
    dt = t_cross / steps_per_cross
    half_dt = 0.5 * dt
    total_steps = crossings * steps_per_cross

    # Diagnostics at t = 0
    _, _, pos_rel0, _, sig0, r_half_init, _, vir0, drift0 = _compute_diagnostics(
        pos, vel, acc, m_star, total_mass_kg, a_p
    )

    ts_rows = [{
        "Target": cfg["name"],
        "Mode": mode,
        "crossing_time": 0,
        "sigma_LOS_kms": round(sig0, 5),
        "r_half_over_init": 1.0,
        "virial_2K_over_W": round(vir0, 5),
        "com_drift_kms": drift0,
    }]

    # Radial projected profile accumulators over t in [20, crossings] t_cross
    r_bin_edges = np.array([0.05, 0.25, 0.45, 0.70, 1.00, 1.35, 1.75, 2.20, 2.75, 3.40], dtype=np.float64)
    r_bin_centers = 0.5 * (r_bin_edges[:-1] + r_bin_edges[1:])
    rad_v2_sum = np.zeros(len(r_bin_centers), dtype=np.float64)
    rad_counts = np.zeros(len(r_bin_centers), dtype=np.int64)

    sigma_window = []
    virial_window = []

    for step in range(1, total_steps + 1):
        # Symplectic KDK Velocity Verlet
        vel += half_dt * acc
        pos += dt * vel
        acc = solver.get_particle_accelerations(pos, node_weights, mode=mode, g_ext_N_x=g_ext_N_x)
        vel += half_dt * acc

        if step % steps_per_cross == 0:
            t_cr = step // steps_per_cross
            _, _, pos_rel, vel_rel, sig_kms, _, rh_ratio, vir_val, drift_kms = _compute_diagnostics(
                pos, vel, acc, m_star, total_mass_kg, a_p, r_half_init=r_half_init
            )
            ts_rows.append({
                "Target": cfg["name"],
                "Mode": mode,
                "crossing_time": t_cr,
                "sigma_LOS_kms": round(sig_kms, 5),
                "r_half_over_init": round(rh_ratio, 5),
                "virial_2K_over_W": round(vir_val, 5),
                "com_drift_kms": drift_kms,
            })

            if t_cr >= min(20, crossings):
                sigma_window.append(sig_kms)
                virial_window.append(vir_val)

                # Accumulate 3 orthogonal line-of-sight projections for radial sigma_LOS(R/r_h)
                rx, ry, rz = pos_rel[:, 0], pos_rel[:, 1], pos_rel[:, 2]
                vx, vy, vz = vel_rel[:, 0], vel_rel[:, 1], vel_rel[:, 2]
                proj_pairs = (
                    (np.hypot(rx, ry) / a_p, vz),
                    (np.hypot(ry, rz) / a_p, vx),
                    (np.hypot(rz, rx) / a_p, vy),
                )
                for r_proj_norm, v_los in proj_pairs:
                    b_idx = np.searchsorted(r_bin_edges, r_proj_norm, side="right") - 1
                    valid_b = (b_idx >= 0) & (b_idx < len(r_bin_centers))
                    rad_counts += np.bincount(b_idx[valid_b], minlength=len(r_bin_centers))
                    rad_v2_sum += np.bincount(b_idx[valid_b], weights=(v_los[valid_b] ** 2), minlength=len(r_bin_centers))

    runtime_s = time.perf_counter() - t_start
    sigma_mean = float(np.mean(sigma_window)) if sigma_window else float(ts_rows[-1]["sigma_LOS_kms"])
    sigma_std = float(np.std(sigma_window)) if sigma_window else 0.0
    vir_mean = float(np.mean(virial_window)) if virial_window else float(ts_rows[-1]["virial_2K_over_W"])
    final_rh_ratio = float(ts_rows[-1]["r_half_over_init"])
    final_drift_kms = float(ts_rows[-1]["com_drift_kms"])
    com_drift_norm = final_drift_kms / max(sigma_mean, 1.0e-12)

    rad_sigma_kms = np.sqrt(rad_v2_sum / np.maximum(rad_counts, 1)) / 1000.0
    rad_sigma_kms = np.where(rad_counts >= 10, rad_sigma_kms, np.nan)

    summary_row = {
        "Target": cfg["name"],
        "Mode": mode,
        "N_stars": n_stars,
        "crossings": crossings,
        "sigma_LOS_mean_kms": round(sigma_mean, 4),
        "sigma_LOS_std_kms": round(sigma_std, 4),
        "observed_sigma_kms": cfg["obs_sigma_kms"],
        "observed_err_kms": cfg["obs_err_kms"],
        "r_half_final_over_init": round(final_rh_ratio, 4),
        "virial_2K_over_W_mean": round(vir_mean, 4),
        "com_drift_kms": final_drift_kms,
        "com_drift_norm": com_drift_norm,
        "runtime_s": round(runtime_s, 2),
    }

    return {
        "target_key": target_key,
        "mode": mode,
        "summary": summary_row,
        "time_series": ts_rows,
        "rad_r_norm": r_bin_centers,
        "rad_sigma_kms": rad_sigma_kms,
    }


def _pairwise_geff_accelerations(pos, nodes, softening_m):
    """
    Naive pairwise g_eff = sqrt(g_N^2 + g_N * a_Omega) solver.
    Violates Newton's 3rd law whenever interacting particles have unequal masses!
    """
    diff = pos[None, :, :] - pos[:, None, :]  # (N, N, 3): vector from i toward j
    r2 = np.sum(diff * diff, axis=-1) + (softening_m ** 2)
    r = np.sqrt(r2)
    g_ij = (nodes[None, :] * K_OMEGA) / r2
    np.fill_diagonal(g_ij, 0.0)
    g_eff_ij = np.sqrt(g_ij * g_ij + g_ij * A_OMEGA)
    return np.sum((g_eff_ij / r)[:, :, None] * diff, axis=1)


def run_momentum_conservation_proof(grid_n=64, n_proof_stars=500, n_steps=100):
    """
    Evolves an identical unequal-mass, spatially asymmetric 500-star cluster for 100 steps under:
      1. Pure K_Omega (3D FFT)
      2. Isolated QUMOND (3D FFT)
      3. Naive Pairwise g_eff
    measuring linear momentum drift |dP(t)| / M_tot [km/s] without any artificial COM recentering.
    """
    rng = np.random.default_rng(77)
    r_h = 25.0 * PC_TO_METER
    box_L = 12.0 * r_h
    solver = QUMONDMeshSolver(box_L, grid_n=grid_n, fft_workers=-1)
    eps = 0.5 * solver.dx

    pos0 = rng.normal(0.0, 0.65 * r_h, size=(n_proof_stars, 3))
    pos0[:100, 0] += 0.9 * r_h  # deliberate dipole asymmetry
    masses = rng.lognormal(0.0, 0.9, size=n_proof_stars)
    masses = masses / np.sum(masses) * (1.0e4 * M_SOLAR_KG)
    M_tot = float(np.sum(masses))
    nodes = masses / M_PROTON

    # Exact zero initial COM position & velocity
    pos0 -= np.sum(masses[:, None] * pos0, axis=0) / M_tot
    vel0 = rng.normal(0.0, 500.0, size=(n_proof_stars, 3))
    vel0 -= np.sum(masses[:, None] * vel0, axis=0) / M_tot
    v_rms_kms = float(np.sqrt(np.sum(masses * np.sum(vel0 * vel0, axis=1)) / M_tot)) / 1000.0

    # Time step scaled to QUMOND crossing time
    acc_ref = solver.get_particle_accelerations(pos0, nodes, mode="Isolated QUMOND")
    W_ref = float(np.sum(masses * np.sum(pos0 * acc_ref, axis=1)))
    dt = (2.0 * r_h / np.sqrt(abs(W_ref) / M_tot)) / 25.0
    half_dt = 0.5 * dt

    records = {"steps": np.arange(0, n_steps + 1)}
    final_stats = {}

    for label in ("Pure K_Omega (3D FFT)", "Isolated QUMOND (3D FFT)", "Naive Pairwise g_eff"):
        pos = pos0.copy()
        vel = vel0.copy()
        if "Pairwise" in label:
            acc = _pairwise_geff_accelerations(pos, nodes, eps)
        elif "Pure" in label:
            acc = solver.get_particle_accelerations(pos, nodes, mode="Pure K_Omega")
        else:
            acc = solver.get_particle_accelerations(pos, nodes, mode="Isolated QUMOND")

        drifts_kms = [float(np.linalg.norm(np.sum(masses[:, None] * vel, axis=0)) / M_tot) / 1000.0]
        for _ in range(n_steps):
            vel += half_dt * acc
            pos += dt * vel
            if "Pairwise" in label:
                acc = _pairwise_geff_accelerations(pos, nodes, eps)
            elif "Pure" in label:
                acc = solver.get_particle_accelerations(pos, nodes, mode="Pure K_Omega")
            else:
                acc = solver.get_particle_accelerations(pos, nodes, mode="Isolated QUMOND")
            vel += half_dt * acc
            p_drift_kms = float(np.linalg.norm(np.sum(masses[:, None] * vel, axis=0)) / M_tot) / 1000.0
            drifts_kms.append(p_drift_kms)

        arr_d = np.array(drifts_kms, dtype=np.float64)
        records[label] = arr_d
        final_stats[label] = {
            "drift_kms": float(arr_d[-1]),
            "drift_norm": float(arr_d[-1] / max(v_rms_kms, 1.0e-12)),
        }

    return records, final_stats, v_rms_kms


def build_benchmark_figure(results, target_keys, proof_records, proof_stats, out_png_path):
    """Renders the 4-panel dark-theme QUMOND N-body scientific benchmark figure."""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 11), facecolor="#0b1120")
    for ax in (ax1, ax2, ax3, ax4):
        ax.set_facecolor("#070c18")
        ax.tick_params(colors="#94a3b8")
        ax.grid(True, alpha=0.2, linestyle=":")

    mode_colors = {
        "Pure K_Omega": "#94a3b8",
        "Isolated QUMOND": "#38bdf8",
        "QUMOND + MW EFE": "#10b981",
    }

    # --- Panel 1 (Top-Left): Predicted vs Observed sigma_LOS across Satellites ---
    x_idx = np.arange(len(target_keys))
    width = 0.22
    for m_i, mode in enumerate(MODES):
        vals = []
        errs = []
        for tkey in target_keys:
            match = [r["summary"] for r in results if r["target_key"] == tkey and r["mode"] == mode]
            vals.append(match[0]["sigma_LOS_mean_kms"] if match else np.nan)
            errs.append(match[0]["sigma_LOS_std_kms"] if match else 0.0)
        offset = (m_i - 1.0) * width
        bars = ax1.bar(
            x_idx + offset, vals, width=width * 0.92, yerr=errs, capsize=3,
            color=mode_colors[mode], alpha=0.85, edgecolor="#f8fafc", linewidth=0.6, label=mode
        )
        for rect, val in zip(bars, vals):
            if np.isfinite(val):
                ax1.text(
                    rect.get_x() + rect.get_width() * 0.5, val * 1.08,
                    f"{val:.2f}", ha="center", va="bottom", color="#e2e8f0", fontsize=7.6
                )

    obs_vals = [SATELLITE_TARGETS[k]["obs_sigma_kms"] for k in target_keys]
    obs_errs = [SATELLITE_TARGETS[k]["obs_err_kms"] for k in target_keys]
    ax1.errorbar(
        x_idx + 1.55 * width, obs_vals, yerr=obs_errs, fmt="D", color="#facc15",
        ecolor="#facc15", elinewidth=2.0, capsize=5, markersize=7, zorder=6,
        label=r"Observed $\sigma_{\mathrm{LOS}} \pm 1\sigma$"
    )
    for xi, ov, oe in zip(x_idx, obs_vals, obs_errs):
        ax1.text(xi + 1.55 * width, (ov + oe) * 1.12, f"{ov:.2f}", ha="center", va="bottom", color="#facc15", fontsize=7.8, fontweight="bold")

    ax1.set_yscale("log")
    ax1.set_ylim(0.15, 22.0)
    ax1.set_xticks(x_idx)
    ax1.set_xticklabels([SATELLITE_TARGETS[k]["name"] for k in target_keys], fontsize=10, color="#f8fafc")
    ax1.set_ylabel(r"Velocity Dispersion $\sigma_{\mathrm{LOS}}$ [km/s]", color="#cbd5e1", fontsize=10.5)
    ax1.set_title(r"10,000-Star Equilibrium $\sigma_{\mathrm{LOS}}$ (Avg $t = 20\text{--}50\,t_{\mathrm{cross}}$)", color="#f8fafc", fontsize=11.5)
    ax1.legend(loc="upper left", facecolor="#0f172a", edgecolor="none", labelcolor="#e2e8f0", fontsize=8.0)

    # --- Panel 2 (Top-Right): Palomar 14 (or first target) Time Evolution over 50 Crossings ---
    p2_key = "palomar14" if "palomar14" in target_keys else target_keys[0]
    p2_name = SATELLITE_TARGETS[p2_key]["name"]
    ax2_r = ax2.twinx()
    ax2_r.tick_params(colors="#cbd5e1")
    for mode in MODES:
        match = [r["time_series"] for r in results if r["target_key"] == p2_key and r["mode"] == mode]
        if not match:
            continue
        ts_df = pd.DataFrame(match[0])
        col = mode_colors[mode]
        ax2.plot(ts_df["crossing_time"], ts_df["sigma_LOS_kms"], color=col, linewidth=2.0, label=rf"{mode} $\sigma_{{\mathrm{{LOS}}}}$")
        ax2_r.plot(ts_df["crossing_time"], ts_df["r_half_over_init"], color=col, linewidth=1.3, linestyle="--", alpha=0.80)

    obs_p2 = SATELLITE_TARGETS[p2_key]["obs_sigma_kms"]
    err_p2 = SATELLITE_TARGETS[p2_key]["obs_err_kms"]
    ax2.axhspan(obs_p2 - err_p2, obs_p2 + err_p2, color="#facc15", alpha=0.14, label=rf"Observed ({obs_p2:.2f}$\pm${err_p2:.2f} km/s)")
    ax2_r.axhline(1.0, color="#64748b", linestyle=":", linewidth=1.1)

    ax2.set_xlabel(r"Dynamical Crossing Time ($t / t_{\mathrm{cross}}$)", color="#cbd5e1", fontsize=10.5)
    ax2.set_ylabel(r"$\sigma_{\mathrm{LOS}}(t)$ [km/s] (Solid)", color="#38bdf8", fontsize=10.5)
    ax2_r.set_ylabel(r"Half-Mass Ratio $r_h(t)/r_h(0)$ (Dashed)", color="#cbd5e1", fontsize=10.5)
    ax2_r.set_ylim(0.75, 1.35)
    ax2.set_title(rf"{p2_name} Long-Term Virial Stability (50 Crossing Times)", color="#f8fafc", fontsize=11.5)
    ax2.legend(loc="upper right", facecolor="#0f172a", edgecolor="none", labelcolor="#e2e8f0", fontsize=7.6)

    # --- Panel 3 (Bottom-Left): Radial Projected Dispersion Profiles (20-50 t_cross Avg) ---
    radial_targets = [k for k in ("palomar14", "draco") if k in target_keys]
    if not radial_targets:
        radial_targets = target_keys[:2]
    t_styles = {radial_targets[0]: "-"}
    if len(radial_targets) > 1:
        t_styles[radial_targets[1]] = "-."

    for tk in radial_targets:
        ls = t_styles[tk]
        tname = SATELLITE_TARGETS[tk]["name"]
        for mode in MODES:
            match = [r for r in results if r["target_key"] == tk and r["mode"] == mode]
            if not match:
                continue
            r_c = match[0]["rad_r_norm"]
            sig_r = match[0]["rad_sigma_kms"]
            valid = np.isfinite(sig_r)
            ax3.plot(
                r_c[valid], sig_r[valid],
                color=mode_colors[mode], linestyle=ls, linewidth=2.1, marker="o", markersize=4,
                label=f"{tname} — {mode}"
            )

    ax3.set_xlabel(r"Projected Radius $R / r_h$", color="#cbd5e1", fontsize=10.5)
    ax3.set_ylabel(r"Radial $\sigma_{\mathrm{LOS}}(R)$ [km/s]", color="#cbd5e1", fontsize=10.5)
    ax3.set_title(r"Projected Dispersion Profiles (Averaged $20\text{--}50\,t_{\mathrm{cross}}$)", color="#f8fafc", fontsize=11.5)
    ax3.legend(loc="upper right", ncol=2 if len(radial_targets) > 1 else 1, facecolor="#0f172a", edgecolor="none", labelcolor="#e2e8f0", fontsize=7.0)

    # --- Panel 4 (Bottom-Right): Momentum Conservation Proof (3D FFT QUMOND vs Pairwise g_eff) ---
    steps = proof_records["steps"]
    proof_styles = [
        ("Naive Pairwise g_eff", "#f43f5e", "-", 2.4),
        ("Isolated QUMOND (3D FFT)", "#38bdf8", "-", 2.2),
        ("Pure K_Omega (3D FFT)", "#10b981", "--", 2.0),
    ]
    for lbl, col, ls, lw in proof_styles:
        vals = np.maximum(proof_records[lbl], 1.0e-18)
        ax4.semilogy(steps, vals, color=col, linestyle=ls, linewidth=lw, label=lbl)

    d_pair = proof_stats["Naive Pairwise g_eff"]["drift_kms"]
    d_qum = proof_stats["Isolated QUMOND (3D FFT)"]["drift_kms"]
    d_newt = proof_stats["Pure K_Omega (3D FFT)"]["drift_kms"]
    ratio_pq = d_pair / max(d_qum, 1.0e-18)
    badge = (
        f"Asymmetric 500-Star Unequal-Mass Test (100 Steps):\n"
        f"• Naive Pairwise |dP|/M = {d_pair:.3e} km/s (Self-Accelerates!)\n"
        f"• 3D FFT QUMOND  |dP|/M = {d_qum:.3e} km/s ({ratio_pq:,.0f}x lower)\n"
        f"• Pure K_Omega   |dP|/M = {d_newt:.3e} km/s (Machine Precision)"
    )
    ax4.text(
        0.97, 0.30, badge, transform=ax4.transAxes,
        ha="right", va="bottom", fontsize=7.8, color="#e2e8f0",
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#0f172a", edgecolor="#f43f5e", alpha=0.90),
    )

    ax4.set_xlabel("Integration Step (Unequal-Mass Asymmetric Cluster)", color="#cbd5e1", fontsize=10.5)
    ax4.set_ylabel(r"Linear Momentum Drift $|\Delta \mathbf{P}| / M_{\mathrm{tot}}$ [km/s]", color="#cbd5e1", fontsize=10.5)
    ax4.set_title("Linear Momentum Conservation: 3D FFT QUMOND vs. Pairwise", color="#f8fafc", fontsize=11.5)
    ax4.legend(loc="lower left", facecolor="#0f172a", edgecolor="none", labelcolor="#e2e8f0", fontsize=8.0)

    plt.tight_layout()
    fig.savefig(out_png_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="3D FFT QUMOND Particle-Mesh N-Body Testbed")
    parser.add_argument("--target", choices=["palomar14", "ngc2419", "draco", "all"], default="all")
    parser.add_argument("--particles", type=int, default=10000)
    parser.add_argument("--grid", type=int, default=64)
    parser.add_argument("--crossings", type=int, default=50)
    parser.add_argument("--mw-mass", type=float, default=8.8e10)
    parser.add_argument("--jobs", type=int, default=0, help="Worker processes (0=auto, 1=sequential)")
    parser.add_argument("--out-dir", type=str, default="qumond_audit")
    args = parser.parse_args()

    t_all0 = time.perf_counter()
    os.makedirs(args.out_dir, exist_ok=True)

    target_keys = list(SATELLITE_TARGETS.keys()) if args.target == "all" else [args.target]
    n_jobs_total = len(target_keys) * len(MODES)
    n_workers = args.jobs if args.jobs > 0 else min(n_jobs_total, max(1, (os.cpu_count() or 4) // 2), 6)
    fft_workers = -1 if n_workers == 1 else 1

    tasks = [
        (tkey, mode, args.particles, args.grid, args.crossings, args.mw_mass, fft_workers)
        for tkey in target_keys
        for mode in MODES
    ]

    print(
        f"Running 3D FFT QUMOND N-Body Testbed: targets={target_keys}, N={args.particles:,}, "
        f"grid={args.grid}^3, crossings={args.crossings}, M_MW={args.mw_mass:.2e} M_sun, workers={n_workers}"
    )

    if n_workers > 1:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            results = list(pool.map(run_single_simulation, tasks))
    else:
        results = [run_single_simulation(t) for t in tasks]

    # Run unequal-mass pairwise vs QUMOND momentum conservation proof
    proof_records, proof_stats, v_rms_proof = run_momentum_conservation_proof(grid_n=args.grid)

    # Save CSVs & validate zero NaN/Inf
    df_summary = pd.DataFrame([r["summary"] for r in results])
    df_ts = pd.DataFrame([row for r in results for row in r["time_series"]])

    for label, df_chk in (("qumond_satellite_summary.csv", df_summary), ("qumond_time_series.csv", df_ts)):
        num_arr = df_chk.select_dtypes(include=[np.number]).to_numpy(dtype=np.float64)
        assert np.all(np.isfinite(num_arr)), f"Non-finite values detected in {label}"

    summary_csv = os.path.join(args.out_dir, "qumond_satellite_summary.csv")
    ts_csv = os.path.join(args.out_dir, "qumond_time_series.csv")
    fig_png = os.path.join(args.out_dir, "qumond_nbody_benchmark.png")

    df_summary.to_csv(summary_csv, index=False)
    df_ts.to_csv(ts_csv, index=False)
    build_benchmark_figure(results, target_keys, proof_records, proof_stats, fig_png)

    wall_s = time.perf_counter() - t_all0

    print("\n## 1. Satellite Kinematics Summary (20-50 t_cross Average)\n")
    print(df_summary.to_string(index=False))

    print(f"\n## 2. Momentum Conservation Proof (N=500 Unequal-Mass Asymmetric Cluster, v_rms={v_rms_proof:.3f} km/s, 100 steps)\n")
    for lbl in ("Pure K_Omega (3D FFT)", "Isolated QUMOND (3D FFT)", "Naive Pairwise g_eff"):
        st = proof_stats[lbl]
        print(f"  - {lbl:26s}: |dP|/M = {st['drift_kms']:.4e} km/s  (|dP|/(M*v_rms) = {st['drift_norm']:.4e})")

    print(f"\nArtifacts written to '{args.out_dir}/' in {wall_s:.2f} s.")


if __name__ == "__main__":
    main()
