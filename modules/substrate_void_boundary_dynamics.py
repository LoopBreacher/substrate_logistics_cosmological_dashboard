"""
Substrate Void Boundary Dynamics & Boundary Velocity Shear Module
-----------------------------------------------------------------
Substrate Logistics Cosmological Dashboard - Domain 2
Module: substrate_void_boundary_dynamics.py

Physical Mechanics:
- Evaluates continuous intergalactic background density (rho_macro) and filament density (rho_filament).
- Derives localized expansion unspooling boosts H_local(r) from 3D compiled catalog node density.
- Computes radial recession velocity v_rec(r) and boundary velocity shear spikes dv_rec/dr.
- Demonstrates vacuum floor locking (67.42 km/s/Mpc) in void cores and wall unspooling peaks.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

from .data_loader import (
    K_OMEGA,
    MPC_TO_METER,
    UNIT_CONV,
    H_GLOBAL
)

# ==========================================
# 1. DUAL-SCALE INTERGALACTIC DENSITY FIELD
# ==========================================
def get_dual_scale_density(pos_mpc, tree, gal_positions, nodes_per_galaxy, sigma_macro=1.5, sigma_filament=0.5):
    """
    Evaluates continuous intergalactic background density (rho_macro) 
    and localized filament structural density (rho_filament) along the ray trajectory.
    """
    # 1. Macro intergalactic volume density
    idx_macro = tree.query_ball_point(pos_mpc, r=3.0 * sigma_macro)
    if not idx_macro:
        rho_macro = 0.0
    else:
        dists_m = np.linalg.norm(gal_positions[idx_macro] - pos_mpc, axis=1)
        kernel_m = np.exp(-0.5 * (dists_m / sigma_macro)**2) / ((2.0 * np.pi * sigma_macro**2)**1.5)
        rho_macro = np.sum(nodes_per_galaxy[idx_macro] * kernel_m) / (MPC_TO_METER**3)

    # 2. Localized filament core density
    idx_filament = tree.query_ball_point(pos_mpc, r=3.0 * sigma_filament)
    if not idx_filament:
        rho_filament = 0.0
    else:
        dists_f = np.linalg.norm(gal_positions[idx_filament] - pos_mpc, axis=1)
        kernel_f = np.exp(-0.5 * (dists_f / sigma_filament)**2) / ((2.0 * np.pi * sigma_filament**2)**1.5)
        rho_filament = np.sum(nodes_per_galaxy[idx_filament] * kernel_f) / (MPC_TO_METER**3)

    return rho_macro, rho_filament

# ==========================================
# 2. MAIN UI EXECUTION & PLOTTING ENGINE
# ==========================================
def run_void_boundary_analysis(
    tree, gal_positions, nodes_per_galaxy,
    void_center_mpc=None, max_r_mpc=18.0, n_pts=200, h_global=H_GLOBAL
):
    """
    Main UI entry point for Void Wall Dynamics.
    Renders clean density profiles, radial unspooling curves, and boundary velocity shear.
    """
    if void_center_mpc is None:
        void_center_mpc = np.array([-2.5, -6.0, -1.5])  # Local Void Core (Mpc)

    # Ray trajectory from Local Void core toward Virgo Cluster wall
    virgo_coord = SkyCoord(ra=187.7*u.deg, dec=12.4*u.deg, distance=16.5*u.Mpc, frame='icrs')
    virgo_pos = virgo_coord.cartesian.get_xyz().value
    ray_dir = virgo_pos - void_center_mpc
    sheet_dir = ray_dir / np.linalg.norm(ray_dir)

    r_pts = np.linspace(0.0, max_r_mpc, n_pts)
    h_prof = np.zeros(n_pts)
    rho_macro_prof = np.zeros(n_pts)
    rho_filament_prof = np.zeros(n_pts)

    for i, r in enumerate(r_pts):
        p = void_center_mpc + r * sheet_dir
        r_m, r_f = get_dual_scale_density(p, tree, gal_positions, nodes_per_galaxy)

        # Unspooling boost driven directly by compiled intergalactic matter density
        dH_unspool = np.sqrt((8.0 * np.pi * K_OMEGA * r_m) / 3.0) / UNIT_CONV if r_m > 0 else 0.0

        h_prof[i] = h_global + dH_unspool
        rho_macro_prof[i] = r_m
        rho_filament_prof[i] = r_f

    v_rec = r_pts * h_prof
    v_shear = np.gradient(v_rec, r_pts)

    # Render Visual Plot
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 9), sharex=True, facecolor='#0b1120')
    for ax in (ax1, ax2, ax3):
        ax.set_facecolor('#070c18')
        ax.tick_params(colors='#94a3b8')

    # Panel 1: Dual-Scale Continuous Intergalactic Density
    ax1.plot(r_pts, rho_filament_prof, color='#ec4899', linewidth=2.0, label=r'Filament Core Density $\rho_{\mathrm{filament}}(r)$')
    ax1.plot(r_pts, rho_macro_prof, color='#06b6d4', linewidth=2.2, label=r'Macro Field Density $\rho_{\mathrm{macro}}(r)$')
    ax1.set_ylabel(r'$\rho_N$ [nodes/$\mathrm{m}^3$]', color='#cbd5e1')
    ax1.set_title('Void-to-Wall Continuous Intergalactic Density Profiles', color='#f8fafc', fontsize=12)
    ax1.grid(True, alpha=0.2, linestyle=':')
    ax1.legend(loc='upper left', facecolor='#0f172a', edgecolor='none')

    # Panel 2: Radial Local Expansion Field H_local(r)
    ax2.plot(r_pts, h_prof, color='#f59e0b', linewidth=2.5, label=r'Local Expansion Rate $H_{\mathrm{local}}(r)$')
    ax2.axhline(h_global, color='#38bdf8', linestyle=':', label=f'Global Vacuum Floor ({h_global:.2f} km/s/Mpc)')
    ax2.set_ylabel(r'$H_{\mathrm{local}}$ [km/s/Mpc]', color='#cbd5e1')
    ax2.set_title('Radial Expansion Rate Unspooling Field', color='#f8fafc', fontsize=12)
    ax2.grid(True, alpha=0.2, linestyle=':')
    ax2.legend(loc='upper left', facecolor='#0f172a', edgecolor='none')

    # Panel 3: Boundary Velocity Shear Spike dv_rec/dr
    ax3.plot(r_pts, v_shear, color='#10b981', linewidth=2.2, label=r'Boundary Velocity Shear $dv_{\mathrm{rec}}/dr$')
    ax3.axhline(h_global, color='#38bdf8', linestyle=':', label='Homogeneous Baseline')
    ax3.set_xlabel('Radial Distance from Void Center $r$ [Mpc]', color='#cbd5e1')
    ax3.set_ylabel('Shear [km/s/Mpc]', color='#cbd5e1')
    ax3.set_title('Cosmic Void Wall Velocity Shear Dynamics', color='#f8fafc', fontsize=12)
    ax3.grid(True, alpha=0.2, linestyle=':')
    ax3.legend(loc='upper left', facecolor='#0f172a', edgecolor='none')

    plt.tight_layout()

    # Construct Clean Output DataFrame
    df_table = pd.DataFrame({
        "r_mpc": np.round(r_pts, 3),
        "rho_macro_nodes_m3": rho_macro_prof,
        "rho_filament_nodes_m3": rho_filament_prof,
        "H_local_kmsMpc": np.round(h_prof, 3),
        "v_recession_kms": np.round(v_rec, 2),
        "v_shear_kmsMpc": np.round(v_shear, 3)
    })

    metrics = {
        "Void Core Expansion": f"{h_prof[0]:.2f} km/s/Mpc",
        "Peak Wall Expansion": f"{np.max(h_prof):.2f} km/s/Mpc",
        "Max Velocity Shear": f"{np.max(v_shear):.2f} km/s/Mpc",
        "Shear Peak Distance": f"{r_pts[np.argmax(v_shear)]:.1f} Mpc"
    }

    return fig, df_table, metrics