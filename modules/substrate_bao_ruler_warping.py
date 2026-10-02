import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

from .data_loader import (
    get_hlocal,
    R_BAO_MPC,
    H_GLOBAL
)

# ==========================================
# 1. ALCOCK-PACZYNSKI BAO WARPING SOLVER
# ==========================================
def calculate_bao_warping(
    ra_deg, dec_deg,
    tree, gal_positions, nodes_per_galaxy,
    path_depth_mpc=R_BAO_MPC, n_steps=120, sigma_mpc=1.8, h_global=H_GLOBAL
):
    """
    Computes path-integrated expansion H_eff and Alcock-Paczynski 
    distortion factors integrated across the physical BAO sound horizon scale (default 149.21 Mpc).
    """
    coord = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
    cart_dir = coord.cartesian.get_xyz().value
    cart_dir = cart_dir / np.linalg.norm(cart_dir)

    s_steps = np.linspace(0.1, path_depth_mpc, n_steps)
    h_samples = [
        get_hlocal(tree, gal_positions, nodes_per_galaxy, s * cart_dir, sigma_mpc=sigma_mpc, h_global=h_global)
        for s in s_steps
    ]
    h_eff = float(np.mean(h_samples))

    C_factor = h_eff / h_global
    alpha_parallel = 1.0 / C_factor
    alpha_perp = np.sqrt(C_factor)
    F_AP = alpha_perp / alpha_parallel

    r_parallel_mpc = R_BAO_MPC * alpha_parallel
    r_perp_mpc = R_BAO_MPC * alpha_perp

    return h_eff, C_factor, alpha_parallel, alpha_perp, F_AP, r_parallel_mpc, r_perp_mpc

# ==========================================
# 2. MAIN UI EXECUTION & PLOTTING ENGINE
# ==========================================
def run_bao_warping_analysis(
    tree, gal_positions, nodes_per_galaxy,
    path_depth_mpc=R_BAO_MPC, sigma_mpc=1.8, h_global=H_GLOBAL, dec_scan_deg=0.0
):
    """
    Main UI entry point for Tab BAO Ruler Warping. 
    Renders BAO deformed standard ruler shells, AP scan angle, and returns (fig, df_results, metrics).
    """
    target_environments = [
        {"name": "Virgo Cluster Core",      "ra": 187.7, "dec": 12.4,  "env": "Overdense Cluster Core"},
        {"name": "Centaurus Filament Peak", "ra": 201.0, "dec": -29.8, "env": "Supercluster Filament Axis"},
        {"name": "Perseus-Pisces Ridge",   "ra": 49.9,  "dec": 41.5,  "env": "Major Overdense Wall"},
        {"name": "Coma Cluster Node",       "ra": 194.9, "dec": 27.9,  "env": "Rich Massive Cluster"},
        {"name": "Eridanus Void Core",      "ra": 50.0,  "dec": -20.0, "env": "Underdense Void Wall"},
        {"name": "Southern Local Void",    "ra": 280.0, "dec": -20.0, "env": "Deep Cosmic Void Floor"},
    ]

    results = []
    for target in target_environments:
        h_obs, C_fac, a_par, a_perp, F_AP, r_par, r_perp = calculate_bao_warping(
            target["ra"], target["dec"],
            tree, gal_positions, nodes_per_galaxy,
            path_depth_mpc=path_depth_mpc, sigma_mpc=sigma_mpc, h_global=h_global
        )
        results.append({
            "Environment": target["name"],
            "Target Sky Region": target["env"],
            "Effective H_eff (km/s/Mpc)": round(h_obs, 2),
            "Correction Factor C": round(C_fac, 4),
            "Alpha Parallel (Radial)": round(a_par, 4),
            "Alpha Perp (Transverse)": round(a_perp, 4),
            "Alcock-Paczynski F_AP": round(F_AP, 4),
            "Radial Ruler r_par (Mpc)": round(r_par, 2),
            "Transverse Ruler r_perp (Mpc)": round(r_perp, 2),
            "Dilation Anomaly [%]": round((F_AP - 1.0) * 100.0, 2)
        })
    df_results = pd.DataFrame(results)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), facecolor='#0b1120')
    ax1.set_facecolor('#070c18')
    ax2.set_facecolor('#070c18')

    theta_grid = np.linspace(0, 2 * np.pi, 200)

    # Panel 1: Real-Space Deformed BAO Standard Ruler Shells
    x_ideal = R_BAO_MPC * np.cos(theta_grid)
    y_ideal = R_BAO_MPC * np.sin(theta_grid)
    ax1.plot(x_ideal, y_ideal, color='#94a3b8', linestyle=':', linewidth=2.0, label=f'Isotropic Baseline ($r_s = {R_BAO_MPC:.1f}$ Mpc)')

    # Representative Filament Shell
    _, _, _, _, f_ap_fil, r_par_fil, r_perp_fil = calculate_bao_warping(
        201.0, -29.8, tree, gal_positions, nodes_per_galaxy,
        path_depth_mpc=path_depth_mpc, sigma_mpc=sigma_mpc, h_global=h_global
    )
    x_fil = r_perp_fil * np.cos(theta_grid)
    y_fil = r_par_fil * np.sin(theta_grid)
    ax1.plot(x_fil, y_fil, color='#ef4444', linewidth=2.5, label=f'Filament Warped ($F_{{\\mathrm{{AP}}}} = {f_ap_fil:.3f}$)')

    # Representative Void Floor Shell
    _, _, _, _, f_ap_void, r_par_void, r_perp_void = calculate_bao_warping(
        280.0, -20.0, tree, gal_positions, nodes_per_galaxy,
        path_depth_mpc=path_depth_mpc, sigma_mpc=sigma_mpc, h_global=h_global
    )
    x_void = r_perp_void * np.cos(theta_grid)
    y_void = r_par_void * np.sin(theta_grid)
    ax1.plot(x_void, y_void, color='#38bdf8', linewidth=2.2, linestyle='--', label=f'Void Floor ($F_{{\\mathrm{{AP}}}} = {f_ap_void:.3f}$)')

    ax1.set_aspect('equal')
    ax1.set_xlabel(r'Transverse Extent $r_{\perp}$ [Mpc]', color='#cbd5e1', fontsize=11)
    ax1.set_ylabel(r'Radial Line-of-Sight Extent $r_{\parallel}$ [Mpc]', color='#cbd5e1', fontsize=11)
    ax1.set_title(f'Real-Space BAO Ruler Deformation (Ray Depth: {path_depth_mpc:.1f} Mpc)', color='#f8fafc', fontsize=12)
    ax1.tick_params(colors='#94a3b8')
    ax1.grid(True, alpha=0.2, linestyle=':')
    ax1.legend(loc='lower right', facecolor='#0f172a', edgecolor='none', labelcolor='#e2e8f0')

    # Panel 2: Alcock-Paczynski Distortion Parameter Scan
    ra_scan = np.linspace(0, 360, 60)
    f_ap_scan = [
        calculate_bao_warping(
            ra, dec_scan_deg, tree, gal_positions, nodes_per_galaxy,
            path_depth_mpc=path_depth_mpc, sigma_mpc=sigma_mpc, h_global=h_global
        )[4]
        for ra in ra_scan
    ]

    ax2.plot(ra_scan, f_ap_scan, color='#a855f7', linewidth=2.5, label=f'Alcock-Paczynski Factor $F_\\mathrm{{AP}}(\\theta, \\delta={dec_scan_deg:.0f}^\\circ)$')
    ax2.axhline(1.0, color='#94a3b8', linestyle=':', label=r'Isotropic Baseline ($F_{\mathrm{AP}} = 1.0$)')
    ax2.set_xlabel('Right Ascension Scan Angle [Degrees]', color='#cbd5e1', fontsize=11)
    ax2.set_ylabel(r'AP Warping Parameter $F_{\mathrm{AP}} = \alpha_{\perp} / \alpha_{\parallel}$', color='#cbd5e1', fontsize=11)
    ax2.set_title(r'DESI / Euclid Anisotropic BAO Dilation Forecast', color='#f8fafc', fontsize=12)
    ax2.tick_params(colors='#94a3b8')
    ax2.grid(True, alpha=0.2, linestyle=':')
    ax2.legend(loc='upper right', facecolor='#0f172a', edgecolor='none', labelcolor='#e2e8f0')

    plt.tight_layout()

    metrics = {
        "BAO Natural Ruler $r_s$": f"{R_BAO_MPC:.2f} Mpc",
        "Ray Integration Depth": f"{path_depth_mpc:.1f} Mpc",
        "Filament Peak $F_{\\mathrm{AP}}$": f"{f_ap_fil:.4f}",
        "Void Floor $F_{\\mathrm{AP}}$": f"{f_ap_void:.4f}"
    }

    return fig, df_results, metrics