import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

from .data_loader import (
    get_hlocal_batch,
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
    distortion factors integrated across the physical BAO sound horizon scale (default 149.25 Mpc).
    """
    coord = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
    cart_dir = coord.cartesian.get_xyz().value
    cart_dir = cart_dir / np.linalg.norm(cart_dir)

    s_steps = np.linspace(0.1, path_depth_mpc, n_steps)
    ray_pts = np.outer(s_steps, cart_dir)
    h_samples = get_hlocal_batch(
        tree, gal_positions, nodes_per_galaxy, ray_pts,
        sigma_mpc=sigma_mpc, h_global=h_global
    )
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

    n_steps = 120
    s_steps = np.linspace(0.1, path_depth_mpc, n_steps)

    # Batched evaluation for all 6 target sightlines
    env_coords = SkyCoord(
        ra=[t["ra"] for t in target_environments] * u.deg,
        dec=[t["dec"] for t in target_environments] * u.deg,
        frame='icrs'
    )
    env_dirs = env_coords.cartesian.get_xyz().value.T
    env_dirs = env_dirs / np.linalg.norm(env_dirs, axis=1, keepdims=True)
    env_ray_pts = (env_dirs[:, np.newaxis, :] * s_steps[np.newaxis, :, np.newaxis]).reshape(-1, 3)

    h_env_batch = get_hlocal_batch(
        tree, gal_positions, nodes_per_galaxy, env_ray_pts,
        sigma_mpc=sigma_mpc, h_global=h_global
    ).reshape(len(target_environments), n_steps)
    h_eff_env = np.mean(h_env_batch, axis=1)

    results = []
    env_lookup = {}
    for idx_t, target in enumerate(target_environments):
        h_obs = float(h_eff_env[idx_t])
        C_fac = h_obs / h_global
        a_par = 1.0 / C_fac
        a_perp = np.sqrt(C_fac)
        F_AP = a_perp / a_par
        r_par = R_BAO_MPC * a_par
        r_perp = R_BAO_MPC * a_perp

        env_lookup[target["name"]] = (F_AP, r_par, r_perp)
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

    # Reuse already-computed Representative Filament & Void Floor Shells
    f_ap_fil, r_par_fil, r_perp_fil = env_lookup["Centaurus Filament Peak"]
    x_fil = r_perp_fil * np.cos(theta_grid)
    y_fil = r_par_fil * np.sin(theta_grid)
    ax1.plot(x_fil, y_fil, color='#ef4444', linewidth=2.5, label=f'Filament Warped ($F_{{\\mathrm{{AP}}}} = {f_ap_fil:.3f}$)')

    f_ap_void, r_par_void, r_perp_void = env_lookup["Southern Local Void"]
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

    # Panel 2: Batched Alcock-Paczynski Distortion Parameter Scan (60 directions)
    ra_scan = np.linspace(0, 360, 60)
    scan_coords = SkyCoord(ra=ra_scan * u.deg, dec=np.full_like(ra_scan, dec_scan_deg) * u.deg, frame='icrs')
    scan_dirs = scan_coords.cartesian.get_xyz().value.T
    scan_dirs = scan_dirs / np.linalg.norm(scan_dirs, axis=1, keepdims=True)
    scan_ray_pts = (scan_dirs[:, np.newaxis, :] * s_steps[np.newaxis, :, np.newaxis]).reshape(-1, 3)

    h_scan_batch = get_hlocal_batch(
        tree, gal_positions, nodes_per_galaxy, scan_ray_pts,
        sigma_mpc=sigma_mpc, h_global=h_global
    ).reshape(len(ra_scan), n_steps)
    h_eff_scan = np.mean(h_scan_batch, axis=1)
    f_ap_scan = (h_eff_scan / h_global) ** 1.5

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