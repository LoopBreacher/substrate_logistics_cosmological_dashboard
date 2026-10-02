"""
Line-of-Sight Kinematics & Peculiar Velocity Solver
--------------------------------------------------
Substrate Logistics Cosmological Dashboard - Domain 1
Module: line_of_sight_kinematics.py

Physical Mechanics:
- Takes gold-standard physical benchmark distances (d_TRGB or d_Cepheid) for local targets.
- Integrates local unspooling rate H_local(s) along photon sightline rays.
- Decouples true cosmological expansion cz_exp from peculiar velocity flow:
    v_pec = cz_obs - cz_exp
- Demonstrates why naive redshift-distance inversion fails in local volume regimes.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from astropy.coordinates import SkyCoord
import astropy.units as u

from .data_loader import (
    get_hlocal,
    H_GLOBAL
)
from .substrate_bulk_flow_field import compute_line_of_sight_peculiar_velocity

def solve_line_of_sight_kinematics(
    ra_deg, dec_deg, cz_obs,
    tree, gal_positions, nodes_per_galaxy,
    d_known=None, sigma_mpc=1.8, h_global=H_GLOBAL, d_max_search=35.0, n_steps=300
):
    """
    Evaluates line-of-sight ray integration and kinematic decomposition.

    Parameters:
    -----------
    ra_deg, dec_deg : float
        Target sky coordinates in Right Ascension and Declination (degrees).
    cz_obs : float
        Observed radial recessional velocity (km/s).
    tree : scipy.spatial.cKDTree
        Spatial catalog index.
    gal_positions : np.ndarray
        3D catalog galaxy positions (Mpc).
    nodes_per_galaxy : np.ndarray
        Baryonic node count per galaxy.
    d_known : float, optional
        Gold-standard physical benchmark distance (TRGB / Cepheid) in Mpc.
    sigma_mpc : float
        Gaussian kernel smoothing scale (Mpc).
    h_global : float
        Global vacuum expansion floor (km/s/Mpc).
    d_max_search : float
        Maximum ray tracing integration depth (Mpc).
    n_steps : int
        Number of radial evaluation steps along the line of sight.

    Returns:
    --------
    dict : Computed kinematic path profile metrics and arrays.
    """
    coord = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
    cart_dir = coord.cartesian.get_xyz().value
    cart_dir = cart_dir / np.linalg.norm(cart_dir)
    
    # 1. Integrate H_local along ray trajectory out to d_max_search
    s_array = np.linspace(0.001, d_max_search, n_steps)
    h_samples = np.array([
        get_hlocal(tree, gal_positions, nodes_per_galaxy, s * cart_dir, sigma_mpc=sigma_mpc, h_global=h_global)
        for s in s_array
    ])
    
    cz_accumulated = np.array([simpson(h_samples[:i], x=s_array[:i]) for i in range(2, n_steps + 1)])
    s_eval = s_array[1:]
    
    # 2. Forward Kinematics (Benchmark Distance Evaluation)
    if d_known is not None and d_known > 0:
        cz_exp_at_d = float(np.interp(d_known, s_eval, cz_accumulated))
        v_pec_derived = cz_obs - cz_exp_at_d
    else:
        cz_exp_at_d = None
        v_pec_derived = None

    # 3. Model-Coupled Inversion
    d_approx = cz_obs / h_global
    target_pos = d_approx * cart_dir
    v_pec_model = compute_line_of_sight_peculiar_velocity(
        target_pos, tree, gal_positions, nodes_per_galaxy, sigma_mpc=sigma_mpc, h_global=h_global
    )
    
    v_expansion_target = max(cz_obs - v_pec_model, 10.0)
    
    if v_expansion_target > cz_accumulated[-1]:
        d_kinematic_inverted = v_expansion_target / h_global
    else:
        d_kinematic_inverted = float(np.interp(v_expansion_target, cz_accumulated, s_eval))
        
    d_linear_hubble = cz_obs / h_global
    
    return {
        "d_kinematic_inverted": d_kinematic_inverted,
        "d_known": d_known,
        "cz_obs": cz_obs,
        "cz_exp_at_d": cz_exp_at_d,
        "v_pec_derived": v_pec_derived,
        "v_pec_model": v_pec_model,
        "v_expansion_target": v_expansion_target,
        "d_linear_hubble": d_linear_hubble,
        "s_array": s_array,
        "h_samples": h_samples,
        "s_eval": s_eval,
        "cz_accumulated": cz_accumulated
    }


def run_line_of_sight_kinematics_analysis(
    ra_deg, dec_deg, cz_obs,
    tree, gal_positions, nodes_per_galaxy,
    d_known=None, sigma_mpc=1.8, h_global=H_GLOBAL, d_max_search=35.0
):
    """
    Main Streamlit UI entry point for Line-of-Sight Kinematics.

    Returns:
    --------
    fig : matplotlib.figure.Figure
    df_bm : pandas.DataFrame
    metrics_dict : dict
    """
    res = solve_line_of_sight_kinematics(
        ra_deg, dec_deg, cz_obs,
        tree, gal_positions, nodes_per_galaxy,
        d_known=d_known, sigma_mpc=sigma_mpc, h_global=h_global, d_max_search=d_max_search
    )
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), facecolor='#0b1120')
    ax1.set_facecolor('#070c18')
    ax2.set_facecolor('#070c18')

    # Panel 1: Ray Local Expansion Profile H_local(s)
    ax1.plot(res["s_array"], res["h_samples"], color='#38bdf8', linewidth=2.5, label=r'Local Expansion $H_{\mathrm{local}}(s)$')
    ax1.axhline(h_global, color='#94a3b8', linestyle=':', label=f'Global Vacuum Floor ({h_global:.2f})')
    
    if res["d_known"] is not None:
        ax1.axvline(res["d_known"], color='#10b981', linestyle='-', linewidth=2.0, label=f'TRGB/Cepheid Benchmark ($d={res["d_known"]:.2f}$ Mpc)')

    ax1.axvline(res["d_kinematic_inverted"], color='#f43f5e', linestyle='--', linewidth=1.8, label=f'Kinematic Inverted Dist ($d={res["d_kinematic_inverted"]:.2f}$ Mpc)')

    ax1.set_xlabel('Ray Path Distance $s$ [Mpc]', color='#cbd5e1')
    ax1.set_ylabel(r'$H_{\mathrm{local}}$ [km/s/Mpc]', color='#cbd5e1')
    ax1.set_title('Line-of-Sight Expansion Profile', color='#f8fafc', fontsize=12)
    ax1.tick_params(colors='#94a3b8')
    ax1.grid(True, alpha=0.2, linestyle=':')
    ax1.legend(loc='upper right', facecolor='#0f172a', edgecolor='none')

    # Panel 2: Accumulated Expansion & Kinematic Decomposition
    ax2.plot(res["s_eval"], res["cz_accumulated"], color='#a855f7', linewidth=2.5, label=r'Substrate Integrated Expansion $cz_{\mathrm{exp}}$')
    ax2.plot(res["s_eval"], res["s_eval"] * h_global, color='#94a3b8', linestyle=':', label=r'Linear Baseline ($cz = H_{\mathrm{global}} \cdot d$)')
    ax2.axhline(cz_obs, color='#f59e0b', linestyle=':', label=f'Observed cz ({cz_obs:.0f} km/s)')
    
    if res["d_known"] is not None and res["cz_exp_at_d"] is not None:
        ax2.scatter([res["d_known"]], [res["cz_exp_at_d"]], color='#10b981', s=100, zorder=6, label=f'Expansion cz_exp ({res["cz_exp_at_d"]:.1f} km/s)')
        ax2.annotate(f'Derived v_pec = {res["v_pec_derived"]:+.1f} km/s', 
                     xy=(res["d_known"], res["cz_exp_at_d"]),
                     xytext=(res["d_known"] + 1.2, res["cz_exp_at_d"] - 80),
                     arrowprops=dict(facecolor='#10b981', shrink=0.05, width=1, headwidth=6),
                     color='#10b981', fontweight='bold', fontsize=9)

    ax2.set_xlabel('Physical Distance $d$ [Mpc]', color='#cbd5e1')
    ax2.set_ylabel('Recession Velocity [km/s]', color='#cbd5e1')
    ax2.set_title('Line-of-Sight Kinematic Decomposition', color='#f8fafc', fontsize=12)
    ax2.tick_params(colors='#94a3b8')
    ax2.grid(True, alpha=0.2, linestyle=':')
    ax2.legend(loc='lower right', facecolor='#0f172a', edgecolor='none')

    plt.tight_layout()

    # Benchmark targets table
    benchmarks = [
        {"Target": "NGC 6946 (Fireworks)", "RA": 308.71, "Dec": 60.15, "cz": 300.0, "d_benchmark": 7.72, "Env": "Local Group Wall"},
        {"Target": "NGC 1365 (Fornax)", "RA": 53.40, "Dec": -36.14, "cz": 1636.0, "d_benchmark": 17.17, "Env": "Overdense Cluster"},
        {"Target": "M87 (Virgo Core)", "RA": 187.70, "Dec": 12.39, "cz": 1284.0, "d_benchmark": 16.50, "Env": "Dense Cluster Core"},
        {"Target": "M83 (Southern Pinwheel)", "RA": 204.25, "Dec": -29.87, "cz": 513.0, "d_benchmark": 4.90, "Env": "Filament Core"},
        {"Target": "KK246 (Deep Void)", "RA": 300.00, "Dec": -28.50, "cz": 416.0, "d_benchmark": 7.83, "Env": "Deep Local Void"},
    ]
    
    bm_rows = []
    
    user_row = {
        "Galaxy Target": f"Active Query (RA: {ra_deg:.1f}°, Dec: {dec_deg:.1f}°)",
        "Environment": "User Target",
        "Observed cz [km/s]": round(cz_obs, 1),
        "TRGB/Cepheid Benchmark [Mpc]": round(d_known, 2) if d_known else "N/A",
        "Substrate Expansion cz_exp [km/s]": round(res["cz_exp_at_d"], 1) if res["cz_exp_at_d"] else round(res["v_expansion_target"], 1),
        "Derived v_pec [km/s]": round(res["v_pec_derived"], 1) if res["v_pec_derived"] is not None else round(res["v_pec_model"], 1),
        "Linear Hubble Distance [Mpc]": round(res["d_linear_hubble"], 2)
    }
    bm_rows.append(user_row)

    for b in benchmarks:
        r_b = solve_line_of_sight_kinematics(
            b['RA'], b['Dec'], b['cz'],
            tree, gal_positions, nodes_per_galaxy,
            d_known=b['d_benchmark'], sigma_mpc=sigma_mpc, h_global=h_global, d_max_search=d_max_search
        )
        bm_rows.append({
            "Galaxy Target": b['Target'],
            "Environment": b['Env'],
            "Observed cz [km/s]": round(b['cz'], 1),
            "TRGB/Cepheid Benchmark [Mpc]": round(b['d_benchmark'], 2),
            "Substrate Expansion cz_exp [km/s]": round(r_b["cz_exp_at_d"], 1),
            "Derived v_pec [km/s]": round(r_b["v_pec_derived"], 1),
            "Linear Hubble Distance [Mpc]": round(r_b["d_linear_hubble"], 2)
        })
        
    df_bm = pd.DataFrame(bm_rows)

    metrics = {
        "Substrate Expansion (cz_exp)": f"{res['cz_exp_at_d']:.1f} km/s" if res['cz_exp_at_d'] else f"{res['v_expansion_target']:.1f} km/s",
        "Derived Peculiar Velocity (v_pec)": f"{res['v_pec_derived']:+.1f} km/s" if res['v_pec_derived'] is not None else f"{res['v_pec_model']:+.1f} km/s",
        "TRGB/Cepheid Benchmark": f"{d_known:.2f} Mpc" if d_known else "N/A",
        "Linear Hubble Distance": f"{res['d_linear_hubble']:.2f} Mpc"
    }

    return fig, df_bm, metrics