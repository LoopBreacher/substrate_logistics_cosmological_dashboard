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

# ==========================================
# 1. REDSHIFT PATH INTEGRATOR ENGINE
# ==========================================
def compute_redshift_path(
    ra_deg, dec_deg,
    tree, gal_positions, nodes_per_galaxy,
    max_dist_mpc=25.0, sigma_mpc=1.8, h_global=H_GLOBAL, n_steps=250
):
    """
    Computes line-of-sight velocity accumulation cz(d) and magnitude residual delta_mu(d)
    along a specific celestial direction vector.
    """
    coord = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
    cart_dir = coord.cartesian.get_xyz().value
    cart_dir = cart_dir / np.linalg.norm(cart_dir)
    
    s_array = np.linspace(0.001, max_dist_mpc, n_steps)
    h_local_samples = np.array([
        get_hlocal(tree, gal_positions, nodes_per_galaxy, s * cart_dir, sigma_mpc=sigma_mpc, h_global=h_global)
        for s in s_array
    ])
    
    cz_accumulated = np.array([simpson(h_local_samples[:i], x=s_array[:i]) for i in range(2, n_steps + 1)])
    s_eval = s_array[1:]
    
    d_naive_global = cz_accumulated / h_global
    mu_true = 5.0 * np.log10(s_eval * 1e5)
    mu_naive = 5.0 * np.log10(d_naive_global * 1e5)
    delta_mu = mu_true - mu_naive
    
    return s_eval, cz_accumulated, h_local_samples[1:], delta_mu

# ==========================================
# 2. MAIN UI EXECUTION & PLOTTING ENGINE
# ==========================================
def run_redshift_path_analysis(
    fil_ra, fil_dec, void_ra, void_dec,
    tree, gal_positions, nodes_per_galaxy,
    max_dist_mpc=25.0, sigma_mpc=1.8, h_global=H_GLOBAL
):
    """
    Main UI entry point for Tab 3 (Redshift Path Integral). Computes dual corridor paths,
    renders 2-panel diagnostic plots, and returns (fig, df_cp, metrics).
    """
    s_fil, cz_fil, h_fil, dmu_fil = compute_redshift_path(
        fil_ra, fil_dec, tree, gal_positions, nodes_per_galaxy,
        max_dist_mpc=max_dist_mpc, sigma_mpc=sigma_mpc, h_global=h_global
    )
    s_void, cz_void, h_void, dmu_void = compute_redshift_path(
        void_ra, void_dec, tree, gal_positions, nodes_per_galaxy,
        max_dist_mpc=max_dist_mpc, sigma_mpc=sigma_mpc, h_global=h_global
    )

    cz_global = s_fil * h_global
    cz_shoes = s_fil * 73.04

    # Render Diagnostic Figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), facecolor='#0b1120')
    ax1.set_facecolor('#070c18')
    ax2.set_facecolor('#070c18')

    # Panel 1: Redshift Accumulation
    ax1.plot(s_fil, cz_fil, color='#ef4444', linewidth=2.5, label='Filament Corridor (Centaurus/Local Sheet)')
    ax1.plot(s_void, cz_void, color='#38bdf8', linewidth=2.5, linestyle='--', label='Void Corridor (Local Void)')
    ax1.plot(s_fil, cz_global, color='#94a3b8', linestyle=':', linewidth=1.8, label=rf'Planck Baseline ($H_0={h_global:.2f}$)')
    ax1.plot(s_fil, cz_shoes, color='#10b981', linestyle='-.', linewidth=1.8, label=r'SH0ES Baseline ($H_0=73.04$)')
    ax1.set_xlabel('True Physical Distance $d$ [Mpc]', color='#cbd5e1', fontsize=11)
    ax1.set_ylabel(r'Accumulated Velocity $cz = \int H_{\mathrm{local}} ds$ [km/s]', color='#cbd5e1', fontsize=11)
    ax1.set_title('Redshift Accumulation Path Integral', color='#f8fafc', fontsize=12)
    ax1.tick_params(colors='#94a3b8')
    ax1.grid(True, alpha=0.2, linestyle=':')
    ax1.legend(loc='upper left', facecolor='#0f172a', edgecolor='none')

    # Panel 2: Distance Modulus Residual
    ax2.plot(s_fil, dmu_fil, color='#ef4444', linewidth=2.5, label='Filament Direction Residual')
    ax2.plot(s_void, dmu_void, color='#38bdf8', linewidth=2.5, linestyle='--', label='Void Direction Residual')
    ax2.axhline(0, color='#94a3b8', linestyle=':', label=r'Homogeneous Baseline ($\Delta\mu = 0$)')
    ax2.set_xlabel('True Physical Distance $d$ [Mpc]', color='#cbd5e1', fontsize=11)
    ax2.set_ylabel(r'SNe Ia Distance Modulus Residual $\Delta\mu$ [mag]', color='#cbd5e1', fontsize=11)
    ax2.set_title('Euclid / Rubin (LSST) Predicted Residual Anomaly', color='#f8fafc', fontsize=12)
    ax2.tick_params(colors='#94a3b8')
    ax2.grid(True, alpha=0.2, linestyle=':')
    ax2.legend(loc='lower left', facecolor='#0f172a', edgecolor='none')

    plt.tight_layout()

    # Diagnostic Summary Table with Exact Distance Interpolation
    eval_checkpoints = [5.0, 10.0, 15.0, 20.0]
    eval_checkpoints = [cp for cp in eval_checkpoints if cp <= max_dist_mpc]
    cp_rows = []
    
    for cp in eval_checkpoints:
        # Interpolate exact values at checkpoint distance cp
        cz_f_cp = float(np.interp(cp, s_fil, cz_fil))
        cz_v_cp = float(np.interp(cp, s_void, cz_void))
        dmu_f_cp = float(np.interp(cp, s_fil, dmu_fil))
        dmu_v_cp = float(np.interp(cp, s_void, dmu_void))
        
        cp_rows.append({
            "Distance [Mpc]": round(cp, 1),
            "Filament cz [km/s]": round(cz_f_cp, 1),
            "Filament H_eff [km/s/Mpc]": round(cz_f_cp / cp, 2),
            "Filament Δμ [mag]": round(dmu_f_cp, 4),
            "Void cz [km/s]": round(cz_v_cp, 1),
            "Void H_eff [km/s/Mpc]": round(cz_v_cp / cp, 2),
            "Void Δμ [mag]": round(dmu_v_cp, 4),
            "Velocity Contrast [km/s]": round(cz_f_cp - cz_v_cp, 1)
        })
        
    df_cp = pd.DataFrame(cp_rows)

    # Executive Metrics Card Data
    h_eff_fil = cz_fil / s_fil
    h_eff_void = cz_void / s_void
    metrics = {
        "Filament Peak H_eff": f"{h_eff_fil.max():.2f} km/s/Mpc",
        "Void Floor H_eff": f"{h_eff_void.min():.2f} km/s/Mpc",
        "Max Filament Δμ": f"{dmu_fil.min():.4f} mag",
        "Max Path Velocity Δcz": f"{(cz_fil - cz_void).max():+.1f} km/s"
    }

    return fig, df_cp, metrics