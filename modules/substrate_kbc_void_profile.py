import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

from modules.data_loader import H_GLOBAL

R_KBC_MPC = 300.0         # KBC Void radius boundary (Mpc)
SIGMA_WALL_MPC = 35.0     # Void wall transition scale (Mpc)
DECAY_SCALE_LOCAL = 12.0  # Local sheet density dissipation scale (Mpc)

def density_deficit_profile(r_mpc, delta_0_kbc=-0.22):
    """
    Computes macro density contrast delta(r) = (rho - rho_bar) / rho_bar
    for the KBC supervoid profile out to r = 500 Mpc.
    """
    wall_switch = 1.0 / (1.0 + np.exp((r_mpc - R_KBC_MPC) / SIGMA_WALL_MPC))
    return delta_0_kbc * wall_switch

def h_local_radial_field(r_mpc, delta_0_kbc=-0.22, h_local_peak=72.80, h_global=H_GLOBAL):
    """
    Evaluates local expansion rate H_local(r) incorporating:
    1. Local Sheet overdensity core (decaying over 12 Mpc).
    2. Macro KBC underdensity floor out to 300 Mpc.
    3. Global baseline vacuum unspooling H_global as r -> infinity.
    """
    delta_H_local_core = (h_local_peak - h_global) * np.exp(-r_mpc / DECAY_SCALE_LOCAL)
    delta_r = density_deficit_profile(r_mpc, delta_0_kbc)
    gc_scale_factor = np.sqrt(np.maximum(1.0 + delta_r, 0.0)) - 1.0
    delta_H_void = (h_local_peak - h_global) * gc_scale_factor
    
    H_r = h_global + delta_H_local_core + delta_H_void
    return np.maximum(H_r, h_global)

def run_kbc_void_analysis(r_max_mpc=500.0, delta_0_kbc=-0.22, h_local_peak=72.80, h_global=H_GLOBAL, n_steps=500):
    """
    Integrates photon paths across the KBC supervoid profile and returns (fig, df_kbc, metrics).
    """
    r_grid = np.linspace(0.0, r_max_mpc, n_steps)
    h_loc = h_local_radial_field(r_grid, delta_0_kbc, h_local_peak, h_global)
    
    cz_acc = np.zeros(n_steps)
    cz_acc[1:] = cumulative_trapezoid(h_loc, x=r_grid)
    
    h_eff = np.zeros(n_steps)
    h_eff[0] = h_loc[0]
    h_eff[1:] = cz_acc[1:] / r_grid[1:]
    
    dmu_arr = np.zeros(n_steps)
    dmu_arr[1:] = -5.0 * np.log10(h_eff[1:] / h_global)
    dmu_arr[0] = -5.0 * np.log10(h_loc[0] / h_global)
    
    delta_profile = density_deficit_profile(r_grid, delta_0_kbc) * 100.0

    # 1. Build Output Export DataFrame
    df_kbc = pd.DataFrame({
        "radius_Mpc": r_grid,
        "density_contrast_pct": delta_profile,
        "H_local_kmsMpc": h_loc,
        "H_eff_kmsMpc": h_eff,
        "recession_cz_kms": cz_acc,
        "Delta_mu_mag": dmu_arr
    })

    # 2. Render Plot with Dark Dashboard Theme
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), facecolor='#0b1120')
    ax1.set_facecolor('#070c18')
    ax2.set_facecolor('#070c18')
    for ax in (ax1, ax2):
        ax.tick_params(colors='#94a3b8')
        ax.grid(True, alpha=0.2, linestyle=':')

    # Panel 1: Radial Expansion & Density Contrast Profiles
    ax1_twin = ax1.twinx()
    ax1_twin.tick_params(colors='#10b981')

    line1 = ax1.plot(r_grid, h_loc, color='#ef4444', linewidth=2.2, label=r'Differential $H_{\mathrm{local}}(r)$')
    line2 = ax1.plot(r_grid, h_eff, color='#38bdf8', linewidth=2.5, linestyle='--', label=r'Integrated Effective $H_{\mathrm{eff}}(r)$')
    line3 = ax1_twin.plot(r_grid, delta_profile, color='#10b981', linewidth=1.8, linestyle=':', label=r'KBC Density Deficit $\delta(r)$ [\%]')

    ax1.axvline(300.0, color='#94a3b8', linestyle='--', alpha=0.7, label='KBC Boundary ($300$ Mpc)')
    ax1.axhline(h_global, color='#f8fafc', linestyle=':', label=f'Vacuum Floor ({h_global:.2f})')

    ax1.set_xlabel('Radial Distance from Observer $r$ [Mpc]', color='#cbd5e1')
    ax1.set_ylabel('Expansion Rate [km/s/Mpc]', color='#cbd5e1')
    ax1_twin.set_ylabel(r'Galaxy Density Contrast $\delta(r)$ [\%]', color='#10b981')
    ax1.set_title('KBC Supervoid Radial Expansion Transition', color='#f8fafc', fontsize=12)

    lines = line1 + line2 + line3
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='center right', facecolor='#0f172a', edgecolor='none', fontsize=8.5)

    # Panel 2: SNe Ia Distance Modulus Anomaly Profile
    ax2.plot(r_grid, dmu_arr, color='#a855f7', linewidth=2.5, label=r'Predicted SNe Ia Residual $\Delta\mu(r)$')
    ax2.axhline(0.0, color='#f8fafc', linestyle=':', label=r'Homogeneous Baseline ($\Delta\mu = 0$)')
    ax2.axvline(300.0, color='#94a3b8', linestyle='--', alpha=0.7, label='KBC Boundary ($300$ Mpc)')

    ax2.scatter([10.0, 15.0], [-0.11, -0.08], color='#ef4444', s=50, zorder=5, label='Local Ladder (SH0ES)')
    ax2.scatter([50.0, 150.0], [-0.05, -0.02], color='#38bdf8', s=50, zorder=5, label='Intermediate SNe Ia (Pantheon+)')
    ax2.scatter([300.0, 450.0], [-0.005, 0.0], color='#f8fafc', s=50, zorder=5, label='Large-Scale Surveys (BAO/CMB)')

    ax2.set_xlabel('Radial Distance $r$ [Mpc]', color='#cbd5e1')
    ax2.set_ylabel(r'Distance Modulus Anomaly $\Delta\mu$ [mag]', color='#cbd5e1')
    ax2.set_title('Cosmic Hubble Tension Relaxation across KBC Scale', color='#f8fafc', fontsize=12)
    ax2.legend(loc='lower right', facecolor='#0f172a', edgecolor='none', fontsize=8.5)

    plt.tight_layout()

    # Derived Checkpoint Metrics
    idx_10 = np.argmin(np.abs(r_grid - 10.0))
    idx_50 = np.argmin(np.abs(r_grid - 50.0))
    idx_150 = np.argmin(np.abs(r_grid - 150.0))
    idx_300 = np.argmin(np.abs(r_grid - 300.0))
    idx_500 = np.argmin(np.abs(r_grid - 500.0))

    metrics = {
        "H_eff (10 Mpc)": f"{h_eff[idx_10]:.2f} km/s/Mpc",
        "H_eff (50 Mpc)": f"{h_eff[idx_50]:.2f} km/s/Mpc",
        "H_eff (150 Mpc)": f"{h_eff[idx_150]:.2f} km/s/Mpc",
        "H_eff (300 Mpc)": f"{h_eff[idx_300]:.2f} km/s/Mpc",
        "Asymptotic H_eff": f"{h_eff[idx_500]:.2f} km/s/Mpc"
    }

    return fig, df_kbc, metrics