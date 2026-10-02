import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .data_loader import (
    C_LIGHT,
    H_GLOBAL,
    K_OMEGA,
    M_PROTON,
    M_SOLAR_KG,
    KPC_TO_METER,
    A_OMEGA,
    RAD_TO_ARCSEC,
    compute_a_omega,
    get_compiled_node_density,
    get_hlocal
)

# Preset Benchmark Lens Targets
LENS_PRESETS = {
    "Abell 1689 Core (Cluster Lens)": {
        "M_baryon_solar": 1.2e13,      # 12.0 x 10^12 M_sun
        "scale_radius_kpc": 50.0,
        "desc": "Rich cluster core showing massive arc distortion."
    },
    "Bullet Cluster Centroid (1E 0657-56)": {
        "M_baryon_solar": 0.8e13,      # 8.0 x 10^12 M_sun
        "scale_radius_kpc": 35.0,
        "desc": "Compact galaxy centroid offset from gas cloud."
    },
    "Milky Way (Total Baryonic Footprint)": {
        "M_baryon_solar": 1.5e11,      # 0.15 x 10^12 M_sun
        "scale_radius_kpc": 8.0,
        "desc": "Total MW baryonic mass."
    }
}

def _compute_zero_parameter_lensing_bounds(M_baryon_solar, custom_scale_radius_kpc=None, a_omega=A_OMEGA):
    """
    Computes zero-parameter physical crossover radius b_cross where theta_baryon == theta_scaffold,
    and 1:1 total mass horizon b_1to1_total where M_DM_inferred == M_baryon_total.
    """
    C_scaff = (a_omega * M_PROTON) / (K_OMEGA * M_SOLAR_KG)  # M_sun / m^2
    C_scaff_kpc2 = C_scaff * (KPC_TO_METER**2)               # M_sun / kpc^2
    
    b_1to1_total = np.sqrt(M_baryon_solar / C_scaff_kpc2)
    
    if custom_scale_radius_kpc is None or custom_scale_radius_kpc <= 0:
        scale_radius_kpc = np.clip(8.0 * ((M_baryon_solar / 1.5e11) ** (1/3)), 0.05, 50.0)
    else:
        scale_radius_kpc = custom_scale_radius_kpc
        
    b_cross = b_1to1_total - scale_radius_kpc
    
    if b_cross > 0.05 * b_1to1_total:
        b_focus = b_cross
    else:
        b_focus = max(b_1to1_total, scale_radius_kpc, 0.1)
        
    b_min = max(0.01, 0.05 * b_focus)
    b_max = max(10.0, 3.5 * max(b_focus, b_1to1_total))
    
    return b_min, b_max, b_focus, b_1to1_total, scale_radius_kpc

def _format_mass_label_matplotlib(mass_solar):
    """Returns formatted mass string using LaTeX mathtext for clean rendering."""
    if mass_solar >= 1e11:
        return f"{mass_solar / 1e12:.2f} × $10^{{12}} M_\\odot$"
    elif mass_solar >= 1e8:
        return f"{mass_solar / 1e9:.2f} × $10^{{9}} M_\\odot$"
    elif mass_solar >= 1e5:
        return f"{mass_solar / 1e6:.2f} × $10^{{6}} M_\\odot$"
    else:
        return f"{mass_solar:.2e} $M_\\odot$"

def calculate_lensing_profile(
    M_baryon_solar, 
    b_kpc_array, 
    a_omega_active=A_OMEGA, 
    scale_radius_kpc=None,
    lens_distance_mpc=670.0
):
    b_kpc = np.maximum(np.atleast_1d(b_kpc_array), 1e-4)
    b_meters = b_kpc * KPC_TO_METER
    M_baryon_kg = M_baryon_solar * M_SOLAR_KG
    
    # 1. Enclosed Baryonic Mass Profile M_enc(b)
    if scale_radius_kpc is not None and scale_radius_kpc > 0:
        a_meters = scale_radius_kpc * KPC_TO_METER
        M_enc_kg = M_baryon_kg * (b_meters**2) / ((b_meters + a_meters)**2)
    else:
        M_enc_kg = M_baryon_kg

    N_B_enc = M_enc_kg / M_PROTON

    # 2. Deflection Angles [radians]
    theta_baryon_rad = (4.0 * N_B_enc * K_OMEGA) / ((C_LIGHT**2) * b_meters)
    theta_scaffold_rad = (4.0 * a_omega_active * b_meters) / (C_LIGHT**2)
    theta_total_rad = theta_baryon_rad + theta_scaffold_rad
    
    # Convert to arcseconds
    theta_baryon_arcsec = theta_baryon_rad * RAD_TO_ARCSEC
    theta_scaffold_arcsec = theta_scaffold_rad * RAD_TO_ARCSEC
    theta_total_arcsec = theta_total_rad * RAD_TO_ARCSEC
    
    # 3. Inferred Dark Matter Halo Illusion
    N_DM = (a_omega_active * (b_meters**2)) / K_OMEGA
    M_DM_inferred_solar = (N_DM * M_PROTON) / M_SOLAR_KG
    
    # 4. Dimensionless Angular Convergence Profile kappa(theta)
    kpc_per_arcsec = lens_distance_mpc * (np.pi / 648000.0) * 1000.0
    theta_b_arcsec = b_kpc / kpc_per_arcsec
    
    kappa_baryon = np.gradient(theta_b_arcsec * theta_baryon_arcsec, theta_b_arcsec) / (2.0 * theta_b_arcsec)
    kappa_scaffold = (4.0 * a_omega_active * KPC_TO_METER * kpc_per_arcsec / (C_LIGHT**2)) * RAD_TO_ARCSEC
    kappa_total = kappa_baryon + kappa_scaffold
    
    return theta_total_arcsec, theta_baryon_arcsec, theta_scaffold_arcsec, M_DM_inferred_solar, kappa_total

def run_lensing_analysis(
    lens_target_name="Abell 1689 Core (Cluster Lens)",
    custom_M_baryon=None,
    custom_scale_radius_kpc=None,
    h_global=H_GLOBAL,
    use_log_scale=False,
    use_catalog_query=False,
    target_pos_mpc=None,
    sigma_mpc=1.8,
    tree=None,
    gal_positions=None,
    nodes_per_galaxy=None
):
    if use_catalog_query and target_pos_mpc is not None and tree is not None:
        rho_N_local = get_compiled_node_density(tree, gal_positions, nodes_per_galaxy, target_pos_mpc, sigma_mpc=sigma_mpc)
        h_local = get_hlocal(tree, gal_positions, nodes_per_galaxy, target_pos_mpc, sigma_mpc=sigma_mpc, h_global=h_global)
        a_omega_used = compute_a_omega(h_local)
        M_baryon = custom_M_baryon if custom_M_baryon is not None else 1.0e11
        target_title = f"{lens_target_name}\n(X={target_pos_mpc[0]:.1f}, Y={target_pos_mpc[1]:.1f}, Z={target_pos_mpc[2]:.1f} Mpc)"
    else:
        if lens_target_name in LENS_PRESETS:
            preset = LENS_PRESETS[lens_target_name]
            M_baryon = custom_M_baryon if custom_M_baryon is not None else preset["M_baryon_solar"]
            custom_scale_radius_kpc = custom_scale_radius_kpc if custom_scale_radius_kpc is not None else preset.get("scale_radius_kpc")
        else:
            M_baryon = custom_M_baryon if custom_M_baryon is not None else 1.0e12

        a_omega_used = compute_a_omega(h_global)
        h_local = h_global
        rho_N_local = 0.0
        target_title = lens_target_name

    # Compute dynamic physical crossover radius & impact bounds
    b_min, b_max, b_focus, b_1to1_total, scale_radius = _compute_zero_parameter_lensing_bounds(
        M_baryon_solar=M_baryon,
        custom_scale_radius_kpc=custom_scale_radius_kpc,
        a_omega=a_omega_used
    )

    b_array = np.geomspace(b_min, b_max, 300) if use_log_scale else np.linspace(b_min, b_max, 300)
    
    theta_tot, theta_bar, theta_scaff, M_DM_inf, kappa_tot = calculate_lensing_profile(
        M_baryon_solar=M_baryon,
        b_kpc_array=b_array,
        a_omega_active=a_omega_used,
        scale_radius_kpc=scale_radius
    )
    
    # Adaptive plot mass scaling
    if M_baryon >= 1e11:
        mass_factor, mass_unit = 1e12, "10¹² M_\\odot"
    elif M_baryon >= 1e8:
        mass_factor, mass_unit = 1e9, "10⁹ M_\\odot"
    elif M_baryon >= 1e5:
        mass_factor, mass_unit = 1e6, "10⁶ M_\\odot"
    else:
        mass_factor, mass_unit = 1.0, "M_\\odot"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), facecolor='#0b1120')
    ax1.set_facecolor('#070c18')
    ax2.set_facecolor('#070c18')
    
    for ax in (ax1, ax2):
        ax.tick_params(colors='#94a3b8')
        ax.grid(True, alpha=0.2, linestyle=':')

    # Panel 1: Deflection Angle Profile theta(b)
    plot_fn = ax1.loglog if use_log_scale else ax1.plot
    plot_fn(b_array, theta_tot, color='#ef4444', linewidth=2.5, label=r'Total Substrate Deflection $\theta_{\mathrm{total}}(b)$')
    plot_fn(b_array, theta_bar, color='#38bdf8', linewidth=1.8, linestyle='--', label=r'Baryonic Component $\theta_{\mathrm{baryon}}$')
    plot_fn(b_array, theta_scaff, color='#a855f7', linewidth=1.8, linestyle=':', label=r'Scaffolding Tax $\theta_{\mathrm{scaffold}} \propto b$')
    
    ax1.axvline(b_focus, color='#10b981', linestyle=':', label=f'Crossover Focus ($b={b_focus:.2f}$ kpc)')
    
    ax1.set_xlabel('Impact Parameter $b$ [kpc]', color='#cbd5e1', fontsize=11)
    ax1.set_ylabel(r'Optical Deflection Angle $\theta$ [arcsec]', color='#cbd5e1', fontsize=11)
    ax1.set_title(f'{target_title}\nLight Deflection Profile', color='#f8fafc', fontsize=11)
    ax1.legend(loc='upper right', facecolor='#0f172a', edgecolor='none', labelcolor='#e2e8f0')

    # Panel 2: Inferred Halo Mass Illusion M_DM(b)
    plot_fn2 = ax2.loglog if use_log_scale else ax2.plot
    plot_fn2(b_array, M_DM_inf / mass_factor, color='#f59e0b', linewidth=2.5, label=r'Inferred Halo Mass $M_{\mathrm{DM}}(b) = \frac{a_\Omega b^2}{\mathcal{K}_\Omega}$')
    
    baryon_label_str = f"True Baryonic Baseline ({_format_mass_label_matplotlib(M_baryon)})"
    ax2.axhline(M_baryon / mass_factor, color='#38bdf8', linestyle='--', label=baryon_label_str)
    
    # 1:1 Total Mass Horizon Guide Line
    ax2.axvline(b_1to1_total, color='#38bdf8', linestyle=':', label=f'1:1 Total Mass Horizon ($b={b_1to1_total:.2f}$ kpc)')
    
    ax2.set_xlabel('Impact Parameter $b$ [kpc]', color='#cbd5e1', fontsize=11)
    ax2.set_ylabel(f'Inferred Halo Mass [${mass_unit}$]', color='#cbd5e1', fontsize=11)
    ax2.set_title(r'Dark Matter Halo Illusion: $M_{\mathrm{DM}} \propto b^2$', color='#f8fafc', fontsize=11)
    ax2.legend(loc='upper left', facecolor='#0f172a', edgecolor='none', labelcolor='#e2e8f0')

    plt.tight_layout()

    # Data Table Export
    df_lensing = pd.DataFrame({
        "impact_parameter_kpc": np.round(b_array, 3),
        "theta_total_arcsec": [f"{val:.4e}" if val < 1e-3 else f"{val:.4f}" for val in theta_tot],
        "theta_baryon_arcsec": [f"{val:.4e}" if val < 1e-3 else f"{val:.4f}" for val in theta_bar],
        "theta_scaffold_arcsec": [f"{val:.4e}" if val < 1e-3 else f"{val:.4f}" for val in theta_scaff],
        "convergence_kappa": [f"{val:.6e}" for val in kappa_tot],
        "inferred_M_DM_solar": [f"{val:.2e}" for val in M_DM_inf],
        "dark_to_light_ratio": np.round(M_DM_inf / M_baryon, 2)
    })

    # Metrics at dynamic crossover focus distance
    idx_focus = np.argmin(np.abs(b_array - b_focus))
    
    metrics = {
        "Target Baseline Acceleration (a_Ω)": f"{a_omega_used:.3e} m/s²",
        "Local Expansion Rate (H_local)": f"{h_local:.2f} km/s/Mpc",
        "1:1 Enclosed Mass Radius (b_cross)": f"{b_focus:.2f} kpc",
        "1:1 Total Mass Horizon (b_1:1)": f"{b_1to1_total:.2f} kpc",
        "Deflection @ Crossover": f"{theta_tot[idx_focus]:.4e} arcsec" if theta_tot[idx_focus] < 1e-2 else f"{theta_tot[idx_focus]:.2f} arcsec",
        "Total Dark/Light @ b_cross": f"{(M_DM_inf[idx_focus] / M_baryon):.2f} : 1"
    }

    return fig, df_lensing, metrics