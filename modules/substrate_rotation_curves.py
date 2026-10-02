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
    compute_a_omega,
    get_compiled_node_density,
    get_hlocal,
    format_baryonic_mass
)

# Preset Benchmark Spiral & Dwarf Galaxies
ROTATION_PRESETS = {
    "NGC 3198 (Classic Spiral Benchmark)": {
        "M_baryon_solar": 3.2e10,      # 32 Billion M_sun (Disk + Gas)
        "R_d_kpc": 2.72,               # SPARC empirical exponential disk scale length
        "r_min_kpc": 0.0,              # Starts smoothly at 0.0 km/s
        "r_max_kpc": 35.0,
        "r_focus_kpc": 20.0,
        "desc": "Textbook SPARC spiral galaxy exhibiting flat outer rotation out to 35 kpc."
    },
    "Milky Way (Stellar Disk + Bulge + Gas)": {
        "M_baryon_solar": 8.8e10,      # 88 Billion M_sun (Full visible baryonic budget)
        "R_d_kpc": 2.30,               # Empirical Milky Way scale length
        "r_min_kpc": 0.0,
        "r_max_kpc": 30.0,
        "r_focus_kpc": 8.12,
        "desc": "Milky Way rotation curve using full baryonic mass (Bulge + Bar + Disk + Gas Ring)."
},
    "DDO 154 (Gas-Rich Dwarf Galaxy)": {
        "M_baryon_solar": 3.8e8,       # 380 Million M_sun
        "R_d_kpc": 1.20,               # Dwarf galaxy scale length
        "r_min_kpc": 0.0,
        "r_max_kpc": 10.0,
        "r_focus_kpc": 5.0,
        "desc": "Extreme gas-dominated dwarf galaxy where dark matter halo fits fail without extreme tuning."
    }
}

def calculate_rotation_profile(M_baryon_solar, r_kpc_array, a_omega_active=A_OMEGA, R_d_kpc=None):
    """
    Computes Substrate rotation curves with geometric underflow coupling:
    g_eff = sqrt(g_baryon^2 + g_baryon * a_Omega)
    
    Supports optional exponential disk scale length R_d_kpc for extended mass M(r),
    with exact central acceleration limit g(0) = GM / (2 R_d^2).
    """
    N_total = (M_baryon_solar * M_SOLAR_KG) / M_PROTON
    
    if R_d_kpc is not None and R_d_kpc > 0:
        x = r_kpc_array / R_d_kpc
        R_d_m = R_d_kpc * KPC_TO_METER
        
        # Exact ratio = [1 - (1 + x) e^-x] / x^2
        # Taylor expansion around x=0: 1/2 - x/3 + x^2/8 - ...
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = (1.0 - (1.0 + x) * np.exp(-x)) / (x**2)
        ratio = np.where(x < 1e-5, 0.5 - x/3.0 + x**2/8.0, ratio)
        
        g_baryon_ms2 = (N_total * K_OMEGA / R_d_m**2) * ratio
    else:
        r_meters = r_kpc_array * KPC_TO_METER
        safe_r_meters = np.maximum(r_meters, 1e-10)
        g_baryon_ms2 = (N_total * K_OMEGA) / (safe_r_meters**2)

    # 1. Classical Keplerian Orbital Velocity: v_kepler = sqrt(r * g_baryon)
    r_meters = r_kpc_array * KPC_TO_METER
    v_kepler_ms = np.sqrt(r_meters * g_baryon_ms2)
    v_kepler_kms = v_kepler_ms / 1000.0
    
    # 2. Geometric Underflow Coupling: g_eff = sqrt(g_baryon^2 + g_baryon * a_Omega)
    g_eff_ms2 = np.sqrt(g_baryon_ms2**2 + g_baryon_ms2 * a_omega_active)
    
    # 3. Substrate Clamped Orbital Velocity: v_substrate = sqrt(r * g_eff) [m/s]
    v_substrate_ms = np.sqrt(r_meters * g_eff_ms2)
    v_substrate_kms = v_substrate_ms / 1000.0
    
    # Force exact 0.0 at r = 0 coordinate origin
    v_kepler_kms = np.where(r_kpc_array <= 0, 0.0, v_kepler_kms)
    v_substrate_kms = np.where(r_kpc_array <= 0, 0.0, v_substrate_kms)
    
    # 4. Theoretical Asymptotic Flat Floor: v_flat = (N_total * K_Omega * a_Omega)^(1/4)
    v_flat_ms = (N_total * K_OMEGA * a_omega_active)**(0.25)
    v_flat_kms = float(v_flat_ms / 1000.0)
    
    # 5. Transition Radius R*: Radius where Effective Substrate Acceleration crosses floor (g_eff = a_Omega)
    idx_star = np.argmin(np.abs(g_eff_ms2 - a_omega_active))
    transition_radius_kpc = float(r_kpc_array[idx_star])
    
    return v_substrate_kms, v_kepler_kms, g_baryon_ms2, g_eff_ms2, transition_radius_kpc, v_flat_kms

def run_rotation_curve_analysis(
    galaxy_target_name="NGC 3198 (Classic Spiral Benchmark)",
    custom_M_baryon=None,
    custom_R_d=None,
    h_global=H_GLOBAL,
    use_catalog_query=False,
    target_pos_mpc=None,
    sigma_mpc=1.8,
    tree=None,
    gal_positions=None,
    nodes_per_galaxy=None
):
    """
    Main Streamlit UI Execution Engine for Galactic Rotation Curves.
    """
    if use_catalog_query and target_pos_mpc is not None and tree is not None:
        rho_N_local = get_compiled_node_density(tree, gal_positions, nodes_per_galaxy, target_pos_mpc, sigma_mpc=sigma_mpc)
        h_local = get_hlocal(tree, gal_positions, nodes_per_galaxy, target_pos_mpc, sigma_mpc=sigma_mpc, h_global=h_global)
        a_omega_used = compute_a_omega(h_local)
        
        M_baryon = custom_M_baryon if custom_M_baryon else 5.0e10
        R_d_val = custom_R_d if custom_R_d else 2.5
        r_min, r_max, r_focus = 0.0, 30.0, 10.0
        target_title = f"{galaxy_target_name}\n(X={target_pos_mpc[0]:.1f}, Y={target_pos_mpc[1]:.1f}, Z={target_pos_mpc[2]:.1f} Mpc)"
    else:
        if galaxy_target_name in ROTATION_PRESETS:
            preset = ROTATION_PRESETS[galaxy_target_name]
            M_baryon = custom_M_baryon if custom_M_baryon else preset["M_baryon_solar"]
            R_d_val = custom_R_d if custom_R_d else preset.get("R_d_kpc", 2.5)
            r_min = preset.get("r_min_kpc", 0.0)
            r_max = preset["r_max_kpc"]
            r_focus = preset["r_focus_kpc"]
        else:
            M_baryon = custom_M_baryon if custom_M_baryon else 5.0e10
            R_d_val = custom_R_d if custom_R_d else 2.5
            r_min, r_max, r_focus = 0.0, 30.0, 10.0

        a_omega_used = compute_a_omega(h_global)
        h_local = h_global
        rho_N_local = 0.0
        target_title = galaxy_target_name

    r_array = np.linspace(r_min, r_max, 500)
    
    v_sub, v_kep, g_bar, g_eff, R_star_kpc, v_flat_kms = calculate_rotation_profile(
        M_baryon_solar=M_baryon,
        r_kpc_array=r_array,
        a_omega_active=a_omega_used,
        R_d_kpc=R_d_val
    )
    
    # Render Diagnostic Plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), facecolor='#0b1120')
    for ax in (ax1, ax2):
        ax.set_facecolor('#070c18')
        ax.tick_params(colors='#94a3b8')
        ax.grid(True, alpha=0.2, linestyle=':')

    # Panel 1: Rotation Curve v(r)
    ax1.plot(r_array, v_sub, color='#10b981', linewidth=2.5, label=r'Substrate Rotation Floor $v_{\mathrm{rot}}(r)$')
    ax1.plot(r_array, v_kep, color='#ef4444', linewidth=1.8, linestyle='--', label=r'Keplerian Baryonic Decay $v \propto 1/\sqrt{r}$')
    ax1.axhline(v_flat_kms, color='#f59e0b', linestyle=':', label=f'Asymptotic Flat Floor ({v_flat_kms:.1f} km/s)')
    ax1.axvline(R_star_kpc, color='#38bdf8', linestyle='-.', label=f'Transition Radius R* ({R_star_kpc:.1f} kpc)')
    
    ax1.set_xlabel('Galactocentric Radius $r$ [kpc]', color='#cbd5e1', fontsize=11)
    ax1.set_ylabel('Orbital Rotation Velocity $v$ [km/s]', color='#cbd5e1', fontsize=11)
    ax1.set_title(f'{target_title}\nGalactic Rotation Curve (Zero Dark Matter)', color='#f8fafc', fontsize=11)
    ax1.legend(loc='lower right', facecolor='#0f172a', edgecolor='none', labelcolor='#e2e8f0')

    # Panel 2: Acceleration Profile g(r) vs a_Omega Underflow Floor
    ax2.plot(r_array, g_eff, color='#10b981', linewidth=2.5, label=r'Substrate Effective Acceleration $g_{\mathrm{eff}}(r)$')
    ax2.plot(r_array, g_bar, color='#ef4444', linewidth=1.8, linestyle='--', label=r'Baryonic Acceleration $g_{\mathrm{baryon}}(r)$')
    ax2.axhline(a_omega_used, color='#a855f7', linestyle='--', linewidth=2.0, label=f'Kinematic Underflow Floor $a_\\Omega$ ({a_omega_used:.3e} m/s²)')
    
    lbl_R_star = r'R* Boundary ($g_{\mathrm{eff}} = a_\Omega$)'
    ax2.axvline(R_star_kpc, color='#38bdf8', linestyle='-.', label=lbl_R_star)

    ax2.set_yscale('log')
    ax2.set_xlabel('Galactocentric Radius $r$ [kpc]', color='#cbd5e1', fontsize=11)
    ax2.set_ylabel(r'Acceleration $g$ [m/s²]', color='#cbd5e1', fontsize=11)
    ax2.set_title(r'Kinematic Underflow Threshold: $g_{\mathrm{eff}}$ vs $g_{\mathrm{baryon}}$', color='#f8fafc', fontsize=11)
    ax2.legend(loc='upper right', facecolor='#0f172a', edgecolor='none', labelcolor='#e2e8f0')

    plt.tight_layout()

    # Export DataFrame with full double-precision floating point numbers
    df_rotation = pd.DataFrame({
        "radius_kpc": np.round(r_array, 2),
        "v_substrate_kms": np.round(v_sub, 2),
        "v_kepler_kms": np.round(v_kep, 2),
        "g_eff_ms2": g_eff,
        "g_baryon_ms2": g_bar,
        "a_omega_floor_ms2": a_omega_used
    })

    # Formatted copy for UI rendering to prevent zero-truncation in tables
    df_ui_display = df_rotation.copy()
    for col in ["g_eff_ms2", "g_baryon_ms2", "a_omega_floor_ms2"]:
        df_ui_display[col] = df_ui_display[col].map(lambda x: f"{x:.3e}")

    idx_focus = np.argmin(np.abs(r_array - r_focus))

    metrics = {
        "True Baryonic Mass": format_baryonic_mass(M_baryon),
        "Disk Scale Length (R_d)": f"{R_d_val:.2f} kpc",
        "Kinematic Underflow Floor (a_Ω)": f"{a_omega_used:.3e} m/s²",
        "Transition Radius (R*)": f"{R_star_kpc:.2f} kpc",
        "Flat Asymptotic Velocity": f"{v_flat_kms:.1f} km/s"
    }

    return fig, df_rotation, metrics