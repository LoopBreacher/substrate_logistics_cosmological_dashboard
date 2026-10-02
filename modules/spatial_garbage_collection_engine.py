import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .data_loader import (
    C_LIGHT,
    K_OMEGA,
    M_PROTON,
    MPC_TO_METER
)

# Bare-metal sub-node spatial pixel floor anchor (meters)
LAMBDA_P_BAR = 2.10309e-16
# Sub-node spatial voxel volume V_node = lambda_p_bar^3 (m^3)
V_NODE_M3 = LAMBDA_P_BAR ** 3  # 9.301e-48 m^3

# Speed of light Kernel Lock threshold ratio: v_lock = 1.0 * c
V_LOCK_RATIO = 1.0
V_LOCK_M_S = V_LOCK_RATIO * C_LIGHT  # 299,792,458 m/s

# ==========================================
# 1. PRESET CELESTIAL HARDWARE BENCHMARKS
# ==========================================
# Helper to compute exact Kernel Lock Horizon from mass
def _get_kernel_lock_radius(mass_kg):
    node_count_N = mass_kg / M_PROTON
    return (2.0 * node_count_N * K_OMEGA) / (C_LIGHT ** 2)

CELESTIAL_PRESETS = {
    "Earth": {
        "mass_kg": 5.9722e24,
        "radius_m": 6.3712e6,
        "description": "Terrestrial Boundary Layer (N = 3.570e+51 nodes)"
    },
    "Moon": {
        "mass_kg": 7.342e22,
        "radius_m": 1.7374e6,
        "description": "Lunar Boundary Layer (N = 4.389e+49 nodes)"
    },
    "Sun": {
        "mass_kg": 1.98847e30,
        "radius_m": 6.9634e8,
        "description": "Solar Boundary Layer (N = 1.189e+57 nodes)"
    },
    "Cygnus X-1": {
        "mass_kg": 4.2155e31,  # ~21.2 M_sun
        "radius_m": _get_kernel_lock_radius(4.2155e31),  # Exact 62.583 km
        "description": "Stellar-Mass Kernel Lock (N = 2.520e+58 nodes, 21.2 M_sun)"
    },
    "Sagittarius A*": {
        "mass_kg": 8.55e36,    # ~4.3 million M_sun
        "radius_m": _get_kernel_lock_radius(8.55e36),    # Exact 12.693 million km
        "description": "Galactic Core Kernel Lock (N = 5.109e+63 nodes, 4.3M M_sun)"
    },
    "M87*": {
        "mass_kg": 1.2925e40,  # ~6.5 billion M_sun
        "radius_m": _get_kernel_lock_radius(1.2925e40),  # Exact 128.26 AU (~19.188 billion km)
        "description": "EHT Supermassive Target (N = 7.727e+66 nodes, 6.5B M_sun)"
    }
}

# ==========================================
# 2. CORE SPATIAL GARBAGE COLLECTION SOLVER
# ==========================================
def compute_garbage_collection_field(
    mass_kg, surface_radius_m,
    r_max_ratio=10.0, n_pts=250
):
    """
    Computes spatial address garbage collection dynamics:
    1. Integer Informational Load N = M / m_p
    2. Surface and radial freefall acceleration g_omega(r) = N * K_omega / r^2
    3. Inflow / Escape velocity v_escape(r) = sqrt(2 * g_omega * r)
    4. Circular Low-Orbit Bypass speed v_orbit(r) = sqrt(g_omega * r) = v_escape / sqrt(2)
    5. Volumetric address purge rate dV_deleted/dt = 4 * pi * r^2 * v_inflow
    6. Discrete pixel purge rate dN_pixels/dt = dV_deleted / V_node
    7. Kernel Lock Safety Margin = max(0, 1.0 - v_surf/c) * 100%
    """
    node_count_N = mass_kg / M_PROTON
    
    # Kernel Lock Boundary Radius where v_inflow = c
    r_kernel_lock_m = (2.0 * node_count_N * K_OMEGA) / (V_LOCK_M_S ** 2)
    
    r_min = max(r_kernel_lock_m * 0.5, surface_radius_m * 0.1)
    r_max = surface_radius_m * r_max_ratio
    r_grid_m = np.geomspace(r_min, r_max, n_pts)
    
    g_omega_prof = (node_count_N * K_OMEGA) / (r_grid_m ** 2)
    v_inflow_m_s = np.sqrt(2.0 * g_omega_prof * r_grid_m)  # v_escape = v_inflow
    v_orbit_m_s = np.sqrt(g_omega_prof * r_grid_m)         # Circular orbit speed
    v_inflow_c_ratio = v_inflow_m_s / C_LIGHT
    
    vol_deleted_m3_s = 4.0 * np.pi * (r_grid_m ** 2) * v_inflow_m_s
    pixel_purge_rate_s = vol_deleted_m3_s / V_NODE_M3
    
    g_surface = (node_count_N * K_OMEGA) / (surface_radius_m ** 2)
    v_inflow_surface_m_s = np.sqrt(2.0 * g_surface * surface_radius_m)
    v_orbit_surface_m_s = np.sqrt(g_surface * surface_radius_m)
    v_surf_c_ratio = v_inflow_surface_m_s / C_LIGHT
    
    # Kernel Lock Safety Margin Calculations
    safety_margin_pct = max(0.0, (1.0 - v_surf_c_ratio) * 100.0)
    capacity_used_pct = min(100.0, v_surf_c_ratio * 100.0)
    
    vol_deleted_surface = 4.0 * np.pi * (surface_radius_m ** 2) * v_inflow_surface_m_s
    pixel_purge_surface = vol_deleted_surface / V_NODE_M3
    
    return {
        "node_count_N": node_count_N,
        "r_kernel_lock_m": r_kernel_lock_m,
        "g_surface": g_surface,
        "v_inflow_surface_m_s": v_inflow_surface_m_s,
        "v_orbit_surface_m_s": v_orbit_surface_m_s,
        "v_inflow_surface_c_ratio": v_surf_c_ratio,
        "safety_margin_pct": safety_margin_pct,
        "capacity_used_pct": capacity_used_pct,
        "vol_deleted_surface": vol_deleted_surface,
        "pixel_purge_surface": pixel_purge_surface,
        "r_grid_m": r_grid_m,
        "g_omega_prof": g_omega_prof,
        "v_inflow_m_s": v_inflow_m_s,
        "v_orbit_m_s": v_orbit_m_s,
        "v_inflow_c_ratio": v_inflow_c_ratio,
        "vol_deleted_m3_s": vol_deleted_m3_s,
        "pixel_purge_rate_s": pixel_purge_rate_s
    }

# ==========================================
# 3. MAIN UI EXECUTION & PLOTTING ENGINE
# ==========================================
def run_spatial_garbage_collection_analysis(
    preset_key="Earth",
    custom_mass_kg=None,
    custom_radius_m=None,
    r_max_ratio=10.0
):
    """
    Main Streamlit UI entry point for Spatial Garbage Collection Engine.
    Returns (fig, df_table, metrics).
    """
    if preset_key in CELESTIAL_PRESETS and custom_mass_kg is None:
        p = CELESTIAL_PRESETS[preset_key]
        mass_kg = p["mass_kg"]
        radius_m = p["radius_m"]
        target_label = f"{preset_key} ({p['description']})"
    else:
        mass_kg = custom_mass_kg if custom_mass_kg else 5.9722e24
        radius_m = custom_radius_m if custom_radius_m else 6.3712e6
        target_label = f"Custom Load (M = {mass_kg:.3e} kg, R = {radius_m:.3e} m)"

    res = compute_garbage_collection_field(mass_kg, radius_m, r_max_ratio=r_max_ratio)
    
    r_grid = res["r_grid_m"]
    r_ratio = r_grid / radius_m

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), facecolor='#0b1120')
    ax1.set_facecolor('#070c18')
    ax2.set_facecolor('#070c18')

    for ax in (ax1, ax2):
        ax.tick_params(colors='#94a3b8')
        ax.grid(True, alpha=0.2, linestyle=':')

    # Panel 1: Dynamic Y-Scaling for Inflow & Escape Velocity Profile
    v_inflow_kms = res["v_inflow_m_s"] / 1000.0
    v_max_kms = float(np.max(v_inflow_kms))
    c_kms = C_LIGHT / 1000.0
    surf_radius_km = radius_m / 1000.0

    ax1.plot(r_ratio, v_inflow_kms, color='#38bdf8', linewidth=2.5, 
             label=r'Inflow / Escape Velocity $v_{\mathrm{escape}}(r)$ [km/s]')
    ax1.axvline(1.0, color='#10b981', linestyle='--', linewidth=1.8, 
                label=f'Physical Surface Boundary ($R = {surf_radius_km:,.1f}$ km)')

    # Dynamic scaling for non-black hole targets vs near-light speed targets
    if v_max_kms < 0.5 * c_kms:
        ax1.set_ylim(0, v_max_kms * 1.18)
        badge_text = (f"Kernel Lock Limit: 1.0 c ({c_kms:,.0f} km/s)\n"
                      f"Safety Margin: {res['safety_margin_pct']:.4f}%\n"
                      f"Capacity Used: {res['capacity_used_pct']:.6f}%")
        ax1.text(0.96, 0.94, badge_text, transform=ax1.transAxes,
                 ha='right', va='top', fontsize=8.5, fontweight='bold',
                 color='#38bdf8',
                 bbox=dict(boxstyle='round,pad=0.5', facecolor='#0f172a', edgecolor='#38bdf8', alpha=0.85))
        legend_loc = 'center right'
    else:
        ax1.set_ylim(0, max(v_max_kms, c_kms) * 1.08)
        ax1.axhline(c_kms, color='#f43f5e', linestyle=':', linewidth=1.8, 
                    label='Kernel Lock Commit Limit ($1.0 \\cdot c$)')
        if r_grid[0] <= res["r_kernel_lock_m"] <= r_grid[-1]:
            r_lock_val = res["r_kernel_lock_m"]
            ax1.axvline(res["r_kernel_lock_m"] / radius_m, color='#f43f5e', linestyle='-.', linewidth=1.8, 
                        label=f'Kernel Lock Horizon ($R = {r_lock_val:.3e}$ m)')
        legend_loc = 'upper right'

    ax1.set_xlabel('Normalized Radial Distance $r / R_{\mathrm{surface}}$', color='#cbd5e1', fontsize=11)
    ax1.set_ylabel(r'Velocity $v$ [km/s]', color='#cbd5e1', fontsize=11)
    ax1.set_title(f'Spatial Coordinate Address Inflow: {target_label}', color='#f8fafc', fontsize=12)
    ax1.legend(loc=legend_loc, facecolor='#0f172a', edgecolor='none', fontsize=8.5)

    # Panel 2: Discrete Spatial Pixel Purge Rate
    ax2.plot(r_ratio, res["pixel_purge_rate_s"], color='#a855f7', linewidth=2.5, 
             label=r'Discrete Address Purge Rate $\dot{N}_{\mathrm{pixels}}(r)$ [pixels/s]')
    ax2.axvline(1.0, color='#10b981', linestyle='--', linewidth=1.8, label='Physical Surface Boundary')
    ax2.set_yscale('log')

    ax2.set_xlabel('Normalized Radial Distance $r / R_{\mathrm{surface}}$', color='#cbd5e1', fontsize=11)
    ax2.set_ylabel(r'Spatial Address Deletion Rate [pixels / s]', color='#cbd5e1', fontsize=11)
    ax2.set_title(r'Substrate Memory Garbage Collection Intensity ($E=0$)', color='#f8fafc', fontsize=12)
    ax2.legend(loc='upper right', facecolor='#0f172a', edgecolor='none', fontsize=8.5)

    plt.tight_layout()

    df_table = pd.DataFrame({
        "r_meters": np.round(r_grid, 2),
        "r_ratio_Rsurface": np.round(r_ratio, 3),
        "g_omega_m_s2": np.round(res["g_omega_prof"], 4),
        "v_escape_kms": np.round(res["v_inflow_m_s"] / 1000.0, 3),
        "v_orbit_kms": np.round(res["v_orbit_m_s"] / 1000.0, 3),
        "v_inflow_c_ratio": res["v_inflow_c_ratio"],
        "vol_deleted_m3_s": res["vol_deleted_m3_s"],
        "N_pixels_purged_per_sec": res["pixel_purge_rate_s"]
    })

    if res['v_inflow_surface_c_ratio'] < 1.0:
        safety_str = f"{res['safety_margin_pct']:.4f}% ({res['capacity_used_pct']:.4f}% Capacity)"
    else:
        safety_str = "0.0000% (KERNEL LOCK ACTIVE)"

    metrics = {
        "Integer Node Load (N)": f"{res['node_count_N']:.4e} nodes",
        "Surface Acceleration (g_Ω)": f"{res['g_surface']:.4f} m/s²",
        "Escape Velocity (v_escape)": f"{res['v_inflow_surface_m_s']/1000.0:.3f} km/s ({res['v_inflow_surface_c_ratio']:.4e} c)",
        "Low-Orbit Speed (v_orbit)": f"{res['v_orbit_surface_m_s']/1000.0:.3f} km/s",
        "Surface Pixel Purge Rate": f"{res['pixel_purge_surface']:.4e} pixels/s",
        "Kernel Lock Safety Margin": safety_str
    }

    return fig, df_table, metrics