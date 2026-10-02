"""
CPU Bandwidth Throttling & Gravitational Time Dilation Module
-------------------------------------------------------------
Substrate Logistics Cosmological Dashboard - Domain 4
Module: cpu_bandwidth_throttling.py

Physical Mechanics:
- Replaces continuous spacetime warping with metric CPU bandwidth throttling.
- Governed by the Pythagorean Bandwidth Allocation: c^2 = v^2 + u^2
- External clock dilation is evaluated from dynamic loads: u_clock = sqrt(1.0 - (g_Omega^2 + v_Omega^2))
- Disambiguates Absolute Deep Space vs. Earth Ground Observer reference frames.
- Calculates relative clock frequency shifts Delta_f / f_0 and daily clock drift (us/day).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from modules.data_loader import C_LIGHT, K_OMEGA, M_PROTON, M_SOLAR_KG

# ==========================================
# CELESTIAL & ORBITAL PRESETS
# ==========================================
_sag_mass = 4.154e6 * M_SOLAR_KG
_sag_nodes = _sag_mass / M_PROTON
_sag_reh = (2.0 * _sag_nodes * K_OMEGA) / (C_LIGHT ** 2)

THROTTLING_PRESETS = {
    "Earth Surface (Ground Observer)": {
        "mass_kg": 5.9722e24,
        "radius_m": 6.3712e6,
        "velocity_m_s": 465.1,  # Earth surface rotation at equator
        "alpha_tension": 0.0,
        "desc": "Terrestrial ground clock at equator (r = 6,371.2 km)"
    },
    "GPS Satellite Orbit (20,180 km Altitude)": {
        "mass_kg": 5.9722e24,
        "radius_m": 26.560e6,  # 6371.2 + 20180 km
        "velocity_m_s": 3874.0,  # Circular MEO orbital speed
        "alpha_tension": 0.0,
        "desc": "MEO atomic clock on GPS constellation satellite (+38.59 us/day vs ground)"
    },
    "ISS Orbit (400 km Altitude)": {
        "mass_kg": 5.9722e24,
        "radius_m": 6.7712e6,  # 6371.2 + 400 km
        "velocity_m_s": 7670.0,  # Low Earth Orbit speed
        "alpha_tension": 0.0,
        "desc": "LEO atomic clock on International Space Station (-24.62 us/day vs ground)"
    },
    "Solar Photosphere (Sun Surface)": {
        "mass_kg": 1.98847e30,
        "radius_m": 6.9634e8,
        "velocity_m_s": 0.0,
        "alpha_tension": 0.0,
        "desc": "Photospheric surface of the Sun (-183.14 ms/day gravitational redshift)"
    },
    "Sirius B (White Dwarf)": {
        "mass_kg": 1.018 * M_SOLAR_KG,
        "radius_m": 5.800e6,
        "velocity_m_s": 0.0,
        "alpha_tension": 0.0,
        "desc": "Degenerate white dwarf high-density surface (-22.39 s/day drift)"
    },
    "Crab Pulsar (Neutron Star)": {
        "mass_kg": 1.4 * M_SOLAR_KG,
        "radius_m": 12.0e3,
        "velocity_m_s": 0.0,
        "alpha_tension": 0.0,
        "desc": "Relativistic nuclear density surface (u_clock = 0.8097)"
    },
    "Sagittarius A* Event Horizon (R_EH)": {
        "mass_kg": _sag_mass,
        "radius_m": _sag_reh,  # Exact Parameter-Free R_EH Boundary
        "velocity_m_s": 0.0,
        "alpha_tension": 0.0,
        "desc": "Supermassive black hole exact Event Horizon (100% CPU Kernel Lock, u_clock = 0.00)"
    }
}


def compute_bandwidth_throttling_point(mass_kg, radius_m, velocity_m_s=0.0, alpha_tension=0.0):
    """
    Evaluates metric system bandwidth allocations and observable clock rate at a coordinate location.
    """
    node_count = mass_kg / M_PROTON
    
    # Normalized external dynamic loads
    v_escape = np.sqrt(2.0 * node_count * K_OMEGA / radius_m)
    g_omega_sq = (v_escape / C_LIGHT) ** 2     # Gravitational address deletion load (2GM / c^2 R)
    v_omega_sq = (velocity_m_s / C_LIGHT) ** 2 # Kinematic routing load (v^2 / c^2)
    alpha_omega_sq = alpha_tension ** 2        # Static internal structural tension load
    
    omega_ext = g_omega_sq + v_omega_sq
    omega_total = omega_ext + alpha_omega_sq
    
    # Observable clock update rate driven strictly by external dynamic loads
    if omega_ext >= 1.0:
        u_clock = 0.0
        is_locked = True
    else:
        u_clock = np.sqrt(1.0 - omega_ext)
        is_locked = False

    u_total = np.sqrt(max(0.0, 1.0 - omega_total))

    # Absolute relative frequency shift Delta_f / f_0 relative to unthrottled deep space (u = 1.0)
    delta_f_space = u_clock - 1.0
    drift_us_space = delta_f_space * 86400.0 * 1e6

    return {
        "N_nodes": node_count,
        "g_omega_sq": g_omega_sq,
        "v_omega_sq": v_omega_sq,
        "alpha_omega_sq": alpha_omega_sq,
        "omega_ext": omega_ext,
        "omega_total": omega_total,
        "u_clock": u_clock,
        "u_total": u_total,
        "delta_f_space": delta_f_space,
        "drift_us_space": drift_us_space,
        "is_locked": is_locked
    }


def compute_bandwidth_throttling_profile(mass_kg, r_surface_m, observer_mode="Static (v = 0)", r_max_ratio=10.0, n_pts=250):
    """
    Generates a 1D radial profile array of bandwidth allocations from surface out to r_max_ratio.
    """
    r_arr = np.geomspace(r_surface_m, r_surface_m * r_max_ratio, n_pts)
    node_count = mass_kg / M_PROTON

    g_sq_arr = (2.0 * node_count * K_OMEGA) / (C_LIGHT**2 * r_arr)
    
    if "Orbital" in observer_mode:
        v_orbit_arr = np.sqrt(node_count * K_OMEGA / r_arr)
        v_sq_arr = (v_orbit_arr / C_LIGHT)**2
    else:
        v_sq_arr = np.zeros_like(r_arr)
    
    omega_ext_arr = g_sq_arr + v_sq_arr
    u_clock_arr = np.sqrt(np.maximum(0.0, 1.0 - omega_ext_arr))
    
    # Shifts & Drifts vs Unthrottled Deep Space (u = 1.0)
    delta_f_space = u_clock_arr - 1.0
    drift_us_space = delta_f_space * 86400.0 * 1e6
    
    # Shifts & Drifts vs Earth Ground Observer (u_earth)
    ground_point = compute_bandwidth_throttling_point(5.9722e24, 6.3712e6, 465.1)
    u_earth = ground_point["u_clock"]
    
    delta_f_earth = (u_clock_arr - u_earth) / u_earth
    drift_us_earth = delta_f_earth * 86400.0 * 1e6

    df_profile = pd.DataFrame({
        "Radius [m]": r_arr,
        "Altitude [km]": (r_arr - r_surface_m) / 1e3,
        "r / R_surface": r_arr / r_surface_m,
        "g_Omega^2": g_sq_arr,
        "v_Omega^2": v_sq_arr,
        "External Load Omega": omega_ext_arr,
        "Observable u_clock": u_clock_arr,
        "Delta_f / f0 (vs Space)": delta_f_space,
        "Drift vs Space [μs/day]": drift_us_space,
        "Delta_f / f0 (vs Earth)": delta_f_earth,
        "Drift vs Earth [μs/day]": drift_us_earth
    })

    return df_profile


def run_cpu_bandwidth_throttling_analysis(
    preset_key="Earth Surface (Ground Observer)",
    custom_mass_kg=None,
    custom_radius_m=None,
    custom_velocity_m_s=None,
    alpha_tension=0.0,
    observer_mode="Static (v = 0)",
    r_max_ratio=10.0
):
    """
    Primary Streamlit UI execution driver function for CPU Bandwidth Throttling.
    """
    if preset_key in THROTTLING_PRESETS and custom_mass_kg is None:
        p = THROTTLING_PRESETS[preset_key]
        mass_kg = p["mass_kg"]
        radius_m = p["radius_m"]
        velocity_m_s = p["velocity_m_s"]
        title_label = preset_key
    else:
        mass_kg = custom_mass_kg if custom_mass_kg is not None else 5.9722e24
        radius_m = custom_radius_m if custom_radius_m is not None else 6.3712e6
        velocity_m_s = custom_velocity_m_s if custom_velocity_m_s is not None else 0.0
        title_label = f"Custom Target (M={mass_kg:.3e} kg, R={radius_m/1e3:.1f} km)"

    pt_res = compute_bandwidth_throttling_point(mass_kg, radius_m, velocity_m_s, alpha_tension)

    # Calculate Earth surface ground clock baseline for comparison
    ground_res = compute_bandwidth_throttling_point(5.9722e24, 6.3712e6, 465.1, 0.0)
    
    # Relative shift vs Earth ground observer: (u_target - u_ground) / u_ground
    if not pt_res['is_locked']:
        rel_delta_f_earth = (pt_res["u_clock"] - ground_res["u_clock"]) / ground_res["u_clock"]
        relative_drift_vs_earth_us = (pt_res["drift_us_space"] - ground_res["drift_us_space"])
    else:
        rel_delta_f_earth = -1.0
        relative_drift_vs_earth_us = -86400.0 * 1e6

    df_profile = compute_bandwidth_throttling_profile(
        mass_kg, radius_m, observer_mode=observer_mode, r_max_ratio=r_max_ratio
    )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), facecolor='#0f172a')
    ax1.set_facecolor('#070b14')
    ax2.set_facecolor('#070b14')

    # PANEL 1: System Resource Budget Allocation (Stacked Bar Chart)
    categories = ['Target System']
    g_val = min(1.0, pt_res['g_omega_sq'])
    v_val = min(1.0 - g_val, pt_res['v_omega_sq'])
    a_val = min(1.0 - g_val - v_val, pt_res['alpha_omega_sq'])
    u_val = pt_res['u_total']
    
    ax1.bar(categories, [g_val], label=r'Gravity Load $g_\Omega^2$', color='#ef4444', edgecolor='none')
    ax1.bar(categories, [v_val], bottom=[g_val], label=r'Kinematic Load $v_\Omega^2$', color='#3b82f6', edgecolor='none')
    ax1.bar(categories, [a_val], bottom=[g_val + v_val], label=r'Internal Tension $\alpha_\Omega^2$', color='#eab308', edgecolor='none')
    ax1.bar(categories, [u_val], bottom=[g_val + v_val + a_val], label=r'Residual CPU Budget $u_\Omega$', color='#10b981', edgecolor='none')

    ax1.set_ylim(0, 1.05)
    ax1.set_ylabel(r"Normalized Bandwidth Allocation ($\Omega \leq 1.0$)", color='#94a3b8', fontsize=11, fontweight='bold')
    ax1.set_title(f"CPU Bandwidth Allocation\n[{title_label}]", color='#f8fafc', fontsize=12, fontweight='bold', pad=12)
    ax1.tick_params(colors='#94a3b8', labelsize=10)
    ax1.grid(True, linestyle=':', alpha=0.2, color='#38bdf8')
    ax1.legend(loc='center right', facecolor='#0f172a', edgecolor='#38bdf8', labelcolor='#f8fafc', fontsize=9)

    if pt_res['is_locked']:
        status_text = (
            "KERNEL LOCK EVENT HORIZON\n"
            "----------------------------\n"
            "• Observable Clock Rate (u): 0.000000\n"
            "• Gravitational Load (g_Ω²): 1.000000\n"
            "• Clock Dephasing: 100% CPU Freeze"
        )
        ax1.text(
            0.05, 0.45, status_text, transform=ax1.transAxes,
            color='#f43f5e', fontsize=9, fontfamily='monospace', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.6', facecolor='#0f172a', edgecolor='#f43f5e', alpha=0.90)
        )
    elif pt_res['omega_ext'] < 0.01:
        annotation_text = (
            f"Bandwidth Breakdown:\n"
            f"-------------------\n"
            f"• Observable Clock Rate (u_clock): {pt_res['u_clock']:.10f}\n"
            f"• Gravity Load (g_Ω²): {g_val:.3e}\n"
            f"• Kinematic Load (v_Ω²): {v_val:.3e}\n"
            f"• Shift vs Earth Ground: {rel_delta_f_earth:+.4e}"
        )
        ax1.text(
            0.05, 0.50, annotation_text, transform=ax1.transAxes,
            color='#38bdf8', fontsize=9, fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.6', facecolor='#0f172a', edgecolor='#38bdf8', alpha=0.85)
        )

    # PANEL 2: Radial Clock Dephasing & Daily Drift Profile
    r_ratios = df_profile['r / R_surface']
    drifts_us = df_profile['Drift vs Space [μs/day]']

    ax2.plot(r_ratios, drifts_us, color='#38bdf8', linewidth=2.5, label=f'Clock Drift vs Space ({observer_mode})')
    ax2.axhline(0, color='#94a3b8', linestyle='--', alpha=0.5, label='Deep Space Unthrottled Baseline')

    target_drift_us = pt_res['drift_us_space']
    ax2.scatter([1.0], [target_drift_us], color='#f43f5e', s=100, zorder=5, label='Evaluation Surface')

    if pt_res['is_locked']:
        ax2.axhline(-8.64e7, color='#f43f5e', linestyle=':', alpha=0.7, label='100% Kernel Lock Freeze Floor')

    ax2.set_xlabel(r"Radial Distance ($r / R_{\mathrm{surface}}$)", color='#94a3b8', fontsize=11, fontweight='bold')
    ax2.set_ylabel(r"Daily Clock Drift [$\mu$s / day]", color='#94a3b8', fontsize=11, fontweight='bold')
    ax2.set_title("Radial Clock Frequency Dephasing Profile", color='#f8fafc', fontsize=12, fontweight='bold', pad=12)
    ax2.tick_params(colors='#94a3b8', labelsize=10)
    ax2.grid(True, linestyle=':', alpha=0.2, color='#38bdf8')
    ax2.legend(loc='lower right', facecolor='#0f172a', edgecolor='#38bdf8', labelcolor='#f8fafc', fontsize=9)

    plt.tight_layout()

    status_str = "🔒 KERNEL LOCK (0.00% CPU)" if pt_res['is_locked'] else f"🟢 OPERATIONAL ({pt_res['u_clock']*100:.6f}%)"
    
    metrics_dict = {
        "System Status": status_str,
        "Integer Node Load (N)": f"{pt_res['N_nodes']:.4e} nodes",
        "Gravitational Load (g_Ω²)": f"{pt_res['g_omega_sq']:.4e}",
        "Kinematic Load (v_Ω²)": f"{pt_res['v_omega_sq']:.4e}",
        "Shift vs Deep Space (Δf/f₀)": f"{pt_res['delta_f_space']:.4e}",
        "Shift vs Earth Ground (Δf/f₀)": f"{rel_delta_f_earth:+.4e}" if not pt_res['is_locked'] else "-1.0000e+00",
        "Relative Drift vs Earth": f"{relative_drift_vs_earth_us:+.2f} μs/day" if not pt_res['is_locked'] else "-86,400.00 s/day (Frozen)"
    }

    return fig, df_profile, metrics_dict