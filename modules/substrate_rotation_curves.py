import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import i0, i1, k0, k1

from .data_loader import (
    C_LIGHT,
    H_GLOBAL,
    K_OMEGA,
    M_PROTON,
    M_SOLAR_KG,
    KPC_TO_METER,
    A_OMEGA,
    DISK_TRUNCATION_RATIO,
    ETA_POCKET,
    G_3PDELTA_MS2,
    SIGMA_POCKET_GALAXY_KGM2,
    SIGMA_POCKET_GALAXY_MSUN_PC2,
    A_BTFR_POCKET,
    compute_a_omega,
    get_compiled_node_density,
    get_hlocal,
    format_baryonic_mass,
    compute_locked_pocket_r0_m,
    burkert_3d_accel,
)

# Preset Benchmark Spiral, Giant Spiral & Dwarf Galaxies
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
    },
    "UGC 2885 (Giant Spiral Benchmark)": {
        "M_baryon_solar": 3.5e11,      # 350 Billion M_sun
        "R_d_kpc": 11.50,              # Rubin+ 1980 / Hunter+ 2013 giant LSB spiral scale length
        "r_min_kpc": 0.0,
        "r_max_kpc": 70.0,
        "r_focus_kpc": 50.0,
        "desc": "Massive giant spiral testing high-mass BTFR pocket lock out to 70 kpc."
    },
}


def _compute_baryon_accel_ms2(M_baryon_solar, r_kpc_array, R_d_kpc=None, disk_geometry="Freeman Thin Disk"):
    """Computes pure linear K_Omega baryonic radial acceleration g_baryon(r) [m/s^2]."""
    N_total = (M_baryon_solar * M_SOLAR_KG) / M_PROTON
    r_meters = np.asarray(r_kpc_array, dtype=float) * KPC_TO_METER

    if R_d_kpc is not None and R_d_kpc > 0:
        x = np.asarray(r_kpc_array, dtype=float) / float(R_d_kpc)
        R_d_m = float(R_d_kpc) * KPC_TO_METER

        if disk_geometry and str(disk_geometry).startswith("Freeman"):
            y = 0.5 * x
            y_safe = np.clip(y, 1e-6, 300.0)
            bessel_fac = 0.5 * x * (i0(y_safe) * k0(y_safe) - i1(y_safe) * k1(y_safe))
            asym_fac = np.where(x > 0.0, 1.0 / (np.maximum(x, 1e-12) ** 2), 0.0)
            ratio = np.where(x < 1e-4, 0.0, np.where(y > 300.0, asym_fac, bessel_fac))
        else:
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = (1.0 - (1.0 + x) * np.exp(-x)) / (x ** 2)
            ratio = np.where(x < 1e-5, 0.5 - x / 3.0 + (x ** 2) / 8.0, ratio)

        return (N_total * K_OMEGA / (R_d_m ** 2)) * ratio

    safe_r_meters = np.maximum(r_meters, 1e-10)
    return (N_total * K_OMEGA) / (safe_r_meters ** 2)


def calculate_rotation_profile(
    M_baryon_solar,
    r_kpc_array,
    a_omega_active=A_OMEGA,
    R_d_kpc=None,
    disk_geometry="Freeman Thin Disk",
    gravity_model="Legacy Underflow (a_Ω)",
):
    """
    Computes Substrate rotation curves, preserving the 6-tuple signature:
        (v_primary_kms, v_kepler_kms, g_baryon_ms2, g_primary_ms2, transition_radius_kpc, v_flat_kms)
    Supports both Legacy Kinematic Underflow (a_Omega) and Topological Cover Pocket
    (Linear K_Omega + 3 P_delta Burkert).
    """
    r_arr = np.asarray(r_kpc_array, dtype=float)
    r_meters = r_arr * KPC_TO_METER
    M_baryon_kg = float(M_baryon_solar) * M_SOLAR_KG
    N_total = M_baryon_kg / M_PROTON

    g_baryon_ms2 = _compute_baryon_accel_ms2(M_baryon_solar, r_arr, R_d_kpc=R_d_kpc, disk_geometry=disk_geometry)
    v_kepler_kms = np.where(r_arr <= 0.0, 0.0, np.sqrt(np.maximum(r_meters * g_baryon_ms2, 0.0)) / 1000.0)

    if gravity_model and "Cover Pocket" in str(gravity_model):
        r0_m = compute_locked_pocket_r0_m(M_baryon_kg, SIGMA_POCKET_GALAXY_KGM2)
        g_pocket_ms2 = burkert_3d_accel(r_meters, r0_m, SIGMA_POCKET_GALAXY_KGM2)
        g_primary_ms2 = g_baryon_ms2 + g_pocket_ms2
        v_primary_kms = np.where(r_arr <= 0.0, 0.0, np.sqrt(np.maximum(r_meters * g_primary_ms2, 0.0)) / 1000.0)
        v_flat_kms = float(((N_total * K_OMEGA * A_BTFR_POCKET) ** 0.25) / 1000.0)
        r0_kpc = r0_m / KPC_TO_METER
        return v_primary_kms, v_kepler_kms, g_baryon_ms2, g_primary_ms2, r0_kpc, v_flat_kms

    # Legacy Kinematic Underflow (a_Omega)
    g_eff_ms2 = np.sqrt((g_baryon_ms2 ** 2) + g_baryon_ms2 * a_omega_active)
    v_substrate_kms = np.where(r_arr <= 0.0, 0.0, np.sqrt(np.maximum(r_meters * g_eff_ms2, 0.0)) / 1000.0)
    v_flat_kms = float(((N_total * K_OMEGA * a_omega_active) ** 0.25) / 1000.0)

    if np.max(g_eff_ms2) < a_omega_active:
        transition_radius_kpc = None
    else:
        idx_peak = int(np.argmax(g_eff_ms2))
        idx_star = idx_peak + int(np.argmin(np.abs(g_eff_ms2[idx_peak:] - a_omega_active)))
        transition_radius_kpc = float(r_arr[idx_star])

    return v_substrate_kms, v_kepler_kms, g_baryon_ms2, g_eff_ms2, transition_radius_kpc, v_flat_kms


def calculate_cover_pocket_rotation_profile(
    M_baryon_solar, r_kpc_array, a_omega_active=A_OMEGA, R_d_kpc=None, disk_geometry="Freeman Thin Disk"
):
    """
    Evaluates BOTH Topological Cover Pocket (Linear K_Omega + 3 P_delta Burkert)
    and Legacy Kinematic Underflow (a_Omega) profiles simultaneously for comparative overlay.
    """
    r_arr = np.asarray(r_kpc_array, dtype=float)
    r_meters = r_arr * KPC_TO_METER
    M_baryon_kg = float(M_baryon_solar) * M_SOLAR_KG
    N_total = M_baryon_kg / M_PROTON

    v_legacy_kms, v_kepler_kms, g_baryon_ms2, g_legacy_eff_ms2, R_star_kpc, v_flat_legacy_kms = calculate_rotation_profile(
        M_baryon_solar=M_baryon_solar,
        r_kpc_array=r_arr,
        a_omega_active=a_omega_active,
        R_d_kpc=R_d_kpc,
        disk_geometry=disk_geometry,
        gravity_model="Legacy Underflow (a_Ω)",
    )

    r0_m = compute_locked_pocket_r0_m(M_baryon_kg, SIGMA_POCKET_GALAXY_KGM2)
    r0_kpc = float(r0_m / KPC_TO_METER)
    g_pocket_ms2 = burkert_3d_accel(r_meters, r0_m, SIGMA_POCKET_GALAXY_KGM2)
    g_pocket_tot_ms2 = g_baryon_ms2 + g_pocket_ms2

    v_pocket_only_kms = np.where(r_arr <= 0.0, 0.0, np.sqrt(np.maximum(r_meters * g_pocket_ms2, 0.0)) / 1000.0)
    v_pocket_kms = np.where(r_arr <= 0.0, 0.0, np.sqrt(np.maximum(r_meters * g_pocket_tot_ms2, 0.0)) / 1000.0)
    v_btfr_pocket_kms = float(((N_total * K_OMEGA * A_BTFR_POCKET) ** 0.25) / 1000.0)

    return {
        "v_pocket_kms": v_pocket_kms,
        "v_pocket_only_kms": v_pocket_only_kms,
        "v_legacy_aOmega_kms": v_legacy_kms,
        "v_kepler_kms": v_kepler_kms,
        "g_baryon_ms2": g_baryon_ms2,
        "g_pocket_ms2": g_pocket_ms2,
        "g_pocket_tot_ms2": g_pocket_tot_ms2,
        "g_legacy_eff_ms2": g_legacy_eff_ms2,
        "r0_kpc": r0_kpc,
        "R_star_kpc": R_star_kpc,
        "v_flat_legacy_kms": v_flat_legacy_kms,
        "v_btfr_pocket_kms": v_btfr_pocket_kms,
    }


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
    nodes_per_galaxy=None,
    disk_geometry="Freeman Thin Disk",
    gravity_model="Cover Pocket (Linear K_Ω + 3 P_δ Burkert)",
):
    """
    Main Streamlit UI Execution Engine for Galactic Rotation Curves.
    Supports both Topological Cover Pocket (Linear K_Omega + 3 P_delta Burkert) and
    Legacy Kinematic Underflow (a_Omega), overlaying the inactive model for direct comparison.
    """
    if use_catalog_query and target_pos_mpc is not None and tree is not None:
        _ = get_compiled_node_density(tree, gal_positions, nodes_per_galaxy, target_pos_mpc, sigma_mpc=sigma_mpc)
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
        target_title = galaxy_target_name

    r_array = np.linspace(r_min, r_max, 500)
    is_cover_pocket = bool(gravity_model and "Cover Pocket" in str(gravity_model))

    prof = calculate_cover_pocket_rotation_profile(
        M_baryon_solar=M_baryon,
        r_kpc_array=r_array,
        a_omega_active=a_omega_used,
        R_d_kpc=R_d_val,
        disk_geometry=disk_geometry,
    )

    v_pocket = prof["v_pocket_kms"]
    v_pocket_only = prof["v_pocket_only_kms"]
    v_legacy = prof["v_legacy_aOmega_kms"]
    v_kep = prof["v_kepler_kms"]
    g_bar = prof["g_baryon_ms2"]
    g_pocket = prof["g_pocket_ms2"]
    g_pocket_tot = prof["g_pocket_tot_ms2"]
    g_legacy_eff = prof["g_legacy_eff_ms2"]
    r0_kpc = prof["r0_kpc"]
    R_star_kpc = prof["R_star_kpc"]
    v_flat_legacy = prof["v_flat_legacy_kms"]
    r_trunc_kpc = DISK_TRUNCATION_RATIO * R_d_val if (R_d_val is not None and R_d_val > 0) else None

    v_primary = v_pocket if is_cover_pocket else v_legacy
    v_alt = v_legacy if is_cover_pocket else v_pocket
    g_primary = g_pocket_tot if is_cover_pocket else g_legacy_eff
    g_alt = g_legacy_eff if is_cover_pocket else g_pocket_tot

    lbl_primary = (
        rf"Cover Pocket Total (Linear $\mathcal{{K}}_\Omega$ + $3P_\delta$, {disk_geometry})"
        if is_cover_pocket
        else rf"Legacy Underflow ($a_\Omega$ Floor, {disk_geometry})"
    )
    lbl_alt = (
        rf"Legacy Underflow Overlay ($a_\Omega$, $v_{{\mathrm{{flat}}}}={v_flat_legacy:.1f}$ km/s)"
        if is_cover_pocket
        else rf"Cover Pocket Overlay ($3P_\delta$ Burkert, $r_0={r0_kpc:.2f}$ kpc)"
    )

    # Render Diagnostic Plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 5.6), facecolor='#0b1120')
    for ax in (ax1, ax2):
        ax.set_facecolor('#070c18')
        ax.tick_params(colors='#94a3b8')
        ax.grid(True, alpha=0.2, linestyle=':')

    # Panel 1: Rotation Curve v(r)
    ax1.plot(r_array, v_primary, color='#ef4444', linewidth=2.6, label=lbl_primary, zorder=5)
    ax1.plot(r_array, v_alt, color='#a855f7', linewidth=1.8, linestyle='--', label=lbl_alt, zorder=4)
    if is_cover_pocket:
        ax1.plot(r_array, v_pocket_only, color='#10b981', linewidth=1.6, linestyle='-.',
                 label=rf'Cover Pocket Alone ($\Sigma_p={SIGMA_POCKET_GALAXY_MSUN_PC2:.1f}\,M_\odot/\mathrm{{pc}}^2$)')
    ax1.plot(r_array, v_kep, color='#38bdf8', linewidth=1.6, linestyle=':', label=r'Baryonic Keplerian ($v \propto 1/\sqrt{r}$)')

    if r0_kpc <= r_max:
        ax1.axvline(r0_kpc, color='#22d3ee', linestyle='-.', linewidth=1.5, label=rf'Locked Pocket Core $r_0$ ({r0_kpc:.2f} kpc)')
    if not is_cover_pocket and R_star_kpc is not None:
        ax1.axvline(R_star_kpc, color='#10b981', linestyle='-.', linewidth=1.4, label=f'Transition Radius R* ({R_star_kpc:.1f} kpc)')
    if r_trunc_kpc is not None and r_trunc_kpc <= r_max:
        ax1.axvline(r_trunc_kpc, color='#fbbf24', linestyle='--', linewidth=1.5, label=rf'2D Disk Truncation Cliff ({r_trunc_kpc:.1f} kpc)')

    ax1.set_xlabel('Galactocentric Radius $r$ [kpc]', color='#cbd5e1', fontsize=11)
    ax1.set_ylabel('Orbital Rotation Velocity $v$ [km/s]', color='#cbd5e1', fontsize=11)
    mode_tag = "Cover Pocket Architecture" if is_cover_pocket else "Legacy Underflow Architecture"
    ax1.set_title(f'{target_title}\nGalactic Rotation Curve — {mode_tag}', color='#f8fafc', fontsize=11)
    ax1.legend(loc='lower right', fontsize=8.0, facecolor='#0f172a', edgecolor='none', labelcolor='#e2e8f0')

    # Panel 2: Radial Acceleration Profile & 3 P_delta / a_BTFR / a_Omega Floors
    pos_mask = r_array > 0.0
    r_p = r_array[pos_mask]
    ax2.plot(
        r_p,
        np.where(g_primary[pos_mask] > 0.0, g_primary[pos_mask], np.nan),
        color='#ef4444',
        linewidth=2.5,
        label=r'Primary Total Acceleration $g_{\mathrm{tot}}(r)$',
    )
    ax2.plot(
        r_p,
        np.where(g_alt[pos_mask] > 0.0, g_alt[pos_mask], np.nan),
        color='#a855f7',
        linewidth=1.8,
        linestyle='--',
        label=r'Alternative Model Acceleration $g_{\mathrm{alt}}(r)$',
    )
    ax2.plot(
        r_p,
        np.where(g_bar[pos_mask] > 0.0, g_bar[pos_mask], np.nan),
        color='#38bdf8',
        linewidth=1.6,
        linestyle=':',
        label=r'Baryonic Acceleration $g_{\mathrm{baryon}}(r)$',
    )
    ax2.axhline(G_3PDELTA_MS2, color='#10b981', linestyle='-', linewidth=1.6,
                label=rf'Peak Burkert Ceiling $g_{{3P_\delta}}$ ({G_3PDELTA_MS2:.3e} m/s²)')
    ax2.axhline(A_BTFR_POCKET, color='#f59e0b', linestyle='-.', linewidth=1.6,
                label=rf'Emergent BTFR Scale $a_{{\mathrm{{BTFR}}}}$ ({A_BTFR_POCKET:.3e} m/s²)')
    ax2.axhline(a_omega_used, color='#94a3b8', linestyle='--', linewidth=1.3,
                label=rf'Legacy Floor $a_\Omega$ ({a_omega_used:.3e} m/s²)')

    if r0_kpc <= r_max:
        ax2.axvline(r0_kpc, color='#22d3ee', linestyle='-.', linewidth=1.4, label=rf'Locked Core $r_0$ ({r0_kpc:.2f} kpc)')
    elif not is_cover_pocket and R_star_kpc is not None:
        ax2.axvline(R_star_kpc, color='#10b981', linestyle='-.', linewidth=1.4, label=r'R* Boundary ($g_{\mathrm{eff}} = a_\Omega$)')

    ax2.set_yscale('log')
    ax2.set_xlabel('Galactocentric Radius $r$ [kpc]', color='#cbd5e1', fontsize=11)
    ax2.set_ylabel(r'Radial Acceleration $g$ [m/s²]', color='#cbd5e1', fontsize=11)
    ax2.set_title(r'Radial Acceleration: $3P_\delta$ Burkert vs. $a_\Omega$ Thresholds', color='#f8fafc', fontsize=11)
    ax2.legend(loc='upper right', fontsize=7.8, facecolor='#0f172a', edgecolor='none', labelcolor='#e2e8f0')

    plt.tight_layout()

    df_rotation = pd.DataFrame({
        "radius_kpc": np.round(r_array, 3),
        "v_pocket_kms": np.round(v_pocket, 3),
        "v_pocket_only_kms": np.round(v_pocket_only, 3),
        "v_legacy_aOmega_kms": np.round(v_legacy, 3),
        "v_kepler_kms": np.round(v_kep, 3),
        "v_substrate_kms": np.round(v_primary, 3),
        "g_pocket_tot_ms2": g_pocket_tot,
        "g_legacy_eff_ms2": g_legacy_eff,
        "g_baryon_ms2": g_bar,
        "g_pocket_ms2": g_pocket,
        "a_omega_floor_ms2": a_omega_used,
    })

    idx_focus = int(np.argmin(np.abs(r_array - r_focus)))

    if is_cover_pocket:
        metrics = {
            "True Baryonic Mass": format_baryonic_mass(M_baryon),
            "Locked Pocket Core (r_0)": f"{r0_kpc:.2f} kpc",
            "Pocket Surface Density (Σ_p)": f"{SIGMA_POCKET_GALAXY_MSUN_PC2:.2f} M☉/pc²",
            "Emergent a_BTFR": f"{A_BTFR_POCKET:.3e} m/s²",
            "Velocity @ Focus (r_focus)": f"{v_primary[idx_focus]:.1f} km/s (@ {r_focus:.1f} kpc)",
            "Disk Truncation (3.829 R_d)": f"{r_trunc_kpc:.2f} kpc" if r_trunc_kpc is not None else "N/A",
        }
    else:
        metrics = {
            "True Baryonic Mass": format_baryonic_mass(M_baryon),
            "Disk Scale Length (R_d)": f"{R_d_val:.2f} kpc",
            "Kinematic Underflow Floor (a_Ω)": f"{a_omega_used:.3e} m/s²",
            "Transition Radius (R*)": f"{R_star_kpc:.2f} kpc" if R_star_kpc is not None else "Deep Underflow (< a_Ω)",
            "Flat Asymptotic Velocity": f"{v_flat_legacy:.1f} km/s",
            "Velocity @ Focus (r_focus)": f"{v_primary[idx_focus]:.1f} km/s (@ {r_focus:.1f} kpc)",
            "Disk Truncation (3.829 R_d)": f"{r_trunc_kpc:.2f} kpc" if r_trunc_kpc is not None else "N/A",
        }

    return fig, df_rotation, metrics
