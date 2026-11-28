import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

from .data_loader import (
    C_LIGHT,
    H_GLOBAL,
    K_OMEGA,
    M_PROTON,
    M_SOLAR_KG,
    MPC_TO_METER,
    KPC_TO_METER,
    A_OMEGA,
    RAD_TO_ARCSEC,
    CLUSTER_BIAS_1MB,
    G_EFF,
    ETA_POCKET,
    MSUN_PC2_TO_KGM2,
    SIGMA_POCKET_GALAXY_KGM2,
    SIGMA_POCKET_GALAXY_MSUN_PC2,
    SIGMA_POCKET_CLUSTER_KGM2,
    SIGMA_POCKET_CLUSTER_MSUN_PC2,
    N_POCKET_CLUSTER_BRIDGE,
    SIGMA_0_LEDGER,
    BURKERT_F_B_MAX,
    compute_a_omega,
    get_compiled_node_density,
    get_hlocal,
    compute_locked_pocket_r0_m,
    burkert_2d_proj_factor,
)

# Preset Benchmark Lens Targets
LENS_PRESETS = {
    "Bullet Cluster Subcluster (Bullet)": {
        "M_baryon_solar": 3.8e13,      # 38.0 x 10^12 M_sun (Bullet baryonic budget)
        "scale_radius_kpc": 35.0,
        "z_l": 0.296,
        "desc": "Bullet subcluster passing through Main cluster at z_l = 0.296 (R_E ~ 89 kpc under Sigma_0/3)."
    },
    "Bullet Cluster Main (1E 0657-56)": {
        "M_baryon_solar": 1.1e14,      # 110.0 x 10^12 M_sun
        "scale_radius_kpc": 50.0,
        "z_l": 0.296,
        "desc": "Primary Bullet Cluster halo at z_l = 0.296 under isotropic Sigma_0/3 Cover ceiling."
    },
    "Abell 1689 Core (Cluster Lens)": {
        "M_baryon_solar": 1.5e14,      # 150.0 x 10^12 M_sun
        "scale_radius_kpc": 60.0,
        "z_l": 0.183,
        "desc": "Rich strong-lensing cluster core at z_l = 0.183 showing giant Einstein rings."
    },
    "Milky Way (Total Baryonic Footprint)": {
        "M_baryon_solar": 8.8e10,      # 0.088 x 10^12 M_sun
        "scale_radius_kpc": 8.0,
        "z_l": 0.296,
        "desc": "Total Milky Way baryonic budget in the 3 P_delta galaxy floor regime."
    },
}

_PRESET_ALIASES = {
    "Bullet Cluster Centroid (1E 0657-56)": "Bullet Cluster Subcluster (Bullet)",
}


def _compute_sigma_crit_kgm2(z_l=0.296, z_s=1.0, h_global=H_GLOBAL, omega_m=0.30):
    """
    Fast vectorized angular diameter distance & critical lensing surface density
    Sigma_crit [kg/m^2] plus sky scale [kpc/arcsec] for flat LambdaCDM reference cosmology.
    """
    z_l = float(max(z_l, 1.0e-4))
    z_s = float(max(z_s, z_l + 0.05))
    omega_l = 1.0 - float(omega_m)
    h0_si = (float(h_global) * 1000.0) / MPC_TO_METER
    c_over_h0 = C_LIGHT / h0_si

    z_grid_l = np.linspace(0.0, z_l, 256)
    z_grid_s = np.linspace(0.0, z_s, 256)
    dc_l = c_over_h0 * float(np.trapezoid(1.0 / np.sqrt(omega_m * ((1.0 + z_grid_l) ** 3) + omega_l), z_grid_l))
    dc_s = c_over_h0 * float(np.trapezoid(1.0 / np.sqrt(omega_m * ((1.0 + z_grid_s) ** 3) + omega_l), z_grid_s))

    d_l = dc_l / (1.0 + z_l)
    d_s = dc_s / (1.0 + z_s)
    d_ls = (dc_s - dc_l) / (1.0 + z_s)
    sigma_crit = (C_LIGHT ** 2 / (4.0 * np.pi * G_EFF)) * (d_s / (d_l * d_ls))
    kpc_per_arcsec = (d_l / KPC_TO_METER) * (np.pi / (180.0 * 3600.0))
    return float(sigma_crit), float(kpc_per_arcsec)


def _compute_zero_parameter_lensing_bounds(M_baryon_solar, custom_scale_radius_kpc=None, a_omega=A_OMEGA):
    """
    Computes zero-parameter physical crossover radius b_cross where theta_baryon == theta_scaffold,
    and 1:1 total mass horizon b_1to1_total where M_DM_inferred == M_baryon_total.
    """
    C_scaff = (a_omega * M_PROTON) / (K_OMEGA * M_SOLAR_KG)
    C_scaff_kpc2 = C_scaff * (KPC_TO_METER ** 2)

    b_1to1_total = np.sqrt(M_baryon_solar / C_scaff_kpc2)

    if custom_scale_radius_kpc is None or custom_scale_radius_kpc <= 0:
        scale_radius_kpc = np.clip(8.0 * ((M_baryon_solar / 1.5e11) ** (1.0 / 3.0)), 0.05, 50.0)
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
    lens_distance_mpc=670.0,
):
    """Legacy linear scaffolding lensing profile returning 5-tuple."""
    b_kpc = np.maximum(np.atleast_1d(b_kpc_array), 1e-4)
    b_meters = b_kpc * KPC_TO_METER
    M_baryon_kg = M_baryon_solar * M_SOLAR_KG

    if scale_radius_kpc is not None and scale_radius_kpc > 0:
        a_meters = scale_radius_kpc * KPC_TO_METER
        M_enc_kg = M_baryon_kg * (b_meters ** 2) / ((b_meters + a_meters) ** 2)
    else:
        M_enc_kg = M_baryon_kg

    N_B_enc = M_enc_kg / M_PROTON

    theta_baryon_rad = (4.0 * N_B_enc * K_OMEGA) / ((C_LIGHT ** 2) * b_meters)
    theta_scaffold_rad = (4.0 * a_omega_active * b_meters) / (C_LIGHT ** 2)
    theta_total_rad = theta_baryon_rad + theta_scaffold_rad

    theta_baryon_arcsec = theta_baryon_rad * RAD_TO_ARCSEC
    theta_scaffold_arcsec = theta_scaffold_rad * RAD_TO_ARCSEC
    theta_total_arcsec = theta_total_rad * RAD_TO_ARCSEC

    N_DM = (a_omega_active * (b_meters ** 2)) / K_OMEGA
    M_DM_inferred_solar = (N_DM * M_PROTON) / M_SOLAR_KG

    kpc_per_arcsec = lens_distance_mpc * (np.pi / 648000.0) * 1000.0
    theta_b_arcsec = b_kpc / kpc_per_arcsec

    kappa_baryon = np.gradient(theta_b_arcsec * theta_baryon_arcsec, theta_b_arcsec) / (2.0 * theta_b_arcsec)
    kappa_scaffold = (4.0 * a_omega_active * KPC_TO_METER * kpc_per_arcsec / (C_LIGHT ** 2)) * RAD_TO_ARCSEC
    kappa_total = kappa_baryon + kappa_scaffold

    return theta_total_arcsec, theta_baryon_arcsec, theta_scaffold_arcsec, M_DM_inferred_solar, kappa_total


def _select_ceiling_kgm2(M_baryon_solar, ceiling_mode):
    """Selects hardware pressure ceiling surface density [kg/m^2] & readable label."""
    cm = str(ceiling_mode or "")
    if "3 P" in cm or "3P" in cm or "Galaxy Floor" in cm or "169" in cm:
        if "Auto" not in cm:
            return SIGMA_POCKET_GALAXY_KGM2, SIGMA_POCKET_GALAXY_MSUN_PC2, "Galaxy Floor (3 P_δ)"
    if ("Σ_0 / 3" in cm or "Sigma_0/3" in cm or "Σ_0/3" in cm or "2,600" in cm or "Cluster Ceiling" in cm) and "Auto" not in cm:
        return SIGMA_POCKET_CLUSTER_KGM2, SIGMA_POCKET_CLUSTER_MSUN_PC2, "Isotropic Cover (Σ_0 / 3)"
    if "Sigma_0" in cm and "Auto" not in cm and "/ 3" not in cm and "/3" not in cm:
        g_cap = np.sqrt(8.0 * np.pi * G_EFF * SIGMA_0_LEDGER)
        sig_full = float(g_cap / (BURKERT_F_B_MAX * np.pi * G_EFF))
        return sig_full, sig_full / MSUN_PC2_TO_KGM2, "Full Cover (Σ_0)"

    # Auto mode: switch at M_baryon >= 5.0e12 M_sun
    if float(M_baryon_solar) >= 5.0e12:
        return SIGMA_POCKET_CLUSTER_KGM2, SIGMA_POCKET_CLUSTER_MSUN_PC2, "Auto -> Σ_0 / 3 Cluster"
    return SIGMA_POCKET_GALAXY_KGM2, SIGMA_POCKET_GALAXY_MSUN_PC2, "Auto -> 3 P_δ Galaxy"


def calculate_cover_pocket_lensing_profile(
    M_baryon_solar,
    b_kpc_array,
    a_omega_active=A_OMEGA,
    scale_radius_kpc=None,
    z_l=0.296,
    z_s=1.0,
    h_global=H_GLOBAL,
    ceiling_mode="Auto (3 P_δ Galaxy / Σ_0/3 Cluster)",
):
    """
    Fast vectorized analytic Burkert Cover Pocket 2D weak/strong lensing profile
    using exact burkert_2d_proj_factor(X) + cumulative_trapezoid over 2,000 log points.
    """
    sigma_p_kgm2, sigma_p_msun_pc2, ceiling_label = _select_ceiling_kgm2(M_baryon_solar, ceiling_mode)
    sigma_crit_kgm2, kpc_per_arcsec = _compute_sigma_crit_kgm2(z_l=z_l, z_s=z_s, h_global=h_global)

    M_baryon_kg = float(M_baryon_solar) * M_SOLAR_KG
    r0_m = compute_locked_pocket_r0_m(M_baryon_kg, sigma_p_kgm2)
    r0_kpc = r0_m / KPC_TO_METER

    # Compiled baryonic Hernquist fraction: 15% galaxies for passive-ICM clusters, 100% for galaxies
    is_cluster = float(M_baryon_solar) >= 5.0e12
    f_compiled = 0.15 if is_cluster else 1.0
    m_bar_compiled_kg = f_compiled * M_baryon_kg
    a_kpc = float(scale_radius_kpc) if (scale_radius_kpc is not None and scale_radius_kpc > 0) else np.clip(
        8.0 * ((M_baryon_solar / 1.5e11) ** (1.0 / 3.0)), 0.05, 50.0
    )
    a_m = a_kpc * KPC_TO_METER

    b_kpc = np.maximum(np.asarray(b_kpc_array, dtype=float), 1.0e-3)
    b_max_fine_kpc = max(float(np.max(b_kpc)), 10.0 * r0_kpc, 500.0)
    b_fine_m = np.geomspace(1.0e-3 * KPC_TO_METER, b_max_fine_kpc * KPC_TO_METER, 2000)
    b_fine_kpc = b_fine_m / KPC_TO_METER

    sig2d_p_fine = sigma_p_kgm2 * burkert_2d_proj_factor(b_fine_m / max(r0_m, 1.0))
    m2d_p_fine = (
        np.pi * (b_fine_m[0] ** 2) * sig2d_p_fine[0]
        + cumulative_trapezoid(2.0 * np.pi * b_fine_m * sig2d_p_fine, b_fine_m, initial=0.0)
    )

    m2d_bar_fine = m_bar_compiled_kg * (b_fine_m ** 2) / ((b_fine_m + a_m) ** 2)
    sig2d_bar_fine = (m_bar_compiled_kg / np.pi) * a_m / ((b_fine_m + a_m) ** 3)

    sig2d_tot_fine = sig2d_p_fine + sig2d_bar_fine
    m2d_tot_fine = m2d_p_fine + m2d_bar_fine
    kappa_loc_fine = sig2d_tot_fine / sigma_crit_kgm2
    kappa_mean_fine = m2d_tot_fine / (np.pi * (b_fine_m ** 2) * sigma_crit_kgm2)

    # Strong-lensing Einstein radius R_E where mean kappa(<R) = 1.0
    if float(np.max(kappa_mean_fine)) < 1.0:
        r_einstein_kpc = 0.0
    else:
        above = np.flatnonzero(kappa_mean_fine >= 1.0)
        last_i = int(above[-1])
        if last_i + 1 < len(b_fine_kpc) and kappa_mean_fine[last_i] >= 1.0 > kappa_mean_fine[last_i + 1]:
            frac = (kappa_mean_fine[last_i] - 1.0) / max(kappa_mean_fine[last_i] - kappa_mean_fine[last_i + 1], 1.0e-12)
            r_einstein_kpc = float(b_fine_kpc[last_i] + frac * (b_fine_kpc[last_i + 1] - b_fine_kpc[last_i]))
        else:
            r_einstein_kpc = float(b_fine_kpc[last_i])

    m2d_250_solar = float(np.interp(250.0, b_fine_kpc, m2d_tot_fine) / M_SOLAR_KG)
    kappa_peak_0 = float((sigma_p_kgm2 * (np.pi / 2.0) + (m_bar_compiled_kg / (np.pi * (a_m ** 2)))) / sigma_crit_kgm2)

    # Interpolate profiles onto requested b_kpc grid
    b_m = b_kpc * KPC_TO_METER
    sig2d_p = sigma_p_kgm2 * burkert_2d_proj_factor(b_m / max(r0_m, 1.0))
    sig2d_bar = (m_bar_compiled_kg / np.pi) * a_m / ((b_m + a_m) ** 3)
    m2d_p = np.interp(b_kpc, b_fine_kpc, m2d_p_fine)
    m2d_bar = m_bar_compiled_kg * (b_m ** 2) / ((b_m + a_m) ** 2)
    m2d_tot = m2d_p + m2d_bar

    theta_pocket_arcsec = ((4.0 * G_EFF * m2d_p) / ((C_LIGHT ** 2) * b_m)) * RAD_TO_ARCSEC
    theta_baryon_arcsec = ((4.0 * G_EFF * m2d_bar) / ((C_LIGHT ** 2) * b_m)) * RAD_TO_ARCSEC
    theta_total_arcsec = theta_pocket_arcsec + theta_baryon_arcsec
    theta_legacy_scaffold_arcsec = ((4.0 * a_omega_active * b_m) / (C_LIGHT ** 2)) * RAD_TO_ARCSEC

    kappa_local = (sig2d_p + sig2d_bar) / sigma_crit_kgm2
    kappa_mean_enclosed = m2d_tot / (np.pi * (b_m ** 2) * sigma_crit_kgm2)

    return {
        "b_kpc": b_kpc,
        "theta_total_arcsec": theta_total_arcsec,
        "theta_baryon_arcsec": theta_baryon_arcsec,
        "theta_pocket_arcsec": theta_pocket_arcsec,
        "theta_legacy_scaffold_arcsec": theta_legacy_scaffold_arcsec,
        "kappa_local": kappa_local,
        "kappa_mean_enclosed": kappa_mean_enclosed,
        "M_2D_pocket_solar": m2d_p / M_SOLAR_KG,
        "M_2D_total_solar": m2d_tot / M_SOLAR_KG,
        "r0_kpc": r0_kpc,
        "sigma_p_msun_pc2": sigma_p_msun_pc2,
        "ceiling_label": ceiling_label,
        "sigma_crit_msun_pc2": sigma_crit_kgm2 / MSUN_PC2_TO_KGM2,
        "kpc_per_arcsec": kpc_per_arcsec,
        "R_E_kpc": r_einstein_kpc,
        "M_2D_250kpc_solar": m2d_250_solar,
        "kappa_peak_0": kappa_peak_0,
        "f_compiled": f_compiled,
    }


def run_lensing_analysis(
    lens_target_name="Bullet Cluster Subcluster (Bullet)",
    custom_M_baryon=None,
    custom_scale_radius_kpc=None,
    h_global=H_GLOBAL,
    use_log_scale=False,
    use_catalog_query=False,
    target_pos_mpc=None,
    sigma_mpc=1.8,
    tree=None,
    gal_positions=None,
    nodes_per_galaxy=None,
    lensing_model="Cover Pocket (Linear K_Ω + Burkert)",
    ceiling_mode="Auto (3 P_δ Galaxy / Σ_0/3 Cluster)",
    z_l=None,
):
    """
    Main Streamlit UI Execution Engine for Gravitational Lensing.
    Supports both the unified Two-Ceiling Cover Pocket mode (3 P_delta / Sigma_0/3)
    with fast analytic Burkert 2D projection & Einstein radius extraction, and
    the Legacy Linear Scaffold (a_Omega) mode.
    """
    resolved_target = _PRESET_ALIASES.get(lens_target_name, lens_target_name)
    preset_z_l = 0.296

    if use_catalog_query and target_pos_mpc is not None and tree is not None:
        _ = get_compiled_node_density(tree, gal_positions, nodes_per_galaxy, target_pos_mpc, sigma_mpc=sigma_mpc)
        h_local = get_hlocal(tree, gal_positions, nodes_per_galaxy, target_pos_mpc, sigma_mpc=sigma_mpc, h_global=h_global)
        a_omega_used = compute_a_omega(h_local)
        M_baryon = custom_M_baryon if custom_M_baryon is not None else 1.0e11
        target_title = f"{resolved_target}\n(X={target_pos_mpc[0]:.1f}, Y={target_pos_mpc[1]:.1f}, Z={target_pos_mpc[2]:.1f} Mpc)"
    else:
        if resolved_target in LENS_PRESETS:
            preset = LENS_PRESETS[resolved_target]
            M_baryon = custom_M_baryon if custom_M_baryon is not None else preset["M_baryon_solar"]
            custom_scale_radius_kpc = custom_scale_radius_kpc if custom_scale_radius_kpc is not None else preset.get("scale_radius_kpc")
            preset_z_l = preset.get("z_l", 0.296)
        else:
            M_baryon = custom_M_baryon if custom_M_baryon is not None else 3.8e13

        a_omega_used = compute_a_omega(h_global)
        h_local = h_global
        target_title = resolved_target

    active_z_l = float(z_l) if z_l is not None else float(preset_z_l)
    is_cover_pocket = bool(lensing_model and "Cover Pocket" in str(lensing_model))

    if is_cover_pocket:
        sigma_p_kgm2, _, _ = _select_ceiling_kgm2(M_baryon, ceiling_mode)
        r0_kpc = compute_locked_pocket_r0_m(M_baryon * M_SOLAR_KG, sigma_p_kgm2) / KPC_TO_METER
        b_min = max(0.5, 0.02 * r0_kpc)
        b_max = max(350.0, 5.5 * r0_kpc)
        b_array = np.geomspace(b_min, b_max, 300) if use_log_scale else np.linspace(b_min, b_max, 300)

        lens = calculate_cover_pocket_lensing_profile(
            M_baryon_solar=M_baryon,
            b_kpc_array=b_array,
            a_omega_active=a_omega_used,
            scale_radius_kpc=custom_scale_radius_kpc,
            z_l=active_z_l,
            h_global=h_global,
            ceiling_mode=ceiling_mode,
        )

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 5.6), facecolor='#0b1120')
        for ax in (ax1, ax2):
            ax.set_facecolor('#070c18')
            ax.tick_params(colors='#94a3b8')
            ax.grid(True, alpha=0.2, linestyle=':')

        # Panel 1: Optical Deflection Angle theta(b)
        plot1 = ax1.loglog if use_log_scale else ax1.plot
        plot1(b_array, lens["theta_total_arcsec"], color='#ef4444', linewidth=2.5,
              label=r'Total Cover Pocket Deflection $\theta_{\mathrm{tot}}(b)$')
        plot1(b_array, lens["theta_pocket_arcsec"], color='#10b981', linewidth=2.0,
              label=rf"Burkert Cover Pocket ($\Sigma_p={lens['sigma_p_msun_pc2']:.1f}\,M_\odot/\mathrm{{pc}}^2$)")
        plot1(b_array, lens["theta_baryon_arcsec"], color='#38bdf8', linewidth=1.7, linestyle='--',
              label=rf"Compiled Baryons ($f_{{\mathrm{{comp}}}}={lens['f_compiled']:.2f}$)")
        plot1(b_array, lens["theta_legacy_scaffold_arcsec"], color='#a855f7', linewidth=1.6, linestyle=':',
              label=r'Legacy Scaffolding Tax $\theta_{\mathrm{scaffold}} \propto b$')

        ax1.axvline(lens["r0_kpc"], color='#22d3ee', linestyle='-.', linewidth=1.4,
                    label=rf"Locked Pocket Core $r_0={lens['r0_kpc']:.1f}$ kpc")
        if lens["R_E_kpc"] > 0.0:
            ax1.axvline(lens["R_E_kpc"], color='#facc15', linestyle='--', linewidth=1.5,
                        label=rf"Einstein Radius $R_E={lens['R_E_kpc']:.1f}$ kpc")

        ax1.set_xlabel('Impact Parameter $b$ [kpc]', color='#cbd5e1', fontsize=11)
        ax1.set_ylabel(r'Optical Deflection Angle $\theta$ [arcsec]', color='#cbd5e1', fontsize=11)
        ax1.set_title(f'{target_title}\nCover Pocket Optical Deflection ({lens["ceiling_label"]})', color='#f8fafc', fontsize=11)
        ax1.legend(loc='upper right', fontsize=8.0, facecolor='#0f172a', edgecolor='none', labelcolor='#e2e8f0')

        # Panel 2: 2D Convergence kappa(b) & Enclosed Mean kappa(<b)
        plot2 = ax2.loglog if use_log_scale else ax2.plot
        plot2(b_array, lens["kappa_local"], color='#38bdf8', linewidth=2.1,
              label=r'Local Convergence $\kappa(b) = \Sigma_{2D}/\Sigma_{\mathrm{crit}}$')
        plot2(b_array, lens["kappa_mean_enclosed"], color='#f59e0b', linewidth=2.5,
              label=r'Enclosed Mean $\bar{\kappa}(<b) = M_{2D}(<b)/(\pi b^2 \Sigma_{\mathrm{crit}})$')
        ax2.axhline(1.0, color='#ffffff', linestyle='--', linewidth=1.4,
                    label=rf'Strong-Lensing Critical Line $\kappa = 1.0$ ($\Sigma_{{\mathrm{{crit}}}}={lens["sigma_crit_msun_pc2"]:.0f}\,M_\odot/\mathrm{{pc}}^2$)')
        if lens["R_E_kpc"] > 0.0:
            ax2.axvline(lens["R_E_kpc"], color='#22c55e', linestyle='-', linewidth=1.8,
                        label=rf'Einstein Radius $R_E = {lens["R_E_kpc"]:.1f}$ kpc')

        badge_2 = (
            rf"$M_{{2D}}(<250\ \mathrm{{kpc}}) = {lens['M_2D_250kpc_solar']:.2e}\ M_\odot$" "\n"
            rf"Cluster SZ-Lensing Bias $(1-b) = {CLUSTER_BIAS_1MB:.4f}$" "\n"
            rf"Cluster/Galaxy Pressure Bridge $1/(9\mathcal{{C}}_\delta) = {N_POCKET_CLUSTER_BRIDGE:.1f}$"
        )
        ax2.text(
            0.97, 0.96, badge_2, transform=ax2.transAxes, va='top', ha='right', fontsize=8.2, color='#f8fafc',
            bbox=dict(boxstyle='round,pad=0.45', facecolor='#0f172a', edgecolor='#38bdf8', alpha=0.9),
        )
        ax2.set_xlabel('Impact Parameter $b$ [kpc]', color='#cbd5e1', fontsize=11)
        ax2.set_ylabel(r'Dimensionless Convergence $\kappa$', color='#cbd5e1', fontsize=11)
        ax2.set_title(r'2D Weak & Strong Lensing Convergence $\kappa(b)$ ($z_s=1.0$)', color='#f8fafc', fontsize=11)
        ax2.legend(loc='lower left', fontsize=8.0, facecolor='#0f172a', edgecolor='none', labelcolor='#e2e8f0')

        plt.tight_layout()

        df_lensing = pd.DataFrame({
            "impact_parameter_kpc": np.round(b_array, 3),
            "theta_total_arcsec": lens["theta_total_arcsec"],
            "theta_baryon_arcsec": lens["theta_baryon_arcsec"],
            "theta_pocket_arcsec": lens["theta_pocket_arcsec"],
            "theta_legacy_scaffold_arcsec": lens["theta_legacy_scaffold_arcsec"],
            "kappa_local": lens["kappa_local"],
            "kappa_mean_enclosed": lens["kappa_mean_enclosed"],
            "M_2D_pocket_solar": lens["M_2D_pocket_solar"],
            "M_2D_total_solar": lens["M_2D_total_solar"],
        })

        re_str = f"{lens['R_E_kpc']:.1f} kpc" if lens["R_E_kpc"] > 0.0 else "Sub-critical (< 1.0)"
        metrics = {
            "Hardware Pressure Ceiling": f"{lens['ceiling_label']} ({lens['sigma_p_msun_pc2']:,.2f} M☉/pc²)",
            "Locked Core Radius (r_0)": f"{lens['r0_kpc']:.2f} kpc",
            "Peak Convergence κ(0)": f"{lens['kappa_peak_0']:.3f}",
            "Einstein Radius (R_E)": re_str,
            "M_2D (< 250 kpc)": f"{lens['M_2D_250kpc_solar']:.2e} M☉",
            "SZ-Lensing Cluster Bias (1-b)": f"{CLUSTER_BIAS_1MB:.4f} (-{(1.0 - CLUSTER_BIAS_1MB) * 100:.2f}%)",
        }
        return fig, df_lensing, metrics

    # ==========================================
    # LEGACY LINEAR SCAFFOLD (a_Omega) MODE
    # ==========================================
    b_min, b_max, b_focus, b_1to1_total, scale_radius = _compute_zero_parameter_lensing_bounds(
        M_baryon_solar=M_baryon,
        custom_scale_radius_kpc=custom_scale_radius_kpc,
        a_omega=a_omega_used,
    )

    b_array = np.geomspace(b_min, b_max, 300) if use_log_scale else np.linspace(b_min, b_max, 300)
    theta_tot, theta_bar, theta_scaff, M_DM_inf, kappa_tot = calculate_lensing_profile(
        M_baryon_solar=M_baryon,
        b_kpc_array=b_array,
        a_omega_active=a_omega_used,
        scale_radius_kpc=scale_radius,
    )

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

    plot_fn = ax1.loglog if use_log_scale else ax1.plot
    plot_fn(b_array, theta_tot, color='#ef4444', linewidth=2.5, label=r'Total Substrate Deflection $\theta_{\mathrm{total}}(b)$')
    plot_fn(b_array, theta_bar, color='#38bdf8', linewidth=1.8, linestyle='--', label=r'Baryonic Component $\theta_{\mathrm{baryon}}$')
    plot_fn(b_array, theta_scaff, color='#a855f7', linewidth=1.8, linestyle=':', label=r'Scaffolding Tax $\theta_{\mathrm{scaffold}} \propto b$')
    ax1.axvline(b_focus, color='#10b981', linestyle=':', label=f'Crossover Focus ($b={b_focus:.2f}$ kpc)')
    ax1.set_xlabel('Impact Parameter $b$ [kpc]', color='#cbd5e1', fontsize=11)
    ax1.set_ylabel(r'Optical Deflection Angle $\theta$ [arcsec]', color='#cbd5e1', fontsize=11)
    ax1.set_title(f'{target_title}\nLight Deflection Profile (Legacy $a_\\Omega$)', color='#f8fafc', fontsize=11)
    ax1.legend(loc='upper right', facecolor='#0f172a', edgecolor='none', labelcolor='#e2e8f0')

    plot_fn2 = ax2.loglog if use_log_scale else ax2.plot
    plot_fn2(b_array, M_DM_inf / mass_factor, color='#f59e0b', linewidth=2.5, label=r'Inferred Halo Mass $M_{\mathrm{DM}}(b) = \frac{a_\Omega b^2}{\mathcal{K}_\Omega}$')
    baryon_label_str = f"True Baryonic Baseline ({_format_mass_label_matplotlib(M_baryon)})"
    ax2.axhline(M_baryon / mass_factor, color='#38bdf8', linestyle='--', label=baryon_label_str)
    if M_baryon >= 1e12:
        ax2.axhline(
            (CLUSTER_BIAS_1MB * M_baryon) / mass_factor,
            color='#ec4899',
            linestyle='-.',
            linewidth=1.6,
            label=rf'WL Hydrostatic Baseline $(1-b)M_B$ ({CLUSTER_BIAS_1MB:.4f})',
        )
    ax2.axvline(b_1to1_total, color='#38bdf8', linestyle=':', label=f'1:1 Total Mass Horizon ($b={b_1to1_total:.2f}$ kpc)')
    ax2.set_xlabel('Impact Parameter $b$ [kpc]', color='#cbd5e1', fontsize=11)
    ax2.set_ylabel(f'Inferred Halo Mass [${mass_unit}$]', color='#cbd5e1', fontsize=11)
    ax2.set_title(r'Dark Matter Halo Illusion: $M_{\mathrm{DM}} \propto b^2$', color='#f8fafc', fontsize=11)
    ax2.legend(loc='upper left', facecolor='#0f172a', edgecolor='none', labelcolor='#e2e8f0')

    plt.tight_layout()

    df_lensing = pd.DataFrame({
        "impact_parameter_kpc": np.round(b_array, 3),
        "theta_total_arcsec": [f"{val:.4e}" if val < 1e-3 else f"{val:.4f}" for val in theta_tot],
        "theta_baryon_arcsec": [f"{val:.4e}" if val < 1e-3 else f"{val:.4f}" for val in theta_bar],
        "theta_scaffold_arcsec": [f"{val:.4e}" if val < 1e-3 else f"{val:.4f}" for val in theta_scaff],
        "convergence_kappa": [f"{val:.6e}" for val in kappa_tot],
        "inferred_M_DM_solar": [f"{val:.2e}" for val in M_DM_inf],
        "dark_to_light_ratio": np.round(M_DM_inf / M_baryon, 2),
    })

    idx_focus = int(np.argmin(np.abs(b_array - b_focus)))
    metrics = {
        "Target Baseline Acceleration (a_Ω)": f"{a_omega_used:.3e} m/s²",
        "Local Expansion Rate (H_local)": f"{h_local:.2f} km/s/Mpc",
        "1:1 Enclosed Mass Radius (b_cross)": f"{b_focus:.2f} kpc",
        "1:1 Total Mass Horizon (b_1:1)": f"{b_1to1_total:.2f} kpc",
        "Deflection @ Crossover": f"{theta_tot[idx_focus]:.4e} arcsec" if theta_tot[idx_focus] < 1e-2 else f"{theta_tot[idx_focus]:.2f} arcsec",
        "Total Dark/Light @ b_cross": f"{(M_DM_inf[idx_focus] / M_baryon):.2f} : 1",
        "SZ-Lensing Cluster Bias (1-b)": f"{CLUSTER_BIAS_1MB:.4f} (-{(1.0 - CLUSTER_BIAS_1MB)*100:.2f}%)",
    }

    return fig, df_lensing, metrics
