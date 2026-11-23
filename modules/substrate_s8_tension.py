import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

from .data_loader import (
    M_PROTON,
    M_SOLAR_KG,
    KPC_TO_METER,
    AU_TO_METER,
    E_RYD_J,
    T_LOCK_RYD_K,
    OMEGA_M_VOID,
    OMEGA_L_VOID,
    OMEGA_M_GLOBAL,
    OMEGA_M_CORE,
    SIGMA_8_VOID,
    S8_VOID,
    S8_GLOBAL,
    S8_CORE,
    CLUSTER_BIAS_1MB,
    RHO_GC_EXACT,
    RHO_GC_THRESHOLD,
    get_compiled_node_density_batch,
    evaluate_two_condition_compilation_gate,
)

# Observational 2020-2026 Weak Lensing, CMB, and X-Ray Survey Benchmarks: (S8_obs, 1-sigma_err)
OBSERVATIONAL_S8_SURVEYS = {
    "KiDS-1000": (0.766, 0.020),
    "Subaru HSC Y3": (0.769, 0.034),
    "DES Y6 3x2pt": (0.789, 0.012),
    "KiDS Legacy": (0.815, 0.018),
    "Planck CMB": (0.832, 0.013),
    "eROSITA X-Ray": (0.860, 0.016),
}


def run_s8_tension_analysis(
    tree, gal_positions, nodes_per_galaxy, gal_names,
    sigma_halo_mpc=0.05,
    activation_mode="Smooth Volume Fraction (Coarse-Grained)",
    slice_extent_mpc=15.0,
    grid_n=80
):
    """
    Evaluates the environmental S_8 weak-lensing tension field S_8(r) driven by
    3D compiled node density rho_N(r) across catalog galaxies and an equatorial plane slice.
    Returns (fig, df_out, metrics_dict).
    """
    n_total = len(gal_positions)
    if n_total > 10000:
        rng = np.random.default_rng(42)
        eval_idx = np.sort(rng.choice(n_total, size=10000, replace=False))
    else:
        eval_idx = np.arange(n_total)

    pos_eval = np.asarray(gal_positions[eval_idx], dtype=float)
    nodes_eval = np.asarray(nodes_per_galaxy[eval_idx], dtype=float)
    names_eval = np.asarray(gal_names[eval_idx], dtype=str)
    dists_eval = np.linalg.norm(pos_eval, axis=1)
    mass_solar_eval = (nodes_eval * M_PROTON) / M_SOLAR_KG

    # 1. Vectorized catalog density evaluation (zero Python loops)
    rho_N = get_compiled_node_density_batch(
        tree, gal_positions, nodes_per_galaxy, pos_eval, sigma_mpc=float(sigma_halo_mpc)
    )

    is_strict = "Strict" in str(activation_mode)
    if is_strict:
        theta_act = (rho_N >= RHO_GC_THRESHOLD).astype(float)
    else:
        theta_act = rho_N / (rho_N + RHO_GC_EXACT)

    omega_m_local = OMEGA_M_VOID + (OMEGA_M_GLOBAL - OMEGA_M_VOID) * theta_act
    s8_local = SIGMA_8_VOID * np.sqrt(omega_m_local / 0.300)

    # 2. Vectorized 2D equatorial slice (Z = 0)
    xs = np.linspace(-slice_extent_mpc, slice_extent_mpc, grid_n)
    ys = np.linspace(-slice_extent_mpc, slice_extent_mpc, grid_n)
    XX, YY = np.meshgrid(xs, ys)
    slice_pts = np.column_stack([XX.ravel(), YY.ravel(), np.zeros(XX.size, dtype=float)])

    dx = (2.0 * slice_extent_mpc) / max(1, grid_n - 1)
    sigma_slice = max(float(sigma_halo_mpc), min(0.40, dx))
    rho_slice = get_compiled_node_density_batch(
        tree, gal_positions, nodes_per_galaxy, slice_pts, sigma_mpc=sigma_slice
    )
    if is_strict:
        theta_slice = (rho_slice >= RHO_GC_THRESHOLD).astype(float)
    else:
        theta_slice = rho_slice / (rho_slice + RHO_GC_EXACT)

    omega_m_slice = OMEGA_M_VOID + (OMEGA_M_GLOBAL - OMEGA_M_VOID) * theta_slice
    s8_slice = SIGMA_8_VOID * np.sqrt(omega_m_slice / 0.300)
    s8_slice_grid = s8_slice.reshape(grid_n, grid_n)
    rho_slice_grid = rho_slice.reshape(grid_n, grid_n)

    # 3. Radial distance shells & summary statistics
    max_dist = float(np.max(dists_eval)) if len(dists_eval) > 0 else float(slice_extent_mpc)
    if max_dist <= 35.0:
        shell_edges = np.array([0.0, 2.0, 4.0, 6.0, 8.0, 12.0])
    else:
        shell_edges = np.linspace(0.0, max(1.0, min(slice_extent_mpc, max_dist)), 6)

    shell_lines = []
    for k in range(len(shell_edges) - 1):
        in_shell = (dists_eval >= shell_edges[k]) & (dists_eval < shell_edges[k + 1])
        cnt = int(np.sum(in_shell))
        if cnt > 0:
            shell_lines.append(
                f"{shell_edges[k]:.0f}-{shell_edges[k+1]:.0f} Mpc: {float(np.mean(s8_local[in_shell])):.4f} (n={cnt})"
            )

    s8_mean = float(np.mean(s8_local))
    s8_median = float(np.median(s8_local))
    s8_min = float(np.min(s8_local))
    s8_max = float(np.max(s8_local))
    s8_std = float(np.std(s8_local))

    # 4. 4-Panel Dark-Theme Diagnostic Figure
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 11), facecolor='#0b1120')
    for ax in (ax1, ax2, ax3, ax4):
        ax.set_facecolor('#070c18')
        ax.tick_params(colors='#94a3b8')
        ax.grid(True, alpha=0.2, linestyle=':')

    # --- Panel 1: 2D Equatorial Slice S_8(X, Y, Z=0) ---
    im1 = ax1.imshow(
        s8_slice_grid,
        origin='lower',
        extent=[-slice_extent_mpc, slice_extent_mpc, -slice_extent_mpc, slice_extent_mpc],
        cmap='magma',
        vmin=S8_VOID,
        vmax=S8_GLOBAL,
        aspect='equal',
    )
    cbar1 = fig.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
    cbar1.ax.tick_params(colors='#94a3b8')
    cbar1.set_label(r'Local Growth Parameter $S_8(\vec{r})$', color='#cbd5e1', fontsize=10)

    if np.max(rho_slice_grid) >= RHO_GC_THRESHOLD and np.min(rho_slice_grid) < RHO_GC_THRESHOLD:
        ax1.contour(XX, YY, rho_slice_grid, levels=[RHO_GC_THRESHOLD], colors=['#22d3ee'], linewidths=1.3)

    slab_mask = (
        (np.abs(pos_eval[:, 2]) <= 0.15 * slice_extent_mpc)
        & (np.abs(pos_eval[:, 0]) <= slice_extent_mpc)
        & (np.abs(pos_eval[:, 1]) <= slice_extent_mpc)
    )
    if np.any(slab_mask):
        slab_pos = pos_eval[slab_mask]
        slab_s8 = s8_local[slab_mask]
        slab_act = rho_N[slab_mask] >= RHO_GC_THRESHOLD
        if np.any(~slab_act):
            ax1.scatter(
                slab_pos[~slab_act, 0], slab_pos[~slab_act, 1],
                c=slab_s8[~slab_act], cmap='magma', vmin=S8_VOID, vmax=S8_GLOBAL,
                s=18, edgecolors='#94a3b8', linewidths=0.5, alpha=0.85,
                label=r'Void Floor ($\rho_N < 4.0$)',
            )
        if np.any(slab_act):
            ax1.scatter(
                slab_pos[slab_act, 0], slab_pos[slab_act, 1],
                c=slab_s8[slab_act], cmap='magma', vmin=S8_VOID, vmax=S8_GLOBAL,
                s=32, edgecolors='#22d3ee', linewidths=0.9, alpha=0.95,
                label=r'Compiled Core ($\rho_N \geq 4.0$)',
            )

    ax1.set_xlabel(r'Supergalactic / ICRS $X$ [Mpc]', color='#cbd5e1', fontsize=10.5)
    ax1.set_ylabel(r'Supergalactic / ICRS $Y$ [Mpc]', color='#cbd5e1', fontsize=10.5)
    ax1.set_title(rf'2D Equatorial Slice $S_8(X, Y, Z=0)$ ($\sigma_{{\mathrm{{slice}}}}={sigma_slice:.2f}$ Mpc)', color='#f8fafc', fontsize=11.5)
    if np.any(slab_mask):
        ax1.legend(loc='upper right', facecolor='#0f172a', edgecolor='none', labelcolor='#e2e8f0', fontsize=8.5)

    # --- Panel 2: Catalog S_8 Distribution vs. 2026 Observational Survey Landscape ---
    ax2.hist(s8_local, bins=35, range=(0.745, 0.870), density=True, color='#38bdf8', alpha=0.42, edgecolor='#0ea5e9', label=r'Catalog $S_8$ PDF')
    if s8_std > 1e-6 and len(s8_local) > 5:
        try:
            kde = gaussian_kde(s8_local)
            x_kde = np.linspace(0.745, 0.870, 300)
            ax2.plot(x_kde, kde(x_kde), color='#f8fafc', linewidth=2.0, label='Catalog KDE')
        except Exception:
            pass

    ref_lines = [
        (S8_VOID, '#22d3ee', '--', rf'$S_{{8,\mathrm{{void}}}}={S8_VOID:.4f}$'),
        (0.766, '#a78bfa', ':', r'KiDS-1000 ($0.766$)'),
        (0.769, '#34d399', ':', r'HSC Y3 ($0.769$)'),
        (0.789, '#60a5fa', '-.', r'DES Y6 ($0.789$)'),
        (0.815, '#fbbf24', '-.', r'KiDS Legacy ($0.815$)'),
        (S8_GLOBAL, '#f59e0b', '--', rf'Planck / $S_{{8,\mathrm{{global}}}}={S8_GLOBAL:.4f}$'),
        (S8_CORE, '#f43f5e', '--', rf'eROSITA / $S_{{8,\mathrm{{core}}}}={S8_CORE:.4f}$'),
        (s8_mean, '#ffffff', '-', rf'Mean $\langle S_8 \rangle={s8_mean:.4f}$'),
    ]
    for x_ref, col_ref, ls_ref, lbl_ref in ref_lines:
        ax2.axvline(x_ref, color=col_ref, linestyle=ls_ref, linewidth=1.5, label=lbl_ref)

    shell_text = "\n".join(shell_lines[:3]) if shell_lines else "All radii aggregated"
    stats_badge = (
        f"Mean: {s8_mean:.4f} | Med: {s8_median:.4f}\n"
        f"Std: {s8_std:.4f} | Range: {s8_min:.4f}-{s8_max:.4f}\n"
        f"{shell_text}"
    )
    ax2.text(
        0.03, 0.96, stats_badge, transform=ax2.transAxes,
        ha='left', va='top', fontsize=7.8, color='#e2e8f0',
        bbox=dict(boxstyle='round,pad=0.45', facecolor='#0f172a', edgecolor='#38bdf8', alpha=0.88),
    )

    ax2.set_xlim(0.742, 0.875)
    ax2.set_xlabel(r'Structure Growth Parameter $S_8$', color='#cbd5e1', fontsize=10.5)
    ax2.set_ylabel('Probability Density', color='#cbd5e1', fontsize=10.5)
    ax2.set_title(r'Catalog $S_8$ Distribution vs. Observational Surveys', color='#f8fafc', fontsize=11.5)
    ax2.legend(loc='upper right', ncol=2, facecolor='#0f172a', edgecolor='none', labelcolor='#e2e8f0', fontsize=7.3)

    # --- Panel 3: Environmental Scaling S_8 vs Compiled Baryonic Mass ---
    log_mass = np.log10(np.maximum(mass_solar_eval, 1.0))
    sc3 = ax3.scatter(log_mass, s8_local, c=dists_eval, cmap='viridis', s=14, alpha=0.72, edgecolors='none')
    cbar3 = fig.colorbar(sc3, ax=ax3, fraction=0.046, pad=0.04)
    cbar3.ax.tick_params(colors='#94a3b8')
    cbar3.set_label(r'Radial Distance $d$ [Mpc]', color='#cbd5e1', fontsize=10)

    ax3.axhline(S8_VOID, color='#22d3ee', linestyle='--', linewidth=1.6, label=rf'Void Floor $S_{{8,\mathrm{{void}}}}={S8_VOID:.4f}$')
    ax3.axhline(S8_GLOBAL, color='#f59e0b', linestyle='--', linewidth=1.6, label=rf'Global Cover $S_{{8,\mathrm{{global}}}}={S8_GLOBAL:.4f}$')
    ax3.axhline(s8_mean, color='#e2e8f0', linestyle=':', linewidth=1.3, label=rf'Catalog Mean $\langle S_8 \rangle={s8_mean:.4f}$')

    ax3.set_ylim(0.745, 0.868)
    ax3.set_xlabel(r'Compiled Baryonic Mass $\log_{10}(M_B / M_\odot)$', color='#cbd5e1', fontsize=10.5)
    ax3.set_ylabel(r'Local Structure Growth $S_8$', color='#cbd5e1', fontsize=10.5)
    ax3.set_title(r'Environmental Compilation Filter ($\sigma_{\mathrm{halo}}=' + f'{sigma_halo_mpc:.2f}$ Mpc)', color='#f8fafc', fontsize=11.5)
    ax3.legend(loc='upper left', facecolor='#0f172a', edgecolor='none', labelcolor='#e2e8f0', fontsize=8.5)

    # --- Panel 4: 4-Tier Hardware Allocation & Hydrodynamic Cluster Mass Bias ---
    labels_p4 = [
        "Void Floor (Substrate)",
        "KiDS-1000 (Weak Lensing)",
        "Subaru HSC Y3",
        "Catalog Mean <S8>",
        "DES Y6 3x2pt",
        "KiDS Legacy",
        "Global Cover (Substrate)",
        "Planck CMB",
        "eROSITA X-Ray",
        "Cluster Core (Substrate)",
    ]
    vals_p4 = [
        S8_VOID,
        OBSERVATIONAL_S8_SURVEYS["KiDS-1000"][0],
        OBSERVATIONAL_S8_SURVEYS["Subaru HSC Y3"][0],
        s8_mean,
        OBSERVATIONAL_S8_SURVEYS["DES Y6 3x2pt"][0],
        OBSERVATIONAL_S8_SURVEYS["KiDS Legacy"][0],
        S8_GLOBAL,
        OBSERVATIONAL_S8_SURVEYS["Planck CMB"][0],
        OBSERVATIONAL_S8_SURVEYS["eROSITA X-Ray"][0],
        S8_CORE,
    ]
    errs_p4 = [
        0.0,
        OBSERVATIONAL_S8_SURVEYS["KiDS-1000"][1],
        OBSERVATIONAL_S8_SURVEYS["Subaru HSC Y3"][1],
        0.0,
        OBSERVATIONAL_S8_SURVEYS["DES Y6 3x2pt"][1],
        OBSERVATIONAL_S8_SURVEYS["KiDS Legacy"][1],
        0.0,
        OBSERVATIONAL_S8_SURVEYS["Planck CMB"][1],
        OBSERVATIONAL_S8_SURVEYS["eROSITA X-Ray"][1],
        0.0,
    ]
    colors_p4 = [
        '#22d3ee', '#94a3b8', '#94a3b8', '#38bdf8', '#94a3b8',
        '#94a3b8', '#f59e0b', '#cbd5e1', '#cbd5e1', '#f43f5e'
    ]

    y_pos = np.arange(len(labels_p4))
    for idx_row, (y_val, x_val, err_val, col_val, lbl_val) in enumerate(zip(y_pos, vals_p4, errs_p4, colors_p4, labels_p4)):
        if "Substrate" in lbl_val or "Catalog Mean" in lbl_val:
            ax4.barh(y_val, x_val - 0.70, left=0.70, height=0.55, color=col_val, alpha=0.80, edgecolor='#f8fafc', linewidth=0.6)
            ax4.text(x_val + 0.003, y_val, f"{x_val:.4f}", va='center', ha='left', color=col_val, fontsize=8.2, fontweight='bold')
        else:
            ax4.errorbar(x_val, y_val, xerr=err_val, fmt='o', color=col_val, ecolor=col_val, elinewidth=1.8, capsize=4, markersize=6)
            ax4.text(x_val + err_val + 0.003, y_val, f"{x_val:.3f}±{err_val:.3f}", va='center', ha='left', color='#cbd5e1', fontsize=7.8)

    ax4.set_yticks(y_pos)
    ax4.set_yticklabels(labels_p4, fontsize=8.8)
    ax4.set_xlim(0.71, 0.905)
    ax4.set_xlabel(r'Amplitude Parameter $S_8 = \sigma_8 \sqrt{\Omega_m / 0.300}$', color='#cbd5e1', fontsize=10.5)
    ax4.set_title(r'4-Tier Hardware Partition vs. Survey Landscape', color='#f8fafc', fontsize=11.5)

    bias_pct = (1.0 - CLUSTER_BIAS_1MB) * 100.0
    bias_badge = (
        rf"SZ-to-Lensing Bias $(1-b) = \Omega_{{m,\mathrm{{void}}}} / \Omega_{{m,\mathrm{{global}}}} = {CLUSTER_BIAS_1MB:.4f}$"
        "\n"
        rf"(-{bias_pct:.2f}% lensing deficit | $\Omega_\Lambda^{{\mathrm{{void}}}}={OMEGA_L_VOID:.4f}, \Omega_m^{{\mathrm{{core}}}}={OMEGA_M_CORE:.4f}$)"
    )
    ax4.text(
        0.97, 0.05, bias_badge, transform=ax4.transAxes,
        ha='right', va='bottom', fontsize=8.2, fontweight='bold', color='#fbbf24',
        bbox=dict(boxstyle='round,pad=0.48', facecolor='#0f172a', edgecolor='#f59e0b', alpha=0.9),
    )

    plt.tight_layout()

    # 5. Output DataFrame & Metrics Contract
    df_out = pd.DataFrame({
        "galaxy_name": names_eval,
        "distance_mpc": np.round(dists_eval, 4),
        "baryonic_mass_solar": mass_solar_eval,
        "rho_N_nodes_m3": rho_N,
        "theta_activation": np.round(theta_act, 6),
        "omega_m_local": np.round(omega_m_local, 6),
        "S8_local": np.round(s8_local, 6),
    })

    metrics_dict = {
        "Void Floor (S8_void)": f"{S8_VOID:.4f}",
        "Catalog Mean <S8>": f"{s8_mean:.4f}",
        "Global Baseline (S8_global)": f"{S8_GLOBAL:.4f}",
        "Cluster Core (S8_core)": f"{S8_CORE:.4f}",
        "SZ-Lensing Bias (1 - b)": f"{CLUSTER_BIAS_1MB:.4f} (-{bias_pct:.2f}%)",
    }

    return fig, df_out, metrics_dict


# ==========================================
# TWO-CONDITION COMPILATION GATE (16 REGIMES & COOL-CORE ICM)
# ==========================================
PHASE_LOCK_BENCHMARKS = [
    # Sub-threshold void floor (C1 = 0 -> PASSIVE)
    {"name": "Cosmic Void",                 "rho_N_m3": 0.05,   "T_K": 1.0e4,   "R_m": 20.0e3 * KPC_TO_METER, "M_enc_kg": 1.0e10 * M_SOLAR_KG, "mu_mol": 0.60, "v_bulk_ms": 0.0},
    {"name": "Diffuse IGM",                 "rho_N_m3": 0.50,   "T_K": 2.0e4,   "R_m": 5.0e3 * KPC_TO_METER,  "M_enc_kg": 1.0e11 * M_SOLAR_KG, "mu_mol": 0.60, "v_bulk_ms": 0.0},
    # Super-density passive I-phase plasma (C1 = 1, C2 = 0 -> PASSIVE)
    {"name": "Dense WHIM Filament",         "rho_N_m3": 20.0,   "T_K": 3.0e6,   "R_m": 1.5e3 * KPC_TO_METER,  "M_enc_kg": 1.0e13 * M_SOLAR_KG, "mu_mol": 0.60, "v_bulk_ms": 0.0},
    {"name": "Cluster Outskirts (r500)",    "rho_N_m3": 2.0e2,  "T_K": 3.5e7,   "R_m": 1.2e3 * KPC_TO_METER,  "M_enc_kg": 6.0e13 * M_SOLAR_KG, "mu_mol": 0.60, "v_bulk_ms": 0.0},
    {"name": "Hot Galactic Halo / CGM",     "rho_N_m3": 1.0e2,  "T_K": 2.0e6,   "R_m": 100.0 * KPC_TO_METER,  "M_enc_kg": 1.0e11 * M_SOLAR_KG, "mu_mol": 0.60, "v_bulk_ms": 0.0},
    {"name": "Non-Cool-Core Coma",          "rho_N_m3": 3.5e3,  "T_K": 9.5e7,   "R_m": 250.0 * KPC_TO_METER,  "M_enc_kg": 1.5e13 * M_SOLAR_KG, "mu_mol": 0.60, "v_bulk_ms": 0.0},
    {"name": "Cool-Core Perseus (20 kpc)",  "rho_N_m3": 4.0e4,  "T_K": 4.5e7,   "R_m": 20.0 * KPC_TO_METER,   "M_enc_kg": 6.0e11 * M_SOLAR_KG, "mu_mol": 0.60, "v_bulk_ms": 0.0},
    {"name": "Solar Wind (1 AU)",           "rho_N_m3": 7.0e6,  "T_K": 1.2e5,   "R_m": 1.0 * AU_TO_METER,     "M_enc_kg": 1.0 * M_SOLAR_KG,    "mu_mol": 0.60, "v_bulk_ms": 440.0e3},
    # Compiled phase-locked Gamma-phase matter (C1 = 1, C2 = 1 -> COMPILED)
    {"name": "Cool-Core Perseus (5 kpc)",   "rho_N_m3": 1.5e5,  "T_K": 1.0e7,   "R_m": 5.0 * KPC_TO_METER,    "M_enc_kg": 4.0e11 * M_SOLAR_KG, "mu_mol": 0.60, "v_bulk_ms": 0.0},
    {"name": "Cool-Core H-alpha Filaments", "rho_N_m3": 1.0e7,  "T_K": 1.0e4,   "R_m": 2.0 * KPC_TO_METER,    "M_enc_kg": 2.5e11 * M_SOLAR_KG, "mu_mol": 1.00, "v_bulk_ms": 0.0},
    {"name": "DDO 154 HI Disk",             "rho_N_m3": 2.0e5,  "T_K": 5.0e3,   "R_m": 3.0 * KPC_TO_METER,    "M_enc_kg": 3.5e8 * M_SOLAR_KG,  "mu_mol": 1.30, "v_bulk_ms": 0.0},
    {"name": "Cold ISM / GMC",              "rho_N_m3": 1.0e8,  "T_K": 20.0,    "R_m": 10.0e-3 * KPC_TO_METER,"M_enc_kg": 1.0e5 * M_SOLAR_KG,  "mu_mol": 2.33, "v_bulk_ms": 0.0},
    {"name": "99942 Apophis",               "rho_N_m3": 1.8e30, "T_K": 270.0,   "R_m": 185.0,                 "M_enc_kg": 6.1e10,              "mu_mol": 28.0, "v_bulk_ms": 0.0},
    {"name": "Earth",                       "rho_N_m3": 3.3e30, "T_K": 288.0,   "R_m": 6.3712e6,              "M_enc_kg": 5.9722e24,           "mu_mol": 28.0, "v_bulk_ms": 0.0},
    {"name": "Sun",                         "rho_N_m3": 8.4e29, "T_K": 5.778e3, "R_m": 6.9634e8,              "M_enc_kg": 1.0 * M_SOLAR_KG,    "mu_mol": 0.60, "v_bulk_ms": 0.0},
    {"name": "Crab Pulsar",                 "rho_N_m3": 3.0e44, "T_K": 1.0e8,   "R_m": 1.0e4,                 "M_enc_kg": 1.4 * M_SOLAR_KG,    "mu_mol": 1.00, "v_bulk_ms": 0.0},
]


def run_two_condition_phase_audit(
    custom_n_m3=5.0e4,
    custom_T_K=1.5e7,
    custom_R_kpc=10.0,
    custom_M_bcg_11=3.0,
    custom_mu=0.60,
    custom_v_bulk_kms=0.0,
):
    """
    Evaluates the zero-parameter Two-Condition Compilation Gate across 16 astrophysical
    benchmarks + a live user-configured cluster/gas parcel, and renders:
      Panel 1: 4-Quadrant Compilation Phase Diagram (log10 rho_N vs log10 Lambda_lock)
      Panel 2: Cool-Core (Perseus) vs Non-Cool-Core (Coma) radial Lambda_lock(r) profiles
    Returns (fig, df_phase, metrics_dict).
    """
    r_custom_m = float(custom_R_kpc) * KPC_TO_METER
    m_gas_custom = (4.0 / 3.0) * np.pi * (r_custom_m ** 3) * float(custom_n_m3) * M_PROTON
    m_enc_custom = float(custom_M_bcg_11) * 1.0e11 * M_SOLAR_KG + m_gas_custom

    entries = list(PHASE_LOCK_BENCHMARKS) + [{
        "name": "Live Custom Probe",
        "rho_N_m3": float(custom_n_m3),
        "T_K": float(custom_T_K),
        "R_m": r_custom_m,
        "M_enc_kg": m_enc_custom,
        "mu_mol": float(custom_mu),
        "v_bulk_ms": float(custom_v_bulk_kms) * 1.0e3,
    }]

    names_all = [e["name"] for e in entries]
    rho_arr = np.array([e["rho_N_m3"] for e in entries], dtype=float)
    temp_arr = np.array([e["T_K"] for e in entries], dtype=float)
    rad_arr = np.array([e["R_m"] for e in entries], dtype=float)
    mass_arr = np.array([e["M_enc_kg"] for e in entries], dtype=float)
    mu_arr = np.array([e["mu_mol"] for e in entries], dtype=float)
    vb_arr = np.array([e["v_bulk_ms"] for e in entries], dtype=float)

    gate_res = evaluate_two_condition_compilation_gate(
        rho_arr, temp_arr, rad_arr, mass_arr, mu_mol=mu_arr, v_bulk_ms=vb_arr
    )

    n_bench = len(PHASE_LOCK_BENCHMARKS)
    n_comp = int(np.sum(gate_res["theta_compiled"][:n_bench]))
    lam_custom = float(gate_res["lambda_lock"][-1])
    th_custom = int(gate_res["theta_compiled"][-1])

    # Panel 2 radial cluster profiles (1 kpc -> 300 kpc, strictly baryonic mass)
    r_kpc = np.geomspace(1.0, 300.0, 300)
    r_m = r_kpc * KPC_TO_METER

    n_pers = 2.0e5 * (1.0 + (r_kpc / 25.0) ** 2) ** (-1.0)
    t_pers = 8.0e6 + (6.0e7 - 8.0e6) * ((r_kpc / 25.0) ** 1.8) / (1.0 + (r_kpc / 25.0) ** 1.8)
    m_pers = (4.0e11 * M_SOLAR_KG * (r_kpc / (r_kpc + 1.0)) ** 2) + 0.5 * ((4.0 / 3.0) * np.pi * (r_m ** 3) * n_pers * M_PROTON)
    gate_pers = evaluate_two_condition_compilation_gate(n_pers, t_pers, r_m, m_pers, mu_mol=0.60, v_bulk_ms=0.0)
    lam_pers = gate_pers["lambda_lock"]

    n_coma = 3.5e3 * (1.0 + (r_kpc / 250.0) ** 2) ** (-1.125)
    t_coma = np.full_like(r_kpc, 9.5e7)
    m_coma = (1.5e11 * M_SOLAR_KG * (r_kpc / (r_kpc + 8.0)) ** 2) + ((4.0 / 3.0) * np.pi * (r_m ** 3) * n_coma * M_PROTON)
    gate_coma = evaluate_two_condition_compilation_gate(n_coma, t_coma, r_m, m_coma, mu_mol=0.60, v_bulk_ms=0.0)
    lam_coma = gate_coma["lambda_lock"]

    compiled_idx = np.where(gate_pers["theta_compiled"] == 1)[0]
    if len(compiled_idx) > 0:
        last_i = int(compiled_idx[-1])
        if last_i + 1 < len(r_kpc) and lam_pers[last_i] > 1.0 > lam_pers[last_i + 1]:
            frac = (lam_pers[last_i] - 1.0) / (lam_pers[last_i] - lam_pers[last_i + 1])
            r_cross_kpc = float(r_kpc[last_i] + frac * (r_kpc[last_i + 1] - r_kpc[last_i]))
        else:
            r_cross_kpc = float(r_kpc[last_i])
    else:
        r_cross_kpc = 0.0

    # ==========================================
    # 2-PANEL PHASE-LOCK DIAGNOSTIC FIGURE
    # ==========================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.2), facecolor='#0b1120')
    for ax in (ax1, ax2):
        ax.set_facecolor('#070c18')
        ax.tick_params(colors='#94a3b8')
        ax.grid(True, alpha=0.2, linestyle=':')

    # --- Panel 1: 4-Quadrant Compilation Phase Diagram ---
    log_rho = np.log10(np.maximum(rho_arr, 1e-12))
    log_lam = np.log10(np.maximum(gate_res["lambda_lock"], 1e-12))
    log_gc = float(np.log10(RHO_GC_THRESHOLD))

    ax1.set_xlim(-2.5, 47.0)
    ax1.set_ylim(-2.4, 5.2)

    ax1.axvspan(-2.5, log_gc, ymin=0.0, ymax=1.0, color='#38bdf8', alpha=0.08)
    ax1.axhspan(-2.4, 0.0, xmin=(log_gc - (-2.5)) / (47.0 - (-2.5)), xmax=1.0, color='#f59e0b', alpha=0.10)
    ax1.axhspan(0.0, 5.2, xmin=(log_gc - (-2.5)) / (47.0 - (-2.5)), xmax=1.0, color='#10b981', alpha=0.10)

    ax1.axvline(log_gc, color='#22d3ee', linestyle='--', linewidth=1.6, label=rf'Cond 1: $\rho_N = 4.0$ ($\log_{{10}}={log_gc:.3f}$)')
    ax1.axhline(0.0, color='#f43f5e', linestyle='--', linewidth=1.6, label=r'Cond 2: $\Lambda_{\mathrm{lock}} = 1.0$')

    ax1.text(-1.1, -2.05, "SUB-THRESHOLD\nVOID FLOOR\n(C1 = 0)", ha='center', va='bottom', fontsize=8.0, fontweight='bold', color='#60a5fa', alpha=0.9)
    ax1.text(25.0, -2.05, "SUPER-DENSITY PASSIVE\n(I-Phase Plasma: C1=1, C2=0)", ha='center', va='bottom', fontsize=8.2, fontweight='bold', color='#fbbf24', alpha=0.9)
    ax1.text(21.0, 4.85, "COMPILED MATTER\n(Γ-Phase: C1=1, C2=1)", ha='center', va='top', fontsize=8.5, fontweight='bold', color='#34d399', alpha=0.9)

    c1_b = gate_res["c1_density"][:n_bench]
    c2_b = gate_res["c2_phase_locked"][:n_bench]
    void_m = ~c1_b
    passive_hi_m = c1_b & (~c2_b)
    comp_m = c1_b & c2_b

    ax1.scatter(log_rho[:n_bench][void_m], log_lam[:n_bench][void_m], color='#60a5fa', s=70, edgecolors='#f8fafc', linewidths=0.8, zorder=5, label='Void / IGM (C1=0)')
    ax1.scatter(log_rho[:n_bench][passive_hi_m], log_lam[:n_bench][passive_hi_m], color='#f59e0b', s=75, marker='s', edgecolors='#f8fafc', linewidths=0.8, zorder=5, label='WHIM / Hot ICM / Wind (PASSIVE)')
    ax1.scatter(log_rho[:n_bench][comp_m], log_lam[:n_bench][comp_m], color='#10b981', s=75, edgecolors='#f8fafc', linewidths=0.8, zorder=5, label='Compiled Matter (Γ-Phase)')

    ax1.scatter([log_rho[-1]], [log_lam[-1]], color='#facc15', marker='*', s=250, edgecolors='#f8fafc', linewidths=0.9, zorder=8, label='Live Custom Probe')

    phase_offsets = {
        "Cosmic Void": (-4, 14, 'center', 'bottom', False),
        "Diffuse IGM": (12, 8, 'left', 'bottom', True),
        "Dense WHIM Filament": (78, 18, 'left', 'bottom', True),
        "Hot Galactic Halo / CGM": (8, -44, 'left', 'top', True),
        "Cluster Outskirts (r500)": (66, -16, 'left', 'top', True),
        "Non-Cool-Core Coma": (28, -64, 'left', 'top', True),
        "Cool-Core Perseus (20 kpc)": (42, -22, 'left', 'top', True),
        "Solar Wind (1 AU)": (16, -4, 'left', 'center', False),
        "Cool-Core Perseus (5 kpc)": (32, 18, 'left', 'bottom', True),
        "DDO 154 HI Disk": (14, 2, 'left', 'center', False),
        "Cool-Core H-alpha Filaments": (18, -14, 'left', 'top', True),
        "Cold ISM / GMC": (16, 12, 'left', 'bottom', True),
        "Sun": (-16, 12, 'right', 'bottom', True),
        "99942 Apophis": (14, -14, 'left', 'top', True),
        "Earth": (14, 12, 'left', 'bottom', True),
        "Crab Pulsar": (-12, -16, 'right', 'top', False),
        "Live Custom Probe": (36, 0, 'left', 'center', True),
    }
    for i_p, nm_p in enumerate(names_all):
        dx_p, dy_p, ha_p, va_p, arrow_p = phase_offsets.get(nm_p, (8, 8, 'left', 'bottom', False))
        arr_kw = dict(arrowstyle='-', color='#64748b', lw=0.7) if arrow_p else None
        lbl_col = '#facc15' if nm_p == "Live Custom Probe" else '#e2e8f0'
        fw = 'bold' if nm_p == "Live Custom Probe" else 'normal'
        ax1.annotate(
            nm_p,
            (log_rho[i_p], log_lam[i_p]),
            xytext=(dx_p, dy_p),
            textcoords='offset points',
            ha=ha_p, va=va_p,
            fontsize=7.2, color=lbl_col, fontweight=fw,
            arrowprops=arr_kw,
        )

    ax1.set_xlabel(r'Compiled Baryon Density $\log_{10}(\rho_N\ [\mathrm{nodes/m}^3])$', color='#cbd5e1', fontsize=10.5)
    ax1.set_ylabel(r'Phase-Lock Ratio $\log_{10}(\Lambda_{\mathrm{lock}})$', color='#cbd5e1', fontsize=10.5)
    ax1.set_title(r'Two-Condition Compilation Gate ($\Theta = \Theta(\rho_N - 4)\Theta(\Lambda_{\mathrm{lock}} - 1)$)', color='#f8fafc', fontsize=11.5)
    ax1.legend(loc='lower right', bbox_to_anchor=(0.985, 0.145), facecolor='#0f172a', edgecolor='none', labelcolor='#e2e8f0', fontsize=7.2)

    # --- Panel 2: Cool-Core vs Non-Cool-Core Cluster Transition Curve ---
    ax2.loglog(r_kpc, lam_pers, color='#10b981', linewidth=2.5, label=r'Cool-Core Perseus ($T_{\mathrm{core}} \to 8\times 10^6$ K + BCG)')
    ax2.loglog(r_kpc, lam_coma, color='#f59e0b', linewidth=2.4, linestyle='--', label=r'Non-Cool-Core Coma ($T = 9.5\times 10^7$ K, flat core)')
    ax2.axhline(1.0, color='#f43f5e', linestyle='--', linewidth=1.8, label=r'Phase-Lock Threshold ($\Lambda_{\mathrm{lock}} = 1.0$)')

    if r_cross_kpc > 1.0:
        ax2.axvspan(1.0, r_cross_kpc, color='#10b981', alpha=0.14, label=rf'Inner Compiled Core ($r \leq {r_cross_kpc:.2f}$ kpc)')
        ax2.axvline(r_cross_kpc, color='#34d399', linestyle=':', linewidth=1.6)
        ax2.scatter([r_cross_kpc], [1.0], color='#34d399', edgecolor='#f8fafc', s=75, zorder=6)

    if 1.0 <= float(custom_R_kpc) <= 300.0:
        ax2.scatter(
            [float(custom_R_kpc)], [max(lam_custom, 1e-4)],
            color='#facc15', marker='*', s=220, edgecolors='#f8fafc', linewidths=0.8, zorder=7,
            label=rf'Live Custom Probe ($\Lambda={lam_custom:.2e}$)'
        )

    icm_badge = (
        rf"Bohr Unlock Scale: $E_{{\mathrm{{Ryd}}}} = {E_RYD_J/1.602176634e-19:.2f}$ eV"
        "\n"
        rf"$T_{{\mathrm{{lock}}}} = \frac{{2 E_{{\mathrm{{Ryd}}}}}}{{3 k_B}} = {T_LOCK_RYD_K:,.0f}$ K"
        "\n"
        rf"Cool-Core Lock ($r \leq {r_cross_kpc:.2f}$ kpc) vs. Non-Cool-Core"
        "\n"
        r"drives falsifiable $(1-b)$ lensing split"
    )
    ax2.text(
        0.97, 0.96, icm_badge, transform=ax2.transAxes,
        ha='right', va='top', fontsize=8.0, color='#e2e8f0',
        bbox=dict(boxstyle='round,pad=0.45', facecolor='#0f172a', edgecolor='#10b981', alpha=0.90),
    )

    ax2.set_xlim(1.0, 300.0)
    ax2.set_ylim(2.0e-3, 10.0)
    ax2.set_xlabel(r'Clustercentric Radius $r$ [kpc]', color='#cbd5e1', fontsize=10.5)
    ax2.set_ylabel(r'Phase-Lock Ratio $\Lambda_{\mathrm{lock}}(r)$', color='#cbd5e1', fontsize=10.5)
    ax2.set_title(r'Cool-Core vs. Non-Cool-Core ICM Phase Transition', color='#f8fafc', fontsize=11.5)
    ax2.legend(loc='lower left', facecolor='#0f172a', edgecolor='none', labelcolor='#e2e8f0', fontsize=7.8)

    plt.tight_layout()

    # Output Audit DataFrame & Metrics
    df_phase = pd.DataFrame({
        "Regime": names_all,
        "rho_N_nodes_m3": rho_arr,
        "T_K": temp_arr,
        "r_grav": gate_res["r_grav"],
        "r_em": gate_res["r_em"],
        "Lambda_lock": gate_res["lambda_lock"],
        "Cond1_rho>=4": gate_res["c1_density"].astype(bool),
        "Cond2_Lock>=1": gate_res["c2_phase_locked"].astype(bool),
        "Verdict": ["COMPILED (1)" if v == 1 else "PASSIVE (0)" for v in gate_res["theta_compiled"]],
    })

    metrics_dict = {
        "Bohr Unlock Temp (T_lock)": f"{T_LOCK_RYD_K:,.0f} K (13.61 eV)",
        "Benchmarks Compiled": f"{n_comp} / {n_bench} Regimes",
        "Cool-Core Crossover (r_lock)": f"{r_cross_kpc:.2f} kpc (Perseus)",
        "Custom Probe Λ_lock": f"{lam_custom:.3e}",
        "Custom Gate Verdict": "COMPILED (Γ-Phase)" if th_custom == 1 else "PASSIVE (I-Phase)",
    }

    return fig, df_phase, metrics_dict
