import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

from .data_loader import (
    M_PROTON,
    M_SOLAR_KG,
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
