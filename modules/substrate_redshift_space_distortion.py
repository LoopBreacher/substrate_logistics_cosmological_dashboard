import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

from .data_loader import (
    H_GLOBAL,
    K_OMEGA,
    UNIT_CONV,
    MPC_TO_METER
)

# ==========================================
# 1. DYNAMIC CATALOG BASELINE HELPER
# ==========================================
def compute_dynamic_catalog_baseline(gal_positions, nodes_per_galaxy, default_fallback=0.03573):
    """
    Computes the true volumetric mean compiled node density from the loaded catalog:
    rho_N_baseline = (Total Catalog Nodes) / (Total Catalog Survey Sphere Volume)
    """
    if gal_positions is None or len(gal_positions) == 0:
        return default_fallback
        
    r_max = float(np.max(np.linalg.norm(gal_positions, axis=1)))
    if r_max <= 0:
        return default_fallback
        
    total_nodes = float(np.sum(nodes_per_galaxy))
    V_catalog_m3 = (4.0 / 3.0) * np.pi * (r_max * MPC_TO_METER)**3
    rho_baseline = total_nodes / V_catalog_m3
    
    return rho_baseline if rho_baseline > 0 else default_fallback

# ==========================================
# 2. MPC-SCALE CONTINUOUS RSD ENGINE
# ==========================================
def compute_catalog_rsd_field(
    tree, gal_positions, nodes_per_galaxy,
    alpha_inflow=120.0, h_global=H_GLOBAL, sigma_mpc=1.8, dr=0.2
):
    """
    Computes Substrate Logistics line-of-sight peculiar velocities directly from 
    the catalog density tree (rho_N) scaled against the dynamic catalog baseline.
    """
    N_gal = len(gal_positions)
    rho_N_baseline = compute_dynamic_catalog_baseline(gal_positions, nodes_per_galaxy)
    
    if N_gal > 5000:
        np.random.seed(42)
        eval_indices = np.random.choice(N_gal, size=5000, replace=False)
    else:
        eval_indices = np.arange(N_gal)
        
    eval_positions = gal_positions[eval_indices]
    N_eval = len(eval_positions)
    
    v_los_total = np.zeros(N_eval)
    v_infall_los = np.zeros(N_eval)
    v_virial_los = np.zeros(N_eval)
    rho_N_array = np.zeros(N_eval)
    
    for i, idx_gal in enumerate(eval_indices):
        pos = gal_positions[idx_gal]
        r_norm = np.linalg.norm(pos)
        if r_norm == 0:
            continue
        los_dir = pos / r_norm
        
        # 1. Evaluate compiled node density rho_N
        idx_ball = tree.query_ball_point(pos, r=3.0 * sigma_mpc)
        if idx_ball:
            dists = np.linalg.norm(gal_positions[idx_ball] - pos, axis=1)
            kernel = np.exp(-0.5 * (dists / sigma_mpc)**2) / ((2.0 * np.pi * sigma_mpc**2)**1.5)
            rho_N = np.sum(nodes_per_galaxy[idx_ball] * kernel) / (MPC_TO_METER**3)
        else:
            rho_N = 0.0
            
        rho_N_array[i] = rho_N
        
        # 2. Coherent Infall (Kaiser Effect): Evaluate grad(H_local) along line of sight
        pos_p = pos + dr * los_dir
        pos_m = pos - dr * los_dir
        
        idx_p = tree.query_ball_point(pos_p, r=3.0 * sigma_mpc)
        rho_p = np.sum(nodes_per_galaxy[idx_p] * np.exp(-0.5 * (np.linalg.norm(gal_positions[idx_p] - pos_p, axis=1) / sigma_mpc)**2) / ((2.0 * np.pi * sigma_mpc**2)**1.5)) / (MPC_TO_METER**3) if idx_p else 0.0
        h_p = h_global + (np.sqrt((8.0 * np.pi * K_OMEGA * rho_p) / 3.0) / UNIT_CONV) if rho_p > 0 else h_global
        
        idx_m = tree.query_ball_point(pos_m, r=3.0 * sigma_mpc)
        rho_m = np.sum(nodes_per_galaxy[idx_m] * np.exp(-0.5 * (np.linalg.norm(gal_positions[idx_m] - pos_m, axis=1) / sigma_mpc)**2) / ((2.0 * np.pi * sigma_mpc**2)**1.5)) / (MPC_TO_METER**3) if idx_m else 0.0
        h_m = h_global + (np.sqrt((8.0 * np.pi * K_OMEGA * rho_m) / 3.0) / UNIT_CONV) if rho_m > 0 else h_global
        
        grad_H_los = (h_p - h_m) / (2.0 * dr)
        v_infall_los[i] = alpha_inflow * grad_H_los
        
        # 3. Virial Dispersion (Fingers-of-God): Overdensity relative to dynamic baseline
        if rho_N > rho_N_baseline:
            np.random.seed(idx_gal + 42)
            overdensity_ratio = rho_N / rho_N_baseline
            v_disp_scale = 250.0 * np.sqrt(overdensity_ratio - 1.0)
            v_virial_los[i] = np.random.normal(0.0, v_disp_scale)
        else:
            v_virial_los[i] = 0.0
            
        v_los_total[i] = v_infall_los[i] + v_virial_los[i]
        
    los_dirs = eval_positions / np.linalg.norm(eval_positions, axis=1, keepdims=True)
    s_positions = eval_positions + (v_los_total[:, np.newaxis] / h_global) * los_dirs
    
    return eval_positions, s_positions, v_los_total, v_infall_los, v_virial_los, rho_N_array, rho_N_baseline

# ==========================================
# 3. 2D PAIR CORRELATION ENGINE \xi(\sigma, \pi)
# ==========================================
def compute_catalog_2d_correlation(s_positions, max_sep_mpc=20.0, n_bins=35):
    """
    Computes true 2-point galaxy pair correlation function DD(sigma, pi) 
    using 3D KD-Tree pair query across full survey volume.
    Returns Sigma, Pi, normalized xi_2d, and raw DD pair counts.
    """
    s_tree = cKDTree(s_positions)
    pairs = s_tree.query_pairs(r=max_sep_mpc, output_type='ndarray')
    
    if len(pairs) == 0:
        sigma_bins = np.linspace(0.1, max_sep_mpc, n_bins)
        pi_bins = np.linspace(-max_sep_mpc, max_sep_mpc, n_bins)
        Sigma, Pi = np.meshgrid(sigma_bins, pi_bins)
        return Sigma, Pi, np.zeros((n_bins, n_bins)), np.zeros((n_bins, n_bins), dtype=int)
        
    i_idx = pairs[:, 0]
    j_idx = pairs[:, 1]
    
    pos_i = s_positions[i_idx]
    pos_j = s_positions[j_idx]
    
    s_vec = pos_i - pos_j
    l_vec = (pos_i + pos_j) / 2.0
    l_norm = np.linalg.norm(l_vec, axis=1, keepdims=True)
    l_norm[l_norm == 0] = 1.0
    l_hat = l_vec / l_norm
    
    pi = np.abs(np.sum(s_vec * l_hat, axis=1))
    s_sq = np.sum(s_vec**2, axis=1)
    sigma = np.sqrt(np.maximum(0, s_sq - pi**2))
    
    sigma_bins = np.linspace(0.1, max_sep_mpc, n_bins)
    pi_bins = np.linspace(-max_sep_mpc, max_sep_mpc, n_bins)
    
    sigma_full = np.concatenate([sigma, sigma])
    pi_full = np.concatenate([pi, -pi])
    
    DD, xedges, yedges = np.histogram2d(sigma_full, pi_full, bins=[sigma_bins, pi_bins])
    xi_2d = DD / (np.max(DD) + 1e-8)
    
    Sigma, Pi = np.meshgrid((xedges[:-1] + xedges[1:])/2.0, (yedges[:-1] + yedges[1:])/2.0)
    return Sigma, Pi, xi_2d.T, DD.T.astype(int)

# ==========================================
# 4. MAIN STREAMLIT UI EXECUTION ENTRY POINT
# ==========================================
def run_redshift_space_distortions(
    tree, gal_positions, nodes_per_galaxy,
    alpha_inflow=120.0, sigma_mpc=1.8, h_global=H_GLOBAL
):
    """
    Main Streamlit UI entry point for Redshift-Space Distortions.
    Calculates Substrate velocity fields using dynamic catalog density baseline.
    """
    eval_pos, s_pos, v_los, v_inf, v_vir, rho_N, rho_baseline = compute_catalog_rsd_field(
        tree, gal_positions, nodes_per_galaxy,
        alpha_inflow=alpha_inflow, h_global=h_global, sigma_mpc=sigma_mpc
    )
    
    Sigma, Pi, xi_2d, DD_counts = compute_catalog_2d_correlation(s_pos, max_sep_mpc=20.0, n_bins=35)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), facecolor='#0b1120')
    for ax in (ax1, ax2):
        ax.set_facecolor('#070c18')
        ax.tick_params(colors='#94a3b8')
        ax.grid(True, alpha=0.2, linestyle=':')
        
    # Panel 1: Catalog Real Space vs Redshift Space Positions
    r_trans = np.linalg.norm(eval_pos[:, :2], axis=1)
    s_trans = np.linalg.norm(s_pos[:, :2], axis=1)
    
    ax1.scatter(r_trans, eval_pos[:, 2], color='#38bdf8', alpha=0.35, s=12, label='Real Space Coordinates')
    ax1.scatter(s_trans, s_pos[:, 2], color='#ef4444', alpha=0.35, s=12, label='Redshift Space Coordinates')
    
    ax1.set_xlabel(r'Transverse Radial Distance $\sigma$ [Mpc]', color='#cbd5e1')
    ax1.set_ylabel(r'Line-of-Sight Distance $Z$ [Mpc]', color='#cbd5e1')
    ax1.set_title(f'Catalog Galaxies (N={len(eval_pos):,}) Real vs. Redshift Space', color='#f8fafc', fontsize=12)
    ax1.legend(loc='upper right', facecolor='#0f172a', edgecolor='none', fontsize=9)

    # Panel 2: 2D Anisotropic Correlation Map \xi(\sigma, \pi)
    contour = ax2.contourf(Sigma, Pi, xi_2d, levels=18, cmap='magma')
    cbar = fig.colorbar(contour, ax=ax2)
    cbar.set_label(r'Pair Density Amplitude $\xi(\sigma, \pi)$', color='#cbd5e1')
    cbar.ax.tick_params(colors='#94a3b8')

    ax2.set_xlabel(r'Transverse Pair Separation $\sigma$ [Mpc]', color='#cbd5e1')
    ax2.set_ylabel(r'Line-of-Sight Pair Separation $\pi$ [Mpc]', color='#cbd5e1')
    ax2.set_title(r'Substrate Catalog 2D Correlation Map $\xi(\sigma, \pi)$', color='#f8fafc', fontsize=12)

    plt.tight_layout()

    # Enriched CSV Output DataFrame with Raw Pair Counts
    df_rsd = pd.DataFrame({
        'sigma_transverse_mpc': np.round(Sigma.ravel(), 3),
        'pi_los_mpc': np.round(Pi.ravel(), 3),
        'pair_count_DD': DD_counts.ravel(),
        'xi_correlation_amplitude': np.round(xi_2d.ravel(), 5)
    })

    overdense_count = int(np.sum(rho_N > rho_baseline))

    metrics = {
        "Evaluated Objects": f"{len(eval_pos):,}",
        "Dynamic Volumetric Baseline": f"{rho_baseline:.5f} nodes/m³",
        "Overdense Regions (ρ_N > ρ_base)": f"{overdense_count:,}",
        "Peak Infall Velocity": f"{float(np.max(np.abs(v_inf))):.1f} km/s",
        "1-σ Cluster Virial Disp": f"{float(np.std(v_vir[v_vir != 0])):.1f} km/s" if np.any(v_vir != 0) else "0.0 km/s",
        "Peak Outlier Dispersion": f"{float(np.max(np.abs(v_vir))):.1f} km/s"
    }

    return fig, df_rsd, metrics