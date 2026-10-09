import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from astropy.coordinates import SkyCoord
import astropy.units as u

from modules.data_loader import (
    H_GLOBAL,
    M_PROTON,
    M_SOLAR_KG,
    get_hlocal,
    get_hlocal_batch
)

# Standalone Physical Density Environment Presets
PHYSICAL_ENVIRONMENTS_DE = {
    "Macro Volume Average (Macro Scale / σ=13.3 Mpc)": {
        "sigma_mpc": 13.3,
        "description": "Void-dominated 40 Mpc sphere average (Cosmic web mean)",
        "expected_range": "w0 ≈ -0.972 (Relaxation Floor)"
    },
    "Intergalactic Filament Corridor (Filament Scale / σ=2.5 Mpc)": {
        "sigma_mpc": 2.5,
        "description": "Smooth cosmic web filament path profile",
        "expected_range": "w0 ≈ -0.940 (Filament Transition)"
    },
    "Galactic Host Environment (Stellar Scale / σ=1.8 Mpc)": {
        "sigma_mpc": 1.8,
        "description": "Standard galactic host environment (Local Sheet)",
        "expected_range": "w0 ≈ -0.916 (Host Galaxy Base)"
    },
    "Supercluster Core Node (High-Density / σ=1.2 Mpc)": {
        "sigma_mpc": 1.2,
        "description": "Dense supercluster filament junction",
        "expected_range": "w0 ≈ -0.820 (Supercluster Boost)"
    },
    "Virial Cluster Core Peak (Peak Scale / σ=0.8 Mpc)": {
        "sigma_mpc": 0.8,
        "description": "Virialized cluster core peak density",
        "expected_range": "w0 ≈ -0.730 (Maximum Virial Shift)"
    }
}

def run_3d_dark_energy_analysis(
    tree,
    gal_positions,
    nodes_per_galaxy,
    gal_names,
    preset_key="Galactic Host Environment (Stellar Scale / σ=1.8 Mpc)",
    target_idx=0,
    target_ra=201.0,
    target_dec=-29.8,
    max_dist_mpc=200.0,
    grid_n=40,
    sigma_mpc=1.8,
    h_global=H_GLOBAL
):
    """
    Module 16: 3D Dark Energy Field Theory Solver (4-Panel Grid).
    Computes point-wise w0(D), wa(D), H_local(D), cumulative w_LOS(D), and 2D spatial maps.
    """
    plt.style.use('dark_background')
    fig = plt.figure(figsize=(16, 11), facecolor='#0f172a')
    gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.28)
    
    ax1 = fig.add_subplot(gs[0, 0]) # Top Left: CPL Redshift Evolution w(z)
    ax2 = fig.add_subplot(gs[0, 1]) # Top Right: Pointwise Sightline Dynamics H_local(D) & w0(D)
    ax3 = fig.add_subplot(gs[1, 0]) # Bottom Left: Cumulative Path-Integrated w_LOS(D)
    ax4 = fig.add_subplot(gs[1, 1]) # Bottom Right: 2D Spatial Map w0(X, Y)
    
    # -------------------------------------------------------------------------
    # PANEL 1: TARGET GALAXY / CLUSTER CPL REDSHIFT TRAJECTORY
    # -------------------------------------------------------------------------
    target_pos = gal_positions[target_idx]
    target_name = gal_names[target_idx]
    
    h_loc_target = get_hlocal(tree, gal_positions, nodes_per_galaxy, target_pos, sigma_mpc=sigma_mpc, h_global=h_global)
    delta_h_target = h_loc_target - h_global
    
    w0_target = -1.0 + (delta_h_target / h_global)
    wa_target = -1.0 - w0_target  # Structural invariant: w0 + wa = -1.0
    
    a_vec = np.linspace(0.1, 1.0, 100)
    z_vec = (1.0 / a_vec) - 1.0
    w_a_curve = w0_target + wa_target * (1.0 - a_vec)
    
    ax1.plot(z_vec, w_a_curve, color='#38bdf8', lw=2.5, label=f"Local Target: {target_name} ($w_0={w0_target:.3f}$, $w_a={wa_target:.3f}$)")
    ax1.axhline(-1.0, color='#94a3b8', ls='--', alpha=0.7, label=r"Uncompiled Void Floor ($w = -1.0$, $w_a = 0.0$)")
    
    ax1.axhspan(-0.95, -0.70, color='#38bdf8', alpha=0.08, label=r"Compiled Density Shift Range ($w_0 > -1.0$)")
    ax1.fill_between(z_vec, -1.0, w_a_curve, color='#38bdf8', alpha=0.15)
    
    ax1.set_title(f"CPL Redshift Trajectory $w(z)$ | {preset_key.split('(')[0].strip()}", fontsize=11, fontweight='bold', color='#f8fafc')
    ax1.set_xlabel("Cosmological Redshift $z$", fontsize=9, color='#cbd5e1')
    ax1.set_ylabel("Equation of State $w(z)$", fontsize=9, color='#cbd5e1')
    ax1.set_xlim(0, 3.0)
    ax1.set_ylim(-1.05, max(-0.65, w0_target + 0.05))
    ax1.grid(True, linestyle=':', alpha=0.3, color='#334155')
    ax1.legend(loc='upper right', facecolor='#1e293b', edgecolor='#334155', fontsize=8)
    
    # -------------------------------------------------------------------------
    # PANEL 2: POINT-WISE SIGHTLINE DYNAMICS (H_local & Local w0)
    # -------------------------------------------------------------------------
    c_target = SkyCoord(ra=target_ra*u.deg, dec=target_dec*u.deg, frame='icrs')
    ray_dir = c_target.cartesian.xyz.value
    ray_dir /= np.linalg.norm(ray_dir)
    
    d_steps = np.linspace(0.5, max_dist_mpc, 100)
    pos_ray = np.outer(d_steps, ray_dir)
    
    h_ray_batch = get_hlocal_batch(tree, gal_positions, nodes_per_galaxy, pos_ray, sigma_mpc=sigma_mpc, h_global=h_global)
    delta_h_ray = h_ray_batch - h_global
    
    w0_pointwise = -1.0 + (delta_h_ray / h_global)
    wa_pointwise = -(delta_h_ray / h_global)
    
    color_h = '#38bdf8'
    color_w = '#f43f5e'
    
    ax2.plot(d_steps, h_ray_batch, color=color_h, lw=2.2, label=r"Local $H_\mathrm{local}(D)$")
    ax2.axhline(h_global, color='#94a3b8', ls='--', alpha=0.6, label=rf"Global Baseline $H_0 = {h_global:.2f}$")
    ax2.set_xlabel("Sightline Depth $D$ [Mpc]", fontsize=9, color='#cbd5e1')
    ax2.set_ylabel(r"Local Expansion $H_\mathrm{local}$ [km/s/Mpc]", fontsize=9, color=color_h)
    ax2.tick_params(axis='y', labelcolor=color_h)
    
    ax2_twin = ax2.twinx()
    ax2_twin.plot(d_steps, w0_pointwise, color=color_w, lw=2.0, ls='-.', label=r"Pointwise $w_0(D)$")
    ax2_twin.set_ylabel(r"Pointwise Equation of State $w_0(D)$", fontsize=9, color=color_w)
    ax2_twin.tick_params(axis='y', labelcolor=color_w)
    ax2_twin.set_ylim(-1.02, max(-0.70, np.max(w0_pointwise) + 0.03))
    
    lines_1, labels_1 = ax2.get_legend_handles_labels()
    lines_2, labels_2 = ax2_twin.get_legend_handles_labels()
    ax2.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right', facecolor='#1e293b', edgecolor='#334155', fontsize=8)
    
    ax2.set_title(rf"Pointwise Sightline Expansion $H_\mathrm{{local}}(D)$ & Local $w_0(D)$", fontsize=11, fontweight='bold', color='#f8fafc')
    ax2.grid(True, linestyle=':', alpha=0.3, color='#334155')
    
    # -------------------------------------------------------------------------
    # PANEL 3: CUMULATIVE PATH-INTEGRATED w_LOS(D) PROFILE
    # -------------------------------------------------------------------------
    cum_delta_h = np.cumsum(delta_h_ray) * (d_steps[1] - d_steps[0])
    w_los_steps = -1.0 + (cum_delta_h / (d_steps * h_global))
    
    ax3.plot(d_steps, w_los_steps, color='#a855f7', lw=2.5, label=f"Path Integral (RA={target_ra:.1f}°, Dec={target_dec:.1f}°)")
    ax3.axhline(-1.0, color='#94a3b8', ls='--', alpha=0.7, label=r"Pure Void Floor ($w = -1.0$)")
    
    ax3.set_title(rf"Path-Integrated Cumulative $w_\mathrm{{LOS}}(D)$ Profile", fontsize=11, fontweight='bold', color='#f8fafc')
    ax3.set_xlabel("Sightline Depth $D$ [Mpc]", fontsize=9, color='#cbd5e1')
    ax3.set_ylabel(rf"Integrated $w_\mathrm{{LOS}}(D)$", fontsize=9, color='#cbd5e1')
    ax3.set_ylim(-1.02, max(-0.65, np.max(w_los_steps) + 0.03))
    ax3.grid(True, linestyle=':', alpha=0.3, color='#334155')
    ax3.legend(loc='upper right', facecolor='#1e293b', edgecolor='#334155', fontsize=8)
    
    # -------------------------------------------------------------------------
    # PANEL 4: DYNAMIC 2D CARTESIAN SLICE MAP OF w0(X, Y)
    # -------------------------------------------------------------------------
    margin = 20.0
    slice_bounds = max(200.0, float(np.abs(target_pos[0])) + margin, float(np.abs(target_pos[1])) + margin)
    slice_bounds = min(slice_bounds, max_dist_mpc)
    
    x_grid = np.linspace(-slice_bounds, slice_bounds, grid_n)
    y_grid = np.linspace(-slice_bounds, slice_bounds, grid_n)
    XX, YY = np.meshgrid(x_grid, y_grid)
    ZZ = np.zeros_like(XX)
    
    grid_coords = np.vstack([XX.ravel(), YY.ravel(), ZZ.ravel()]).T
    h_grid_batch = get_hlocal_batch(tree, gal_positions, nodes_per_galaxy, grid_coords, sigma_mpc=sigma_mpc, h_global=h_global)
    delta_h_grid = h_grid_batch - h_global
    
    W0_grid = (-1.0 + (delta_h_grid / h_global)).reshape(grid_n, grid_n)
    
    im = ax4.imshow(
        W0_grid,
        extent=[-slice_bounds, slice_bounds, -slice_bounds, slice_bounds],
        origin='lower',
        cmap='magma',
        aspect='equal'
    )
    
    cbar = fig.colorbar(im, ax=ax4, orientation='vertical', pad=0.02)
    cbar.set_label(rf"Present-Day Equation of State $w_0(\vec{{r}})$", fontsize=9, color='#cbd5e1')
    cbar.ax.tick_params(colors='#cbd5e1', labelsize=8)
    
    z_mask = np.abs(gal_positions[:, 2]) < (slice_bounds * 0.3)
    xy_mask = (np.abs(gal_positions[:, 0]) < slice_bounds) & (np.abs(gal_positions[:, 1]) < slice_bounds)
    slice_gal_mask = z_mask & xy_mask
    
    ax4.scatter(
        gal_positions[slice_gal_mask, 0],
        gal_positions[slice_gal_mask, 1],
        s=12,
        color='#38bdf8',
        alpha=0.6,
        edgecolors='none',
        label=f"Catalog Galaxies (N={np.sum(slice_gal_mask):,})"
    )
    
    ax4.scatter([target_pos[0]], [target_pos[1]], s=140, color='#22c55e', marker='*', edgecolors='#ffffff', lw=1.5, zorder=5, label=f"Target: {target_name}")
    
    ax4.set_title(rf"2D Equatorial Slice of 3D Dark Energy Field $w_0(X, Y, Z=0)$ | Slice Extent $\pm {slice_bounds:.1f}$ Mpc", fontsize=11, fontweight='bold', color='#f8fafc')
    ax4.set_xlabel("Supergalactic X [Mpc]", fontsize=9, color='#cbd5e1')
    ax4.set_ylabel("Supergalactic Y [Mpc]", fontsize=9, color='#cbd5e1')
    ax4.grid(True, linestyle=':', alpha=0.2, color='#ffffff')
    ax4.legend(loc='upper right', facecolor='#1e293b', edgecolor='#334155', fontsize=8)
    
    for ax in [ax1, ax2, ax3, ax4]:
        ax.tick_params(colors='#cbd5e1', labelsize=8)
        for spine in ax.spines.values():
            spine.set_color('#334155')
            
    # Complete Output Dataframe with Point-wise and Integrated Values
    df_out = pd.DataFrame({
        "Sightline_Depth_Mpc": d_steps,
        "H_local_Pointwise_km_s_Mpc": h_ray_batch,
        "Delta_H_Pointwise_km_s_Mpc": delta_h_ray,
        "w0_Pointwise": w0_pointwise,
        "wa_Pointwise": wa_pointwise,
        "w_LOS_Integrated": w_los_steps,
        "wa_LOS_Integrated": -1.0 - w_los_steps
    })
    
    metrics_dict = {
        "Target Local w0": f"{w0_target:.4f}",
        "Target Local wa": f"{wa_target:.4f}",
        "Local Boost Shift": f"+{delta_h_target / h_global:.4f}",
        "Structural Invariant (w0 + wa)": f"{w0_target + wa_target:.4f}",
        "Sightline Mean w_LOS": f"{np.mean(w_los_steps):.4f}"
    }
    
    return fig, df_out, metrics_dict