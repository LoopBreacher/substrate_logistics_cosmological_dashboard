import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .data_loader import (
    get_hlocal,
    H_GLOBAL
)

ALPHA_VELOCITY_SCALE = 120.0  # Coupling scale (km/s per unit gradient)

# ==========================================
# 1. LINE-OF-SIGHT PECULIAR VELOCITY SOLVER
# ==========================================
def compute_line_of_sight_peculiar_velocity(
    pos_mpc, tree, gal_positions, nodes_per_galaxy,
    sigma_mpc=1.8, h_global=H_GLOBAL, alpha_scale=ALPHA_VELOCITY_SCALE, dr=0.2
):
    """
    Evaluates 3D spatial gradient +grad H_local at pos_mpc to compute 
    the inward peculiar velocity vector along the line-of-sight unit vector.
    """
    r_norm = np.linalg.norm(pos_mpc)
    if r_norm == 0:
        return 0.0
    los_dir = pos_mpc / r_norm
    
    grad_H = np.zeros(3)
    for i in range(3):
        pos_plus = np.array(pos_mpc, dtype=float)
        pos_minus = np.array(pos_mpc, dtype=float)
        pos_plus[i] += dr
        pos_minus[i] -= dr
        
        h_plus = get_hlocal(tree, gal_positions, nodes_per_galaxy, pos_plus, sigma_mpc=sigma_mpc, h_global=h_global)
        h_minus = get_hlocal(tree, gal_positions, nodes_per_galaxy, pos_minus, sigma_mpc=sigma_mpc, h_global=h_global)
        grad_H[i] = (h_plus - h_minus) / (2.0 * dr)
        
    v_vec = alpha_scale * grad_H
    v_pec_los = float(np.dot(v_vec, los_dir))
    return v_pec_los

# ==========================================
# 2. 2D PECULIAR VELOCITY FIELD CALCULATOR
# ==========================================
def compute_bulk_flow_field(
    tree, gal_positions, nodes_per_galaxy,
    grid_size=30, spatial_bounds=12.0, sigma_mpc=1.8, h_global=H_GLOBAL, alpha_scale=ALPHA_VELOCITY_SCALE
):
    """
    Calculates 2D spatial gradients and peculiar inflow velocity vectors across the galactic plane out to spatial_bounds.
    """
    x_range = np.linspace(-spatial_bounds, spatial_bounds, grid_size)
    y_range = np.linspace(-spatial_bounds, spatial_bounds, grid_size)
    X_grid, Y_grid = np.meshgrid(x_range, y_range)

    H_grid = np.zeros((grid_size, grid_size))
    for i in range(grid_size):
        for j in range(grid_size):
            pos = np.array([X_grid[i, j], Y_grid[i, j], 0.0])
            H_grid[i, j] = get_hlocal(tree, gal_positions, nodes_per_galaxy, pos, sigma_mpc=sigma_mpc, h_global=h_global)

    dx = x_range[1] - x_range[0]
    dH_dy, dH_dx = np.gradient(H_grid, dx)

    Vx_inflow = alpha_scale * dH_dx
    Vy_inflow = alpha_scale * dH_dy
    V_mag = np.sqrt(Vx_inflow**2 + Vy_inflow**2)

    return X_grid, Y_grid, H_grid, Vx_inflow, Vy_inflow, V_mag

# ==========================================
# 3. MAIN UI EXECUTION & PLOTTING ENGINE
# ==========================================
def run_bulk_flow_analysis(
    tree, gal_positions, nodes_per_galaxy, X_gal, Y_gal, Z_gal,
    grid_size=30, spatial_bounds=12.0, sigma_mpc=1.8, h_global=H_GLOBAL, alpha_scale=ALPHA_VELOCITY_SCALE
):
    """
    Main UI entry point for Tab 4. Renders 2D Peculiar Velocity Field matching active spatial bounds.
    """
    X_grid, Y_grid, H_grid, Vx, Vy, V_mag = compute_bulk_flow_field(
        tree, gal_positions, nodes_per_galaxy,
        grid_size=grid_size, spatial_bounds=spatial_bounds,
        sigma_mpc=sigma_mpc, h_global=h_global, alpha_scale=alpha_scale
    )
    
    fig, ax = plt.subplots(figsize=(10, 7), facecolor='#0b1120')
    ax.set_facecolor('#070c18')

    contour = ax.contourf(X_grid, Y_grid, H_grid, levels=25, cmap='plasma', alpha=0.85)
    cb = fig.colorbar(contour, ax=ax)
    cb.set_label(r'Local Unspooling Rate $H_{\mathrm{local}}(\vec{r})$ [km/s/Mpc]', color='#cbd5e1')
    cb.ax.tick_params(colors='#94a3b8')

    # Dynamically scale quiver arrows based on spatial extent
    quiver_scale = 150.0 * (spatial_bounds / 12.0)
    ax.quiver(X_grid, Y_grid, Vx, Vy, V_mag, cmap='cool', scale=quiver_scale, width=0.0035, headwidth=4)

    # Slice thickness scales with field size for clear visualization
    z_thickness = max(3.0, 0.05 * spatial_bounds)
    plane_mask = (np.abs(Z_gal) <= z_thickness) & (np.abs(X_gal) <= spatial_bounds) & (np.abs(Y_gal) <= spatial_bounds)
    
    ax.scatter(X_gal[plane_mask], Y_gal[plane_mask], c='white', edgecolors='black', s=20, alpha=0.8,
               label=rf'Catalog Galaxies ($|Z| \leq {z_thickness:.1f}$ Mpc)')

    # Strict axis bounds lock view frame to vector field
    ax.set_xlim(-spatial_bounds, spatial_bounds)
    ax.set_ylim(-spatial_bounds, spatial_bounds)

    ax.set_xlabel('Local X [Mpc]', color='#cbd5e1', fontsize=11)
    ax.set_ylabel('Local Y [Mpc]', color='#cbd5e1', fontsize=11)
    ax.set_title(rf'Substrate Logistics: Bulk Flow Dynamics ($R \leq {spatial_bounds:.0f}$ Mpc)',
                 color='#f8fafc', fontsize=12)
    ax.tick_params(colors='#94a3b8')
    ax.grid(True, alpha=0.2, linestyle=':')
    ax.legend(loc='upper right', facecolor='#0f172a', edgecolor='none')

    plt.tight_layout()

    # Construct exportable DataFrame
    df_table = pd.DataFrame({
        "X_Mpc": np.round(X_grid.ravel(), 2),
        "Y_Mpc": np.round(Y_grid.ravel(), 2),
        "H_local_kmsMpc": np.round(H_grid.ravel(), 3),
        "Vx_inflow_kms": np.round(Vx.ravel(), 2),
        "Vy_inflow_kms": np.round(Vy.ravel(), 2),
        "V_magnitude_kms": np.round(V_mag.ravel(), 2)
    })

    metrics = {
        "Mean Drift Velocity": f"{np.mean(V_mag):.2f} km/s",
        "Peak Filament Infall": f"{np.max(V_mag):.2f} km/s",
        "Void Core Drift": f"{np.min(V_mag):.2f} km/s"
    }

    return fig, df_table, metrics