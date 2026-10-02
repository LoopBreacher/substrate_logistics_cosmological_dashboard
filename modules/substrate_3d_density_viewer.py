import numpy as np
import pandas as pd
import plotly.graph_objects as go

from .data_loader import (
    get_hlocal,
    get_compiled_node_density,
    M_PROTON,
    M_SOLAR_KG,
    H_GLOBAL
)

def build_3d_cosmic_web_figure(
    tree, gal_positions, nodes_per_galaxy, gal_names,
    sigma_mpc=1.8, h_global=H_GLOBAL, grid_res=20, 
    max_galaxies=2500, show_isosurface=True, isosurface_opacity=0.25
):
    """
    Constructs an interactive 3D Plotly figure rendering galaxy scatter nodes 
    and volumetric filament isosurfaces.
    """
    N_total = len(gal_positions)
    
    # 1. Downsample galaxy nodes if catalog exceeds display cap for performance
    if N_total > max_galaxies:
        np.random.seed(42)
        indices = np.random.choice(N_total, size=max_galaxies, replace=False)
    else:
        indices = np.arange(N_total)

    sub_pos = gal_positions[indices]
    sub_nodes = nodes_per_galaxy[indices]
    sub_names = [gal_names[i] for i in indices]

    # Compute masses and distances for galaxy points
    sub_masses_solar = (sub_nodes * M_PROTON) / M_SOLAR_KG
    sub_dists = np.linalg.norm(sub_pos, axis=1)

    # Calculate H_local at each galaxy point
    sub_h_local = np.array([
        get_hlocal(tree, gal_positions, nodes_per_galaxy, pos, sigma_mpc=sigma_mpc, h_global=h_global)
        for pos in sub_pos
    ])

    # Dynamic scaling for galaxy node markers
    marker_sizes = np.clip(np.log10(np.maximum(sub_masses_solar, 1e6)) * 1.5 - 8.0, 3.0, 12.0)

    # 2. Build 3D Scatter Trace for Galaxy Catalog Points
    hover_texts = [
        f"<b>{name}</b><br>"
        f"Distance: {d:.2f} Mpc<br>"
        f"Position: ({pos[0]:.1f}, {pos[1]:.1f}, {pos[2]:.1f}) Mpc<br>"
        f"Baryonic Mass: {m:.2e} M☉<br>"
        f"Local Expansion H_local: {h:.2f} km/s/Mpc"
        for name, d, pos, m, h in zip(sub_names, sub_dists, sub_pos, sub_masses_solar, sub_h_local)
    ]

    scatter_trace = go.Scatter3d(
        x=sub_pos[:, 0],
        y=sub_pos[:, 1],
        z=sub_pos[:, 2],
        mode='markers',
        marker=dict(
            size=marker_sizes,
            color=sub_h_local,
            colorscale='Plasma',
            colorbar=dict(
                title=dict(text="H<sub>local</sub> [km/s/Mpc]", side="top"),
                ticks="outside",
                len=0.75,
                x=1.02
            ),
            opacity=0.85,
            showscale=True
        ),
        text=hover_texts,
        hoverinfo='text',
        name='Catalog Galaxies'
    )

    traces = [scatter_trace]

    # 3. Optional 3D Volumetric Filament Isosurface Mesh
    if show_isosurface and len(sub_pos) > 0:
        max_extent = float(np.max(np.abs(sub_pos))) * 0.95
        grid_x = np.linspace(-max_extent, max_extent, grid_res)
        grid_y = np.linspace(-max_extent, max_extent, grid_res)
        grid_z = np.linspace(-max_extent, max_extent, grid_res)

        X_g, Y_g, Z_g = np.meshgrid(grid_x, grid_y, grid_z, indexing='ij')
        grid_coords = np.vstack([X_g.ravel(), Y_g.ravel(), Z_g.ravel()]).T

        # Evaluate H_local field across 3D grid
        H_grid_vals = np.array([
            get_hlocal(tree, gal_positions, nodes_per_galaxy, pos, sigma_mpc=sigma_mpc, h_global=h_global)
            for pos in grid_coords
        ])

        # Render Isosurface Shells corresponding to overdense web structures
        h_min, h_max = np.min(H_grid_vals), np.max(H_grid_vals)
        if h_max > h_global + 0.1:
            isomin = h_global + 0.5
            isomax = h_max

            isosurface_trace = go.Isosurface(
                x=X_g.ravel(),
                y=Y_g.ravel(),
                z=Z_g.ravel(),
                value=H_grid_vals,
                isomin=isomin,
                isomax=isomax,
                surface_count=3,
                colorscale='Viridis',
                opacity=isosurface_opacity,
                showscale=False,
                caps=dict(x_show=False, y_show=False, z_show=False),
                name='Filament Web Isofields'
            )
            traces.append(isosurface_trace)

    # 4. Layout & Dark Theme Formatting
    fig = go.Figure(data=traces)
    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="#0b1120",
        plot_bgcolor="#070c18",
        margin=dict(l=0, r=0, b=0, t=30),
        scene=dict(
            xaxis=dict(title="X [Mpc]", backgroundcolor="#070c18", gridcolor="#1e293b", zerolinecolor="#38bdf8"),
            yaxis=dict(title="Y [Mpc]", backgroundcolor="#070c18", gridcolor="#1e293b", zerolinecolor="#38bdf8"),
            zaxis=dict(title="Z [Mpc]", backgroundcolor="#070c18", gridcolor="#1e293b", zerolinecolor="#38bdf8"),
            aspectmode='data'
        ),
        legend=dict(x=0.02, y=0.98, bgcolor="rgba(15, 23, 42, 0.7)")
    )

    # Build Summary DataFrame for active point set
    df_summary = pd.DataFrame({
        "galaxy_name": sub_names,
        "x_mpc": np.round(sub_pos[:, 0], 2),
        "y_mpc": np.round(sub_pos[:, 1], 2),
        "z_mpc": np.round(sub_pos[:, 2], 2),
        "distance_mpc": np.round(sub_dists, 2),
        "baryonic_mass_solar": sub_masses_solar,
        "H_local_kmsMpc": np.round(sub_h_local, 2)
    })

    metrics = {
        "Rendered Galaxies": f"{len(sub_pos):,}",
        "Peak Local Unspooling": f"{np.max(sub_h_local):.2f} km/s/Mpc",
        "Void Floor Unspooling": f"{np.min(sub_h_local):.2f} km/s/Mpc",
        "Mean Field Expansion": f"{np.mean(sub_h_local):.2f} km/s/Mpc"
    }

    return fig, df_summary, metrics

def run_3d_density_viewer_analysis(
    tree, gal_positions, nodes_per_galaxy, gal_names,
    sigma_mpc=1.8, h_global=H_GLOBAL, grid_res=20,
    max_galaxies=2500, show_isosurface=True, isosurface_opacity=0.25
):
    """
    Main Streamlit UI entry point for 3D Density Viewer.
    """
    return build_3d_cosmic_web_figure(
        tree, gal_positions, nodes_per_galaxy, gal_names,
        sigma_mpc=sigma_mpc, h_global=h_global, grid_res=grid_res,
        max_galaxies=max_galaxies, show_isosurface=show_isosurface,
        isosurface_opacity=isosurface_opacity
    )