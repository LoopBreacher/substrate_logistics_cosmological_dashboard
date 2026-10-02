import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

try:
    import healpy as hp
except ImportError:
    from .. import healpy as hp

from .data_loader import (
    get_hlocal_batch,
    H_GLOBAL
)

@st.cache_data(show_spinner=False)
def _compute_sky_forecast_cached(
    _tree, _gal_positions, _nodes_per_galaxy,
    nside=32, d_max_mpc=15.0, sigma_smooth_mpc=1.2, h_global=H_GLOBAL, n_steps=40
):
    npix = hp.nside2npix(nside)
    pix_indices = np.arange(npix)
    theta_arr, phi_arr = hp.pix2ang(nside, pix_indices)
    pix_vecs = np.column_stack([
        np.sin(theta_arr) * np.cos(phi_arr),
        np.sin(theta_arr) * np.sin(phi_arr),
        np.cos(theta_arr)
    ])

    s_steps = np.linspace(0.1, d_max_mpc, n_steps)
    h_samples_matrix = np.zeros((npix, n_steps))

    for j, s in enumerate(s_steps):
        pts_at_s = s * pix_vecs
        h_samples_matrix[:, j] = get_hlocal_batch(
            _tree, _gal_positions, _nodes_per_galaxy, pts_at_s,
            sigma_mpc=sigma_smooth_mpc, h_global=h_global, workers=-1
        )

    h_obs_map = np.mean(h_samples_matrix, axis=1)
    C_map = h_obs_map / h_global
    delta_mu_map = -5.0 * np.log10(C_map)

    return npix, h_obs_map, C_map, delta_mu_map


def run_anisotropic_sky_forecast(
    tree, gal_positions, nodes_per_galaxy,
    nside=32, d_max_mpc=15.0, sigma_smooth_mpc=1.2, h_global=H_GLOBAL
):
    npix, h_obs_map, C_map, delta_mu_map = _compute_sky_forecast_cached(
        tree, gal_positions, nodes_per_galaxy,
        nside=nside, d_max_mpc=d_max_mpc, sigma_smooth_mpc=sigma_smooth_mpc, h_global=h_global
    )

    fig = plt.figure(figsize=(14, 11))

    hp.mollview(
        h_obs_map,
        sub=(3, 1, 1),
        title=f"1. Observed Expansion Field $H_{{obs}}(\\theta, \\phi)$ [km/s/Mpc] ($D_{{max}}={d_max_mpc}$ Mpc)",
        unit="km/s/Mpc",
        cmap="plasma",
        min=h_global,
        max=float(np.max(h_obs_map)),
        flip="astro",
        coord=["G"]
    )
    hp.graticule()

    hp.mollview(
        C_map,
        sub=(3, 1, 2),
        title=r"2. Directional Correction Factor $C(\theta, \phi) = H_{obs} / H_{global}$",
        unit="Multiplier",
        cmap="magma",
        min=1.0,
        max=float(np.max(C_map)),
        flip="astro",
        coord=["G"]
    )
    hp.graticule()

    hp.mollview(
        delta_mu_map,
        sub=(3, 1, 3),
        title=r"3. Euclid & Rubin (LSST) SNe Ia Distance Modulus Residual $\Delta \mu(\theta, \phi)$ [mag]",
        unit=r"$\Delta \mu$ [mag]",
        cmap="viridis_r",
        min=float(np.min(delta_mu_map)),
        max=0.0,
        flip="astro",
        coord=["G"]
    )
    hp.graticule()

    plt.tight_layout()

    pix_indices = np.arange(npix)
    angles = hp.pix2ang(nside, pix_indices)
    theta, phi = angles[0], angles[1]
    glon = np.degrees(phi)
    glat = 90.0 - np.degrees(theta)

    lookup_df = pd.DataFrame({
        'pixel_id': pix_indices,
        'glon_deg': np.round(glon, 2),
        'glat_deg': np.round(glat, 2),
        'H_obs_kmsMpc': np.round(h_obs_map, 2),
        'Correction_Factor_C': np.round(C_map, 4),
        'Delta_mu_mag': np.round(delta_mu_map, 4)
    })

    metrics = {
        'mean_h0': float(np.mean(h_obs_map)),
        'min_c': float(np.min(C_map)),
        'max_c': float(np.max(C_map)),
        'min_dmu': float(np.min(delta_mu_map)),
        'max_dmu': float(np.max(delta_mu_map))
    }

    return fig, lookup_df, metrics