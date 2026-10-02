import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st
from scipy.integrate import simpson
from astropy.coordinates import SkyCoord
import astropy.units as u

try:
    import healpy as hp
except ImportError:
    from .. import healpy as hp

from .data_loader import (
    C_LIGHT,
    H_GLOBAL,
    get_hlocal_batch
)

T_CMB_0 = 2.7255  # Kelvin
C_LIGHT_KM_S = C_LIGHT / 1000.0  # 299792.458 km/s

@st.cache_data(show_spinner=False)
def _compute_cmb_dipole_cached(
    _tree, _gal_positions, _nodes_per_galaxy,
    nside=32, d_max_mpc=35.0, sigma_mpc=2.0, h_global=H_GLOBAL,
    v_pec_km_s=369.0, l_kin_deg=264.0, b_kin_deg=48.0, n_steps=35
):
    npix = hp.nside2npix(nside)
    
    # Fully vectorized sky pixel unit vectors (npix, 3)
    pix_indices = np.arange(npix)
    theta_arr, phi_arr = hp.pix2ang(nside, pix_indices)
    pix_vecs = np.column_stack([
        np.sin(theta_arr) * np.cos(phi_arr),
        np.sin(theta_arr) * np.sin(phi_arr),
        np.cos(theta_arr)
    ])

    # 1. Kinematic Doppler Component
    kin_coord = SkyCoord(l=l_kin_deg*u.deg, b=b_kin_deg*u.deg, frame='galactic')
    v_kin_vec = v_pec_km_s * kin_coord.cartesian.get_xyz().value
    v_proj = np.dot(pix_vecs, v_kin_vec)
    delta_T_kin_mK = T_CMB_0 * (v_proj / C_LIGHT_KM_S) * 1000.0

    # 2. Batch Ray Tracing
    s_steps = np.linspace(0.1, d_max_mpc, n_steps)
    h_samples_matrix = np.zeros((npix, n_steps))

    for j, s in enumerate(s_steps):
        pts_at_s = s * pix_vecs
        sigma_adaptive = sigma_mpc * np.sqrt(1.0 + s / 10.0)
        
        h_samples_matrix[:, j] = get_hlocal_batch(
            _tree, _gal_positions, _nodes_per_galaxy, pts_at_s,
            sigma_mpc=sigma_adaptive, h_global=h_global, workers=-1
        )

    # 3. Path Integration across axis=1
    cz_excess = simpson(h_samples_matrix - h_global, x=s_steps, axis=1)
    delta_T_sub_raw_mK = -T_CMB_0 * (cz_excess / C_LIGHT_KM_S) * 1000.0

    # Monopole Subtraction
    delta_T_sub_pure_mK = delta_T_sub_raw_mK - np.mean(delta_T_sub_raw_mK)
    delta_T_total_mK = delta_T_kin_mK + delta_T_sub_pure_mK

    # Monopole + Dipole Fit
    design_matrix = np.hstack([np.ones((npix, 1)), pix_vecs])
    coeffs, _, _, _ = np.linalg.lstsq(design_matrix, delta_T_total_mK, rcond=None)
    fit_dip_total = coeffs[1:]
    amp_total_mK = float(np.linalg.norm(fit_dip_total))

    if amp_total_mK > 0:
        vx, vy, vz = fit_dip_total / amp_total_mK
        glat_tot = float(np.degrees(np.arcsin(np.clip(vz, -1.0, 1.0))))
        glon_tot = float(np.degrees(np.arctan2(vy, vx)) % 360.0)
    else:
        glon_tot, glat_tot = l_kin_deg, b_kin_deg

    return {
        "npix": npix,
        "delta_T_kin_mK": delta_T_kin_mK,
        "delta_T_sub_pure_mK": delta_T_sub_pure_mK,
        "delta_T_total_mK": delta_T_total_mK,
        "amp_total_mK": amp_total_mK,
        "glon_tot": glon_tot,
        "glat_tot": glat_tot
    }


def run_cmb_dipole_analysis(
    tree, gal_positions, nodes_per_galaxy,
    nside=32, d_max_mpc=35.0, sigma_mpc=2.0, h_global=H_GLOBAL,
    v_pec_km_s=369.0, l_kin_deg=264.0, b_kin_deg=48.0
):
    res = _compute_cmb_dipole_cached(
        tree, gal_positions, nodes_per_galaxy,
        nside=nside, d_max_mpc=d_max_mpc, sigma_mpc=sigma_mpc, h_global=h_global,
        v_pec_km_s=v_pec_km_s, l_kin_deg=l_kin_deg, b_kin_deg=b_kin_deg
    )

    npix = res["npix"]
    delta_T_kin_mK = res["delta_T_kin_mK"]
    delta_T_sub_pure_mK = res["delta_T_sub_pure_mK"]
    delta_T_total_mK = res["delta_T_total_mK"]
    amp_total_mK = res["amp_total_mK"]
    glon_tot = res["glon_tot"]
    glat_tot = res["glat_tot"]

    angles = hp.pix2ang(nside, np.arange(npix))
    glat_arr = 90.0 - np.degrees(angles[0])
    glon_arr = np.degrees(angles[1])

    df_table = pd.DataFrame({
        "pixel_id": np.arange(npix),
        "glon_deg": np.round(glon_arr, 2),
        "glat_deg": np.round(glat_arr, 2),
        "delta_T_kin_mK": np.round(delta_T_kin_mK, 4),
        "delta_T_sub_mK": np.round(delta_T_sub_pure_mK, 4),
        "delta_T_total_mK": np.round(delta_T_total_mK, 4)
    })

    fig = plt.figure(figsize=(11, 7), facecolor='#0b1120')
    plt.rcParams['text.color'] = '#e2e8f0'
    plt.rcParams['axes.labelcolor'] = '#94a3b8'

    hp.mollview(
        delta_T_total_mK,
        title=r"Substrate CMB Temperature Dipole Shift $\Delta T(\theta, \phi)$ [mK]",
        unit=r"$\Delta T$ [mK]",
        cmap="coolwarm",
        coord=["G"],
        fig=fig.number
    )
    hp.graticule()

    hp.projscatter(glon_tot, glat_tot, lonlat=True, coord='G', color='yellow', marker='*', s=220, label=f'Total Dipole Axis ({glon_tot:.1f}°, {glat_tot:.1f}°)')
    hp.projscatter(l_kin_deg, b_kin_deg, lonlat=True, coord='G', color='black', marker='o', s=110, label=f'Kinematic Axis ({l_kin_deg:.1f}°, {b_kin_deg:.1f}°)')
    hp.projscatter(280.0, 45.0, lonlat=True, coord='G', color='cyan', marker='^', s=110, label='Centaurus / Local Sheet Peak')

    plt.legend(loc='lower right', facecolor='#0f172a', edgecolor='none', fontsize=9)

    metrics = {
        "Kinematic Peak": f"{np.max(delta_T_kin_mK):.3f} mK",
        "Substrate Anisotropy Range": f"{np.min(delta_T_sub_pure_mK):.3f} to {np.max(delta_T_sub_pure_mK):.3f} mK",
        "Combined Dipole Amplitude": f"{amp_total_mK:.3f} mK",
        "Shifted Dipole Axis (l, b)": f"({glon_tot:.1f}°, {glat_tot:.1f}°)"
    }

    return fig, df_table, metrics