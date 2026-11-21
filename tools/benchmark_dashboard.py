"""
Headless Automated Benchmark & Verification Suite for the Substrate Logistics Dashboard.

Verifies:
1. Exact first-principles closed-form constants (zero CODATA G, zero rounded literals).
2. Orthogonal J2000 Galactic -> ICRS sky alignment matrix (R_GAL_TO_ICRS).
3. Headless execution across all 18 computational modules with real local catalog caches,
   verifying non-empty DataFrames, zero NaN/Inf numeric outputs, and clean Matplotlib renders.
4. Physics spot checks across rotation curves (Freeman vs Spherical + DDO 154 deep underflow),
   garbage collection purge rates & Quantum Radius, Lensing (1 - b) cluster bias, S8 weak-lensing
   tension (Smooth & Strict), P_delta 4-regime confinement, atomic clock drifts, SNe Ia agreement,
   and shifted CMB dipole axis.
"""

import os
import sys
import time
import warnings
import logging

# Ensure repository root is importable regardless of launch working directory
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import cKDTree

# Silence noisy bare-mode Streamlit cache warnings when executed outside `streamlit run`
for logger_name in ("streamlit", "streamlit.runtime.caching", "streamlit.runtime.scriptrunner_utils"):
    logging.getLogger(logger_name).setLevel(logging.ERROR)
warnings.filterwarnings("ignore", message=".*No runtime found.*")
warnings.filterwarnings("ignore", message=".*to view this Streamlit app.*")

from modules.data_loader import (
    C_LIGHT,
    M_PROTON,
    LAMBDA_BAR_P,
    C_DELTA,
    T_OMEGA,
    K_OMEGA,
    SIGMA_0_LEDGER,
    H_GLOBAL,
    R_BAO_MPC,
    A_OMEGA,
    OMEGA_M_GLOBAL,
    S8_VOID,
    S8_GLOBAL,
    S8_CORE,
    CLUSTER_BIAS_1MB,
    P_DELTA_MACRO_PA,
    B_QCD_MEV_FM3,
    M_PION_MEV,
    R_YUKAWA_FM,
    DISK_TRUNCATION_RATIO,
    R_GAL_TO_ICRS,
    galactic_to_icrs_vecs,
)
from modules.substrate_directional_correction_field import run_anisotropic_sky_forecast
from modules.line_of_sight_kinematics import run_line_of_sight_kinematics_analysis
from modules.redshift_path_integral import run_redshift_path_analysis
from modules.substrate_hubble_diagram_fitter import run_hubble_fitter_analysis
from modules.substrate_kbc_void_profile import run_kbc_void_analysis
from modules.substrate_bulk_flow_field import run_bulk_flow_analysis
from modules.substrate_void_boundary_dynamics import run_void_boundary_analysis
from modules.substrate_redshift_space_distortion import run_redshift_space_distortions
from modules.substrate_rotation_curves import run_rotation_curve_analysis
from modules.substrate_3d_density_viewer import run_3d_density_viewer_analysis
from modules.substrate_3d_dark_energy import run_3d_dark_energy_analysis
from modules.substrate_bao_ruler_warping import run_bao_warping_analysis
from modules.substrate_cmb_dipole_anomaly import run_cmb_dipole_analysis
from modules.substrate_galactic_lensing import run_lensing_analysis
from modules.spatial_garbage_collection_engine import run_spatial_garbage_collection_analysis
from modules.cpu_bandwidth_throttling import run_cpu_bandwidth_throttling_analysis
from modules.substrate_s8_tension import run_s8_tension_analysis
from modules.substrate_topological_drag_floor import run_topological_drag_floor_analysis


def load_local_catalog():
    """Loads the largest available local .npz catalog cache without network access."""
    candidates = [
        ("Cosmicflows-4", "cf4_catalog_cache.npz", 150.0, 2.0),
        ("2M++ Infrared", "twompp_catalog_cache.npz", 200.0, 2.5),
        ("UNG 2013", "ung_catalog_cache.npz", 35.0, 1.8),
    ]
    modules_dir = os.path.join(REPO_ROOT, "modules")
    for label, fname, d_max, sigma_default in candidates:
        fpath = os.path.join(modules_dir, fname)
        if os.path.exists(fpath):
            data = np.load(fpath)
            gal_pos_all = data["gal_positions"]
            nodes_all = data["nodes_per_galaxy"]
            X_all, Y_all, Z_all = data["X"], data["Y"], data["Z"]
            if "gal_names" in data:
                names_all = data["gal_names"]
            else:
                names_all = np.array([f"Galaxy #{i}" for i in range(len(gal_pos_all))], dtype=str)

            dists = np.linalg.norm(gal_pos_all, axis=1)
            mask = (dists > 0.1) & (dists <= d_max)
            gal_pos = gal_pos_all[mask]
            nodes = nodes_all[mask]
            X, Y, Z = X_all[mask], Y_all[mask], Z_all[mask]
            names = names_all[mask]
            tree = cKDTree(gal_pos)
            return label, d_max, sigma_default, tree, gal_pos, nodes, X, Y, Z, names

    # Deterministic synthetic catalog fallback (1,000 galaxies) if no local cache exists
    rng = np.random.default_rng(12345)
    n_syn = 1000
    r_syn = rng.uniform(1.0, 35.0, size=n_syn)
    cos_t = rng.uniform(-1.0, 1.0, size=n_syn)
    phi = rng.uniform(0.0, 2.0 * np.pi, size=n_syn)
    sin_t = np.sqrt(1.0 - cos_t**2)
    X = r_syn * sin_t * np.cos(phi)
    Y = r_syn * sin_t * np.sin(phi)
    Z = r_syn * cos_t
    gal_pos = np.column_stack([X, Y, Z])
    nodes = rng.uniform(1e66, 1e69, size=n_syn)
    names = np.array([f"Synthetic-{i}" for i in range(n_syn)], dtype=str)
    tree = cKDTree(gal_pos)
    return "Synthetic Local Volume", 35.0, 1.8, tree, gal_pos, nodes, X, Y, Z, names


def verify_constants_and_frames():
    """Prints closed-form constants and asserts mathematical/astrometric identities."""
    print("## 1. First-Principles Constants & Galactic -> ICRS Alignment\n")
    invariants = [
        ("C_DELTA", f"{C_DELTA:.10e}", "dimensionless"),
        ("T_OMEGA", f"{T_OMEGA:.10e}", "dimensionless"),
        ("K_OMEGA", f"{K_OMEGA:.10e}", "m^3/(node*s^2)"),
        ("SIGMA_0_LEDGER", f"{SIGMA_0_LEDGER:.10e}", "J/m^3"),
        ("H_GLOBAL", f"{H_GLOBAL:.6f}", "km/s/Mpc"),
        ("R_BAO_MPC", f"{R_BAO_MPC:.6f}", "Mpc"),
        ("A_OMEGA", f"{A_OMEGA:.10e}", "m/s^2"),
        ("S8_VOID", f"{S8_VOID:.6f}", "dimensionless"),
        ("S8_GLOBAL", f"{S8_GLOBAL:.6f}", "dimensionless"),
        ("CLUSTER_BIAS_1MB", f"{CLUSTER_BIAS_1MB:.6f}", "ratio"),
        ("P_DELTA_MACRO_PA", f"{P_DELTA_MACRO_PA:.10e}", "Pa"),
        ("B_QCD_MEV_FM3", f"{B_QCD_MEV_FM3:.6f}", "MeV/fm^3"),
        ("M_PION_MEV", f"{M_PION_MEV:.6f}", "MeV/c^2"),
        ("R_YUKAWA_FM", f"{R_YUKAWA_FM:.6f}", "fm"),
    ]
    print("| Constant | Exact Value | Units |")
    print("| :--- | :--- | :--- |")
    for name, val, unit in invariants:
        print(f"| `{name}` | `{val}` | {unit} |")

    orth_err = float(np.max(np.abs(R_GAL_TO_ICRS.T @ R_GAL_TO_ICRS - np.eye(3))))
    assert orth_err < 1e-14, f"R_GAL_TO_ICRS orthogonality error {orth_err}"

    # Verify against Astropy SkyCoord along the Solar kinematic dipole direction
    c_gal = SkyCoord(l=264.0 * u.deg, b=48.0 * u.deg, frame="galactic")
    v_gal = c_gal.cartesian.xyz.value
    v_icrs_astropy = c_gal.icrs.cartesian.xyz.value
    v_icrs_matrix = galactic_to_icrs_vecs(v_gal)
    rot_err = float(np.max(np.abs(v_icrs_matrix - v_icrs_astropy)))
    assert rot_err < 1e-14, f"galactic_to_icrs_vecs error {rot_err}"
    print(f"\n- **Orthogonality Error (`max|R^T R - I|`):** `{orth_err:.3e}` (< 1e-14 PASS)")
    print(f"- **Astropy Vector Transform Error:** `{rot_err:.3e}` (< 1e-14 PASS)\n")


def validate_module_output(name, res):
    """Validates (fig, df, metrics), renders Matplotlib canvas cleanly, and closes figure."""
    fig, df, metrics = res
    assert isinstance(df, pd.DataFrame) and not df.empty, f"{name}: DataFrame is empty or invalid"
    num_df = df.select_dtypes(include=[np.number])
    if not num_df.empty:
        arr = num_df.to_numpy(dtype=float)
        assert np.all(np.isfinite(arr)), f"{name}: NaN or Inf found in numeric columns"

    assert isinstance(metrics, dict) and len(metrics) > 0, f"{name}: metrics dict empty"

    if isinstance(fig, plt.Figure):
        fig.canvas.draw()
        plt.close(fig)
    return df, metrics


def run_benchmark():
    verify_constants_and_frames()

    cat_label, d_max, sigma_default, tree, gal_pos, nodes, X, Y, Z, names = load_local_catalog()
    print(f"## 2. Headless 18-Module Benchmark ({cat_label}, N={len(gal_pos):,}, d_max={d_max:.0f} Mpc, sigma={sigma_default:.1f} Mpc)\n")

    d_depth = min(15.0, d_max)
    d_bao = min(round(R_BAO_MPC, 2), d_max)

    tasks = [
        ("1. Sky Anisotropy Maps", lambda: run_anisotropic_sky_forecast(
            tree, gal_pos, nodes, nside=16, d_max_mpc=d_depth, sigma_smooth_mpc=sigma_default, h_global=H_GLOBAL
        )),
        ("2. Line-of-Sight Kinematics", lambda: run_line_of_sight_kinematics_analysis(
            308.71, 60.15, 300.0, tree, gal_pos, nodes, d_known=7.72, sigma_mpc=sigma_default, h_global=H_GLOBAL, d_max_search=min(35.0, d_max)
        )),
        ("3. Redshift Path Integral", lambda: run_redshift_path_analysis(
            201.0, -29.8, 280.0, -20.0, tree, gal_pos, nodes, max_dist_mpc=d_max, sigma_mpc=sigma_default, h_global=H_GLOBAL
        )),
        ("4. Supernova Fit (Omega_L=0)", lambda: run_hubble_fitter_analysis(
            tree, gal_pos, nodes, sigma_mpc=sigma_default, h_global=H_GLOBAL
        )),
        ("5. KBC Supervoid Profile", lambda: run_kbc_void_analysis(
            r_max_mpc=500.0, delta_0_kbc=-0.22, h_local_peak=72.80, h_global=H_GLOBAL
        )),
        ("6. Bulk Flow Field", lambda: run_bulk_flow_analysis(
            tree, gal_pos, nodes, X, Y, Z, grid_size=30, spatial_bounds=min(25.0, d_max), sigma_mpc=sigma_default, h_global=H_GLOBAL
        )),
        ("7. Void Boundary Dynamics", lambda: run_void_boundary_analysis(
            tree, gal_pos, nodes, max_r_mpc=18.0, n_pts=200, h_global=H_GLOBAL
        )),
        ("8. Redshift Space Distortions", lambda: run_redshift_space_distortions(
            tree=tree, gal_positions=gal_pos, nodes_per_galaxy=nodes, alpha_inflow=120.0, sigma_mpc=sigma_default, h_global=H_GLOBAL
        )),
        ("9. Galactic Rotation Curves", lambda: run_rotation_curve_analysis(
            galaxy_target_name="NGC 3198 (Classic Spiral Benchmark)", h_global=H_GLOBAL, disk_geometry="Freeman Thin Disk"
        )),
        ("10. 3D Cosmic Web Viewer", lambda: run_3d_density_viewer_analysis(
            tree, gal_pos, nodes, names, sigma_mpc=sigma_default, h_global=H_GLOBAL, grid_res=18, max_galaxies=2000
        )),
        ("11. 3D Dark Energy Field", lambda: run_3d_dark_energy_analysis(
            tree, gal_pos, nodes, names, target_idx=0, max_dist_mpc=d_max, grid_n=40, sigma_mpc=sigma_default, h_global=H_GLOBAL
        )),
        ("12. BAO Ruler Warping", lambda: run_bao_warping_analysis(
            tree, gal_pos, nodes, path_depth_mpc=d_bao, sigma_mpc=sigma_default, h_global=H_GLOBAL
        )),
        ("13. CMB Dipole Anomaly", lambda: run_cmb_dipole_analysis(
            tree, gal_pos, nodes, nside=32, d_max_mpc=min(20.0, d_max), sigma_mpc=sigma_default, h_global=H_GLOBAL
        )),
        ("14. Gravitational Lensing", lambda: run_lensing_analysis(
            lens_target_name="Abell 1689 Core (Cluster Lens)", h_global=H_GLOBAL
        )),
        ("15. Spatial Garbage Collection", lambda: run_spatial_garbage_collection_analysis(
            preset_key="Earth", r_max_ratio=10.0
        )),
        ("16. CPU Bandwidth Throttling", lambda: run_cpu_bandwidth_throttling_analysis(
            preset_key="GPS Satellite Orbit (20,180 km Altitude)", observer_mode="Static (v = 0)"
        )),
        ("17. S8 Weak Lensing Tension", lambda: run_s8_tension_analysis(
            tree, gal_pos, nodes, names, sigma_halo_mpc=0.05,
            activation_mode="Smooth Volume Fraction (Coarse-Grained)", slice_extent_mpc=min(15.0, d_max)
        )),
        ("18. Topological Drag Floor (P_δ)", lambda: run_topological_drag_floor_analysis(
            n_sw_1au_cm3=7.0, v_sw_kms=440.0, planet_radius_km=6371.2, mw_disk_rd_kpc=2.85, mw_mass_solar=6.0e10
        )),
    ]

    print("| Module | Time (ms) | Rows x Cols | Status |")
    print("| :--- | ---: | :--- | :--- |")
    total_ms = 0.0
    captured = {}
    for label, fn in tasks:
        t0 = time.perf_counter()
        res = fn()
        dt_ms = (time.perf_counter() - t0) * 1000.0
        total_ms += dt_ms
        df, metrics = validate_module_output(label, res)
        captured[label] = (df, metrics)
        print(f"| {label} | {dt_ms:8.2f} | {len(df):,} x {df.shape[1]} | PASS |")

    print(f"| **TOTAL (All 18 Modules)** | **{total_ms:8.2f}** | — | **18/18 PASS** |\n")

    # ==========================================
    # 3. PHYSICS SPOT CHECKS & INVARIANT TESTS
    # ==========================================
    print("## 3. Physics Spot Checks & Regression Guards\n")

    _, m_freeman = validate_module_output(
        "NGC3198_Freeman", run_rotation_curve_analysis("NGC 3198 (Classic Spiral Benchmark)", h_global=H_GLOBAL, disk_geometry="Freeman Thin Disk")
    )
    _, m_spherical = validate_module_output(
        "NGC3198_Spherical", run_rotation_curve_analysis("NGC 3198 (Classic Spiral Benchmark)", h_global=H_GLOBAL, disk_geometry="Spherical Enclosed Mass")
    )
    _, m_ddo_f = validate_module_output(
        "DDO154_Freeman", run_rotation_curve_analysis("DDO 154 (Gas-Rich Dwarf Galaxy)", h_global=H_GLOBAL, disk_geometry="Freeman Thin Disk")
    )
    _, m_ddo_s = validate_module_output(
        "DDO154_Spherical", run_rotation_curve_analysis("DDO 154 (Gas-Rich Dwarf Galaxy)", h_global=H_GLOBAL, disk_geometry="Spherical Enclosed Mass")
    )
    assert m_ddo_f["Transition Radius (R*)"] == "Deep Underflow (< a_Ω)", m_ddo_f
    assert m_ddo_s["Transition Radius (R*)"] == "Deep Underflow (< a_Ω)", m_ddo_s

    _, m_earth_gc = validate_module_output("Earth_GC", run_spatial_garbage_collection_analysis("Earth"))
    _, m_gps = validate_module_output("GPS_Clock", run_cpu_bandwidth_throttling_analysis("GPS Satellite Orbit (20,180 km Altitude)"))
    _, m_iss = validate_module_output("ISS_Clock", run_cpu_bandwidth_throttling_analysis("ISS Orbit (400 km Altitude)"))
    _, m_sne_std = validate_module_output("SNe_0.315", run_hubble_fitter_analysis(h_global=H_GLOBAL, om=0.315))
    _, m_sne_gl = validate_module_output("SNe_Global", run_hubble_fitter_analysis(h_global=H_GLOBAL, om=OMEGA_M_GLOBAL))

    # CMB dipole check + catalog cache-key invalidation check across two distinct catalogs
    _, m_cmb_1 = validate_module_output(
        "CMB_Check_1", run_cmb_dipole_analysis(tree, gal_pos, nodes, nside=16, d_max_mpc=min(35.0, d_max), sigma_mpc=sigma_default, h_global=H_GLOBAL)
    )
    sub_n = max(10, len(gal_pos) // 4)
    sub_tree = cKDTree(gal_pos[:sub_n])
    _, m_cmb_2 = validate_module_output(
        "CMB_Check_2", run_cmb_dipole_analysis(sub_tree, gal_pos[:sub_n], nodes[:sub_n] * 2.5, nside=16, d_max_mpc=min(35.0, d_max), sigma_mpc=sigma_default, h_global=H_GLOBAL)
    )
    assert m_cmb_1["Substrate Anisotropy Range"] != m_cmb_2["Substrate Anisotropy Range"], (
        "CMB cache did not invalidate across distinct catalog_sig inputs!"
    )

    # Phase 2 spot checks: Lensing (1 - b), Garbage Collection R_quantum & Sgr A* Kernel Lock, S8 Smooth/Strict, P_delta 4-Regime
    _, m_lens = captured["14. Gravitational Lensing"]
    assert m_lens["SZ-Lensing Cluster Bias (1-b)"].startswith("0.8203"), m_lens

    assert m_earth_gc["Top-Down Quantum Radius"] == "6,371.20 km", m_earth_gc
    _, m_sgra_gc = validate_module_output("SgrA_GC", run_spatial_garbage_collection_analysis("Sagittarius A*"))
    assert "KERNEL LOCK" in m_sgra_gc["Kernel Lock Safety Margin"], m_sgra_gc

    _, m_s8_sm = captured["17. S8 Weak Lensing Tension"]
    _, m_s8_st = validate_module_output(
        "S8_Strict",
        run_s8_tension_analysis(
            tree, gal_pos, nodes, names, sigma_halo_mpc=0.05,
            activation_mode="Strict Heaviside Step Θ(ρ_N - 4.0)", slice_extent_mpc=min(15.0, d_max)
        ),
    )
    assert m_s8_sm["Void Floor (S8_void)"] == f"{S8_VOID:.4f}", m_s8_sm
    assert m_s8_sm["Global Baseline (S8_global)"] == f"{S8_GLOBAL:.4f}", m_s8_sm
    assert m_s8_sm["Cluster Core (S8_core)"] == f"{S8_CORE:.4f}", m_s8_sm
    assert 0.7541 <= float(m_s8_sm["Catalog Mean <S8>"]) <= 0.8327, m_s8_sm
    assert 0.7541 <= float(m_s8_st["Catalog Mean <S8>"]) <= 0.8327, m_s8_st
    assert "0.8203" in m_s8_sm["SZ-Lensing Bias (1 - b)"], m_s8_sm

    df18, m_pd = captured["18. Topological Drag Floor (P_δ)"]
    assert m_pd["MIT Bag Floor (B_QCD)"] == "47.62 MeV/fm³", m_pd
    assert m_pd["Derived Pion / Yukawa"] == "138.31 MeV | 1.427 fm", m_pd
    assert m_pd["Macro Drag Floor (P_δ)"] == "0.2834 pPa", m_pd
    assert m_pd["Shock / Heliopause"] == "89.44 / 120.74 AU", m_pd
    assert m_pd["Disk Cliff / Bubbles"] == "3.829 R_d | 8.55/11.40 kpc", m_pd
    r_gc_val = float(df18.loc[df18["System"] == "Circumgalactic Halo Cutoff (R_GC)", "Substrate_Predicted"].iloc[0])
    assert abs(r_gc_val - 525.55) <= 0.05, f"Unexpected R_GC: {r_gc_val}"

    print(f"- **NGC 3198 (Freeman):** `v_flat={m_freeman['Flat Asymptotic Velocity']}`, `R*={m_freeman['Transition Radius (R*)']}`, `Focus={m_freeman['Velocity @ Focus (r_focus)']}`, `R_trunc={m_freeman['Disk Truncation (3.829 R_d)']}`")
    print(f"- **NGC 3198 (Spherical):** `v_flat={m_spherical['Flat Asymptotic Velocity']}`, `R*={m_spherical['Transition Radius (R*)']}`, `Focus={m_spherical['Velocity @ Focus (r_focus)']}`, `R_trunc={m_spherical['Disk Truncation (3.829 R_d)']}`")
    print(f"- **DDO 154 Underflow Check:** Freeman `{m_ddo_f['Transition Radius (R*)']}` | Spherical `{m_ddo_s['Transition Radius (R*)']}` (PASS)")
    print(f"- **Earth GC:** `g_Omega={m_earth_gc['Surface Acceleration (g_Ω)']}` | Purge Rate: `{m_earth_gc['Surface Pixel Purge Rate']}` | `R_quantum={m_earth_gc['Top-Down Quantum Radius']}` | Sgr A*: `{m_sgra_gc['Kernel Lock Safety Margin']}`")
    print(f"- **Lensing Cluster Bias:** `1 - b = {m_lens['SZ-Lensing Cluster Bias (1-b)']}` (PASS)")
    print(f"- **S8 Weak Lensing ({cat_label}):** Smooth `<S8>={m_s8_sm['Catalog Mean <S8>']}` | Strict `<S8>={m_s8_st['Catalog Mean <S8>']}` | Tiers `{m_s8_sm['Void Floor (S8_void)']} / {m_s8_sm['Global Baseline (S8_global)']} / {m_s8_sm['Cluster Core (S8_core)']}`")
    print(f"- **P_delta 4-Regime Floor:** `B_QCD={m_pd['MIT Bag Floor (B_QCD)']}` | `Pion/Yukawa={m_pd['Derived Pion / Yukawa']}` | `TS/HP={m_pd['Shock / Heliopause']}` | `Disk/Bubbles={m_pd['Disk Cliff / Bubbles']}` | `R_GC={r_gc_val:.2f} kpc`")
    print(f"- **Atomic Clock Drift (vs Earth):** GPS `{m_gps['Relative Drift vs Earth']}` | ISS `{m_iss['Relative Drift vs Earth']}`")
    print(f"- **SNe Ia Agreement:** `om=0.315` -> `{m_sne_std['Substrate vs ΛCDM Agreement']}` | `OMEGA_M_GLOBAL={OMEGA_M_GLOBAL:.4f}` -> `{m_sne_gl['Substrate vs ΛCDM Agreement']}`")
    print(f"- **CMB Shifted Dipole ({cat_label}):** Axis `{m_cmb_1['Shifted Dipole Axis (l, b)']}` | Amp `{m_cmb_1['Combined Dipole Amplitude']}` | Sub-catalog cache switch PASS (`{m_cmb_1['Substrate Anisotropy Range']}` vs `{m_cmb_2['Substrate Anisotropy Range']}`)")


if __name__ == "__main__":
    run_benchmark()
