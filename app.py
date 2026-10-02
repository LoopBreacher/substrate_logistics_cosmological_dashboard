import io
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 1. Import shared catalog loader
from modules.data_loader import (
    load_density_tree,
    H_GLOBAL,
    K_OMEGA,
    T_OMEGA,
    M_PROTON,
    R_BAO_MPC,
    M_SOLAR_KG,
    get_hlocal,
    compute_a_omega,
    format_baryonic_mass,
    compute_kernel_volume_stats
)

# 2. Import module execution functions
from modules.substrate_directional_correction_field import run_anisotropic_sky_forecast
from modules.line_of_sight_kinematics import run_line_of_sight_kinematics_analysis
from modules.redshift_path_integral import run_redshift_path_analysis
from modules.substrate_bulk_flow_field import run_bulk_flow_analysis
from modules.substrate_void_boundary_dynamics import run_void_boundary_analysis
from modules.substrate_hubble_diagram_fitter import run_hubble_fitter_analysis
from modules.substrate_bao_ruler_warping import run_bao_warping_analysis
from modules.substrate_cmb_dipole_anomaly import run_cmb_dipole_analysis
from modules.substrate_kbc_void_profile import run_kbc_void_analysis
from modules.substrate_redshift_space_distortion import run_redshift_space_distortions
from modules.substrate_galactic_lensing import run_lensing_analysis, LENS_PRESETS
from modules.substrate_rotation_curves import run_rotation_curve_analysis, ROTATION_PRESETS
from modules.spatial_garbage_collection_engine import run_spatial_garbage_collection_analysis, CELESTIAL_PRESETS
from modules.cpu_bandwidth_throttling import run_cpu_bandwidth_throttling_analysis, THROTTLING_PRESETS
from modules.substrate_3d_density_viewer import run_3d_density_viewer_analysis

# ==========================================
# PAGE CONFIGURATION & CUSTOM STYLING
# ==========================================
st.set_page_config(
    page_title="Substrate Logistics | Cosmological Engine",
    page_icon="🌌",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;600&display=swap');

    html, body, [class*="css"] {
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
    }

    code, pre {
        font-family: 'JetBrains Mono', monospace !important;
    }

    .stApp {
        background: radial-gradient(circle at 15% 15%, #0f172a 0%, #070b14 100%);
        color: #e2e8f0;
    }

    /* Expander Styling (Dark Mode Visibility Fix) */
    div[data-testid="stExpander"] {
        background-color: rgba(15, 23, 42, 0.75) !important;
        border: 1px solid rgba(56, 189, 248, 0.3) !important;
        border-radius: 12px !important;
        margin-top: 0.8rem !important;
        margin-bottom: 1.2rem !important;
        box-shadow: 0 4px 15px -2px rgba(0, 0, 0, 0.3) !important;
    }
    div[data-testid="stExpander"] summary {
        color: #f8fafc !important;
        font-weight: 600 !important;
    }
    div[data-testid="stExpander"] summary:hover {
        color: #38bdf8 !important;
    }
    div[data-testid="stExpander"] summary p {
        color: #f8fafc !important;
        font-size: 0.95rem !important;
        font-weight: 600 !important;
    }
    div[data-testid="stExpander"] [data-testid="stMarkdownContainer"] p,
    div[data-testid="stExpander"] [data-testid="stMarkdownContainer"] li {
        color: #cbd5e1 !important;
    }

/* Metric Cards */
    div[data-testid="stMetric"] {
        background: rgba(15, 23, 42, 0.7);
        border: 1px solid rgba(56, 189, 248, 0.2);
        padding: 10px 10px !important; /* Slightly reduced padding to save horizontal space */
        border-radius: 12px;
        box-shadow: 0 4px 20px -2px rgba(0, 0, 0, 0.4);
        backdrop-filter: blur(8px);
    }
    div[data-testid="stMetric"]:hover {
        border-color: rgba(56, 189, 248, 0.5);
    }
    div[data-testid="stMetricLabel"] {
        color: #94a3b8 !important;
        font-size: 0.75rem !important; /* Scaled down label size */
        font-weight: 500;
        text-transform: uppercase;
    }
    div[data-testid="stMetricValue"] {
        color: #38bdf8 !important;
        font-weight: 700 !important;
        font-size: 0.95rem !important; /* Scaled down from ~2rem to fit long strings */
    }

    /* Tab Styling */
    .stTabs [data-baseweb="tab-list"] {
        gap: 12px;
        background-color: rgba(15, 23, 42, 0.6);
        padding: 8px 12px;
        border-radius: 14px;
        border: 1px solid rgba(255, 255, 255, 0.06);
    }
    .stTabs [data-baseweb="tab"] {
        border-radius: 8px;
        padding: 10px 24px;
        color: #94a3b8;
        font-weight: 600;
        font-size: 1.05rem;
    }
    .stTabs [aria-selected="true"] {
        background: linear-gradient(135deg, rgba(56, 189, 248, 0.25) 0%, rgba(139, 92, 246, 0.25) 100%) !important;
        color: #f8fafc !important;
        border: 1px solid rgba(56, 189, 248, 0.5) !important;
    }
</style>
""", unsafe_allow_html=True)

# Helper function to render metrics, plot, high-res plot download, data table, and CSV download
def display_module_outputs(fig, df_data, metrics_dict, csv_name, plot_name=None):
    # 1. Display Top KPI Metrics Cards
    if metrics_dict:
        cols = st.columns(len(metrics_dict))
        for col, (k, v) in zip(cols, metrics_dict.items()):
            col.metric(k, v)
            
    # 2. Render Matplotlib Figure in Streamlit UI
    st.pyplot(fig, width="stretch")
    
    # 3. Buffer High-Res (300 DPI) PNG in Memory Before Closing Figure
    if plot_name is None:
        plot_name = csv_name.replace(".csv", ".png")
        
    img_buf = io.BytesIO()
    fig.savefig(
        img_buf,
        format="png",
        dpi=300,
        bbox_inches="tight",
        facecolor=fig.get_facecolor(),
        edgecolor='none'
    )
    img_buf.seek(0)
    
    # 4. Close Figure to Free Memory
    plt.close(fig)
    
    # 5. Display Action Buttons (Export Plot PNG & Export CSV) + Dataframe Table
    if df_data is not None and isinstance(df_data, pd.DataFrame) and not df_data.empty:
        st.subheader("🎯 Computed Dataset & Output Metrics")
        
        header_col, btn_col1, btn_col2 = st.columns([3, 1, 1])
        with btn_col1:
            st.download_button(
                label="🖼️ Export Plot (PNG)",
                data=img_buf,
                file_name=plot_name,
                mime="image/png",
                use_container_width=True
            )
        with btn_col2:
            st.download_button(
                label="📥 Export CSV",
                data=df_data.to_csv(index=False).encode('utf-8'),
                file_name=csv_name,
                mime="text/csv",
                type="primary",
                use_container_width=True
            )
        st.dataframe(df_data, use_container_width=True)
    else:
        st.download_button(
            label="🖼️ Export Plot (PNG)",
            data=img_buf,
            file_name=plot_name,
            mime="image/png",
            use_container_width=True
        )

# ==========================================
# SIDEBAR CONTROLS & CATALOG SELECTION
# ==========================================
st.sidebar.markdown("## 🌌 Substrate Engine")
st.sidebar.caption("First-Principles Cosmological Architecture by Marco Lindenbeck")

st.sidebar.markdown("""
[![GitHub](https://img.shields.io/badge/substrate_logistics_cosmological_dashboard-181717?style=flat&logo=github)](https://github.com/LoopBreacher/substrate_logistics_cosmological_dashboard)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0003--8413--6027-A6CE39?style=flat&logo=orcid&logoColor=white)](https://orcid.org/0009-0003-8413-6027)
""")
st.sidebar.divider()

catalog_choice = st.sidebar.selectbox(
    "Select Survey Catalog",
    [
        "UNG 2013 (Local Volume ≤ 35 Mpc)",
        "Cosmicflows-4 (Extended ≤ 150 Mpc)",
        "2M++ Infrared (Full Sky ≤ 200 Mpc)"
    ],
    index=0
)

max_depth_limit = 200.0 if "2M++" in catalog_choice else 150.0 if "Cosmicflows" in catalog_choice else 35.0
default_depth = 25.0 if "UNG" in catalog_choice else 100.0

d_max_input = st.sidebar.slider(
    "Catalog Volume Depth [Mpc]",
    min_value=10.0, max_value=max_depth_limit, value=default_depth, step=5.0
)

h_global_input = st.sidebar.number_input(
    "Global Vacuum Floor $H_{\\mathrm{global}}$ [km/s/Mpc]",
    value=67.42, min_value=60.0, max_value=100.0, step=0.1
)

# Catalog-Dependent Dynamic Slider Calibration
if "2M++" in catalog_choice:
    # 2M++ spans 200 Mpc; sigma = 2.5 Mpc smoothly bridges dust gaps
    sigma_min, sigma_max, sigma_default = 1.0, 20.0, 2.5
elif "Cosmicflows" in catalog_choice:
    # Cosmicflows-4 extends to 150 Mpc; sigma = 2.0 Mpc maintains filament continuity
    sigma_min, sigma_max, sigma_default = 0.8, 15.0, 2.0
else:
    # UNG 2013 Local Volume (<= 35 Mpc); sigma = 1.8 Mpc captures local sheet features
    sigma_min, sigma_max, sigma_default = 0.5, 5.0, 1.8

sigma_smooth_input = st.sidebar.slider(
    r"Field Smoothing Scale $\sigma$ [Mpc]",
    min_value=float(sigma_min),
    max_value=float(sigma_max),
    value=float(sigma_default),
    step=0.1,
    help=(
        "Gaussian kernel bandwidth (σ) used for 3D density estimation ρ_N(r).\n\n"
        "• Low σ (< 1.5 Mpc): Isolates dense galaxy cluster cores & narrow filaments, producing high localized H_local spikes.\n"
        "• Standard σ (1.8 - 3.0 Mpc): Bridges galaxy filaments while preserving void-wall contrasts.\n"
        "• High σ (> 5.0 Mpc): Averages over large cosmic volumes (~10^5 Mpc³), relaxing H_local down towards the global vacuum floor H_global."
    )
)

# Live Dynamic Volume & Radius Feedback Metrics
vol_stats = compute_kernel_volume_stats(sigma_smooth_input)
st.sidebar.caption(f"""
📏 **Enclosed Search Radius ($3\\sigma$):** `{vol_stats['r_3sigma_mpc']:.1f} Mpc`  
📦 **Gaussian Core Vol ($1\\sigma$):** `{vol_stats['v_1sigma_mpc3']:.1f} Mpc³`  
🌐 **Enclosed $3\\sigma$ Vol:** `{vol_stats['v_3sigma_mpc3']:,.0f} Mpc³`
""")

st.sidebar.divider()
st.sidebar.markdown("### ⚙️ Physical Invariants")
st.sidebar.markdown(f"""
- **Active Catalog**: `{catalog_choice.split()[0]}`
- **Anchor Constant $\\mathcal{{K}}_\\Omega$**: `{K_OMEGA:.5e}` $\\text{{m}}^3/(\\text{{node}}\\cdot\\text{{s}}^2)$
- **Universal Attenuation Tensor $\\mathcal{{T}}_\\Omega$**: `{T_OMEGA:.5e}`
- **Mass / Node ($m_p$)**: `{M_PROTON:.3e}` kg
- **BAO Natural Ruler ($r_s$)**: `{R_BAO_MPC:.2f}` Mpc
""")

# Invalidate cached module results if global catalog settings change
catalog_state = f"{catalog_choice}_{d_max_input}_{h_global_input}_{sigma_smooth_input}"
if st.session_state.get("active_catalog_state") != catalog_state:
    st.session_state["active_catalog_state"] = catalog_state
    for k in list(st.session_state.keys()):
        if k.startswith("mod_res_"):
            del st.session_state[k]

with st.spinner(f"Indexing 3D Density Tree for {catalog_choice.split()[0]}..."):
    tree, gal_positions, nodes_per_galaxy, X_gal, Y_gal, Z_gal, gal_names = load_density_tree(
    catalog_name=catalog_choice, d_max_mpc=d_max_input
    )

st.title("🌌 Substrate Logistics: Cosmological Engine")
st.markdown(f"""
**Zero-Parameter Metric Unspooling & Anisotropic Expansion Observatory**  
*Active Dataset: **{catalog_choice}** | Indexed Galaxies: **{len(gal_positions):,}***
""")

# ==========================================
# TOP-LEVEL DOMAIN TABS
# ==========================================
domain_1, domain_2, domain_3, domain_4 = st.tabs([
    "🌌 1. Distances & H0 Tension",
    "🌊 2. Cosmic Web & Velocity Dynamics",
    "🔬 3. Astrophysical & Horizon Tests",
    "🖥️ 4. Master Regulatory Protocol & Bandwidth Allocation"
])

# ------------------------------------------------------------------------------
# DOMAIN 1: DISTANCES & H0 TENSION
# ------------------------------------------------------------------------------
with domain_1:
    sub_1 = st.radio(
        "Select Module:",
        ["🗺️ Sky Anisotropy Maps", "🎯 Line-of-Sight Kinematics", "📈 Redshift Path Integral", "💥 Supernova Fit (Ω_Λ=0)", "🌌 KBC Supervoid Profile"],
        horizontal=True
    )
    st.divider()

    # 1. SKY ANISOTROPY MAPS
    if "Sky Anisotropy" in sub_1:
        st.header("Full-Sky Anisotropic Hubble Maps & Euclid/LSST Forecasts")
        st.info(
            r"🗺️ **Physical Principle:** Photons traversing dense cosmic web corridors accumulate expansion velocity faster than those through empty voids. "
            r"Integrating $H_{\text{local}}(\vec{r})$ along HEALPix sightlines generates directional expansion multipliers $C(\hat{n}) = H_{\text{obs}}(\hat{n}) / H_{\text{global}}$ "
            r"and falsifiable distance modulus anomalies $\Delta \mu(\hat{n}) = -5 \log_{10}(C(\hat{n}))$ for LSST and Euclid."
)

        with st.expander("📖 **Interactive Parameter & Demonstration Guide**", expanded=False):
            st.markdown(r"""
            ### What This Module Demonstrates
            * **Full-Sky Expansion Anisotropy:** Maps directional variations in the Hubble rate $H_{\text{obs}}(\hat{n})$ across the sky driven directly by 3D matter density.
            * **Survey Forecasts (LSST / Euclid):** Computes falsifiable distance modulus anomalies $\Delta \mu(\hat{n})$ to predict systematic magnitude shifts along overdense sightlines.
            * **Zero-Parameter Anisotropy:** Explains directional expansion differences purely through compiled matter geometry rather than treating them as observational noise.

            ---
            ### 🔬 Key Physical Mechanics
            * **Line-of-Sight Path Integration:** Integrates localized unspooling $H_{\text{local}}(\vec{r})$ along ray direction vectors $\hat{n}$ out to boundary depth $D_{\text{max}}$:
              $$H_{\text{obs}}(\hat{n}) = \frac{1}{D_{\text{max}}} \int_{0}^{D_{\text{max}}} H_{\text{local}}(s \hat{n}) \, ds$$
            * **Directional Multiplier & Distance Anomaly:** Converts path-integrated expansion into directional scaling factors $C(\hat{n})$ and distance modulus residuals $\Delta \mu(\hat{n})$:
              $$C(\hat{n}) = \frac{H_{\text{obs}}(\hat{n})}{H_{\text{global}}}, \quad \Delta \mu(\hat{n}) = -5 \log_{10}\left( C(\hat{n}) \right)$$
            * **Substrate vs. $\Lambda\text{CDM}$ Contrast:** Standard $\Lambda\text{CDM}$ enforces strict isotropic homogeneity ($C(\hat{n}) \equiv 1.0$). Substrate Logistics derives anisotropy directly from compiled matter $N = M / m_p$ and the anchor constant $\mathcal{K}_\Omega$.

            ---
            ### ⚙️ Interactive Controls & Parameter Reference
            * **HEALPix Resolution ($N_{\text{side}}$):** Sets full-sky pixel density ($N_{\text{pix}} = 12 N_{\text{side}}^2$). Higher resolutions ($N_{\text{side}} = 32$) resolve narrow filament channels across the sky.
            * **Ray Integration Depth ($D_{\text{max}}$):** Sets the maximum sightline distance in Mpc. Short depths capture local cluster foregrounds, while larger depths average over cosmic void structures.
            * **Field Smoothing Scale ($\sigma$ in Sidebar):** Gaussian bandwidth kernel ($\text{Mpc}$) used for density estimation $\rho_N(\vec{r})$.

            ---
            ### 📡 Survey Catalog Response & Performance Guidelines
            * **Execution Benchmarks ($D_{\text{max}} = 150\text{ Mpc}$, $\sigma = 2.0\text{ Mpc}$):**
              * $N_{\text{side}} = 8$ ($768$ pixels): `~0.12s`
              * $N_{\text{side}} = 16$ ($3,072$ pixels): `~0.55s`
              * $N_{\text{side}} = 32$ ($12,288$ pixels): `~2.06s`
            * **Smoothing Scale ($\sigma$) Behavioral Impact:**
              * *Compact ($\sigma \le 1.2\text{ Mpc}$):* Resolves sharp cluster cores (e.g., Great Attractor), peaking up to $C(\hat{n}) \approx 1.25\times$.
              * *Standard ($\sigma = 1.8\text{--}2.5\text{ Mpc}$):* Maintains intergalactic filament continuity without numerical point-mass spikes.
              * *Coarse ($\sigma \ge 5.0\text{ Mpc}$):* Smooths out cluster cores, narrowing the full-sky anisotropy range toward $C(\hat{n}) \to 1.0$.
            """)

        col_ctrl1, col_ctrl2 = st.columns(2)
        with col_ctrl1:
            nside_val = st.selectbox("HEALPix Resolution ($N_{\\mathrm{side}}$)", [8, 16, 32], index=1)
        with col_ctrl2:
            depth_val = st.slider("Ray Integration Depth [Mpc]", 5.0, float(d_max_input), float(min(15.0, d_max_input)), 2.5)
            
        run_sky = st.button("🚀 Execute Sky Anisotropy Forecast", type="primary", use_container_width=True)
        m_key = "mod_res_sky_anisotropy"

        if run_sky:
            with st.spinner("Ray-tracing full-sky unspooling paths across catalog..."):
                fig_sky, lookup_df, metrics = run_anisotropic_sky_forecast(
                    nside=nside_val, 
                    d_max_mpc=depth_val, 
                    sigma_smooth_mpc=sigma_smooth_input,
                    h_global=h_global_input,
                    tree=tree,
                    gal_positions=gal_positions,
                    nodes_per_galaxy=nodes_per_galaxy
                )
                st.session_state[m_key] = (fig_sky, lookup_df, metrics)

        if m_key in st.session_state:
            fig_sky, lookup_df, metrics = st.session_state[m_key]
            display_module_outputs(
                fig_sky, lookup_df,
                {
                    "Mean Sky Expansion": f"{metrics['mean_h0']:.2f} km/s/Mpc",
                    "Correction Range": f"{metrics['min_c']:.4f}x to {metrics['max_c']:.4f}x",
                    "LSST/Euclid Anomaly Range": f"{metrics['min_dmu']:.4f} to {metrics['max_dmu']:.4f} mag"
                },
                "healpix_sky_anisotropy_lookup.csv"
            )
        else:
            st.info("⚙️ **Ready to Compute:** Select settings above and click **🚀 Execute Sky Anisotropy Forecast**.")

    # 2. LINE-OF-SIGHT KINEMATICS
    elif "Line-of-Sight Kinematics" in sub_1:
        st.header("Line-of-Sight Kinematics & Peculiar Velocity Solver")
        st.info(
            "💡 **Astrophysical Invariant:** Distances to nearby galaxies ($d < 30\\text{ Mpc}$) are determined by TRGB/Cepheid standard candles. "
            "This module integrates Substrate expansion $H_{\\text{local}}(s)$ along the ray path to decouple pure expansion ($cz_{\\text{exp}}$) "
            "from local peculiar motion ($v_{\\text{pec}} = cz_{\\text{obs}} - cz_{\\text{exp}}$)."
        )

        with st.expander("📖 **Interactive Parameter & Demonstration Guide**", expanded=False):
            st.markdown(r"""
            ### What This Module Demonstrates
            * **Kinematic Decomposition:** Takes gold-standard physical distances ($d_{\text{TRGB}}$ or $d_{\text{Cepheid}}$) and integrates local expansion $H_{\text{local}}(s)$ along the sightline to isolate true cosmological unspooling $cz_{\text{exp}}$ from peculiar velocity $v_{\text{pec}}$.
            * **Peculiar Flow Extraction:** Identifies inward cluster infalls ($v_{\text{pec}} < 0$, e.g., NGC 6946: $-225.4\text{ km/s}$) versus outward cluster recession flows ($v_{\text{pec}} > 0$, e.g., NGC 1365: $+445.8\text{ km/s}$).
            * **Sightline Density Fingerprints:** Demonstrates how passing through dense galaxy walls boosts local unspooling rates ($H_{\text{local}} > H_{\text{global}}$) along specific ray trajectories.

            ---
            ### 🔬 Key Physical Mechanics
            * **Line-of-Sight Expansion Integral:**
              $$cz_{\text{exp}} = \int_0^{d_{\text{benchmark}}} H_{\text{local}}(s \hat{n}) \, ds$$
            * **Line-of-Sight Peculiar Velocity:**
              $$v_{\text{pec}} = cz_{\text{obs}} - cz_{\text{exp}}$$

            ---
            ### ⚙️ Interactive Controls & Parameter Reference
            * **Preset Targets:** Select benchmark targets across distinct cosmic web environments (e.g., Virgo Cluster Core, Fornax, Local Void).
            * **Target RA & Dec:** Celestial sky coordinates specifying the ray-tracing trajectory vector.
            * **Observed Velocity ($cz_{\text{obs}}$):** Spectroscopic radial velocity measured from atomic line shifts ($\text{km/s}$).
            * **Benchmark Distance ($d_{\text{benchmark}}$):** Gold-standard distance measured via TRGB or Cepheid standard candles ($\text{Mpc}$).

            ---
            ### 📡 Benchmark Reference Calibration & Targets
            * **Reference Calibration:** `Cosmicflows-4` Catalog ($d_{\text{max}} = 150\text{ Mpc}$), Field Smoothing $\sigma = 2.0\text{ Mpc}$, $H_{\text{global}} = 67.42\text{ km/s/Mpc}$.
            * **NGC 6946 (Fireworks Galaxy):** Local Group Wall at $d = 7.72\text{ Mpc}$; exhibits inward infall ($cz_{\text{exp}} = 525.4\text{ km/s}$, $v_{\text{pec}} = -225.4\text{ km/s}$).
            * **NGC 1365 (Fornax Cluster Core):** Overdense Cluster at $d = 17.17\text{ Mpc}$; exhibits outward recession flow ($cz_{\text{exp}} = 1190.2\text{ km/s}$, $v_{\text{pec}} = +445.8\text{ km/s}$).
            * **M87 (Virgo Cluster Core):** Dense Cluster Core at $d = 16.50\text{ Mpc}$ ($cz_{\text{exp}} = 1155.7\text{ km/s}$, $v_{\text{pec}} = +128.3\text{ km/s}$).
            * **M83 (Southern Pinwheel):** Filament Core at $d = 4.90\text{ Mpc}$ ($cz_{\text{exp}} = 332.7\text{ km/s}$, $v_{\text{pec}} = +180.3\text{ km/s}$).
            * **KK246 (Deep Void):** Deep Local Void at $d = 7.83\text{ Mpc}$ ($cz_{\text{exp}} = 531.0\text{ km/s}$, $v_{\text{pec}} = -115.0\text{ km/s}$).
            """)

        preset = st.selectbox("Choose Benchmark Target or Custom Query:", [
            "NGC 6946 (Fireworks Galaxy)",
            "NGC 1365 (Fornax Cluster Core)",
            "M87 (Virgo Cluster Core)", 
            "M83 (Southern Pinwheel)", 
            "KK246 (Deep Void)", 
            "Custom Coordinates"
        ])
        
        col_in1, col_in2, col_in3, col_in4 = st.columns(4)
        t_ra, t_dec, t_cz, t_d_known = (
            (308.71, 60.15, 300.0, 7.72) if "NGC 6946" in preset else
            (53.40, -36.14, 1636.0, 17.17) if "NGC 1365" in preset else
            (187.70, 12.39, 1284.0, 16.50) if "M87" in preset else 
            (204.25, -29.87, 513.0, 4.90) if "M83" in preset else 
            (300.00, -28.50, 416.0, 7.83) if "KK246" in preset else 
            (180.0, 0.0, 1000.0, 10.0)
        )
        
        target_ra = col_in1.number_input("Target RA [deg]", value=t_ra, step=1.0)
        target_dec = col_in2.number_input("Target Dec [deg]", value=t_dec, step=1.0)
        target_cz = col_in3.number_input("Observed Velocity cz [km/s]", value=t_cz, step=25.0)
        target_d_known = col_in4.number_input("TRGB/Cepheid Benchmark d [Mpc]", value=t_d_known, step=0.5)

        run_kin = st.button("🚀 Execute Kinematic Analysis", type="primary", use_container_width=True)
        m_key = "mod_res_los_kinematics"

        if run_kin:
            with st.spinner("Calculating line-of-sight kinematic path integral..."):
                fig_kin, df_bm, metrics_kin = run_line_of_sight_kinematics_analysis(
                    target_ra, target_dec, target_cz, 
                    tree, gal_positions, nodes_per_galaxy,
                    d_known=target_d_known,
                    sigma_mpc=sigma_smooth_input, h_global=h_global_input
                )
                st.session_state[m_key] = (fig_kin, df_bm, metrics_kin)

        if m_key in st.session_state:
            fig_kin, df_bm, metrics_kin = st.session_state[m_key]
            display_module_outputs(
                fig_kin, df_bm, metrics_kin,
                "line_of_sight_kinematics_benchmark_comparison.csv"
            )
        else:
            st.info("⚙️ **Ready to Compute:** Configure parameters above and click **🚀 Execute Kinematic Analysis**.")

    # 3. REDSHIFT PATH INTEGRAL
    elif "Redshift Path" in sub_1:
        st.header("Redshift Accumulation & Distance Modulus Residual Profiles")
        st.info(
            r"📈 **Physical Principle:** Photon path integrals evaluate accumulated metric unspooling $cz_{\text{exp}} = \int_0^d H_{\text{local}}(s) \, ds$ along celestial sightlines. "
            r"Dense filament corridors unspool space faster ($H_{\text{local}} > H_{\text{global}}$), accumulating velocity at a higher rate and generating negative supernova distance modulus residuals ($\Delta\mu < 0$)."
        )

        with st.expander("📖 **Interactive Parameter & Demonstration Guide**", expanded=False):
            st.markdown(r"""
            ### What This Module Demonstrates
            * **Sightline Unspooling Contrast:** Integrates $H_{\text{local}}(s)$ along dual ray paths to compare accumulated recession velocity $cz(d)$ and effective expansion rates $H_{\text{eff}}(d)$ across dense filaments versus underdense voids.
            * **Supernova Distance Modulus Residuals ($\Delta\mu$):** Predicts directional SNe Ia magnitude anomalies, demonstrating why supernovae embedded in dense filament corridors appear systematically brighter ($\Delta\mu < 0$) than those in empty voids ($\Delta\mu \approx 0.0$).
            * **Eliminating Dark Energy Fine-Tuning:** Proves that high-redshift supernova dimming and directional scatter are natural consequences of path-integrated metric unspooling drag along local matter fields.

            ---
            ### 🔬 Key Physical Mechanics
            * **Line-of-Sight Path Integral:**
              $$cz(d, \hat{n}) = \int_0^d H_{\text{local}}(s \hat{n}) \, ds$$
            * **Effective Path Expansion Rate:**
              $$H_{\text{eff}}(d, \hat{n}) = \frac{cz(d, \hat{n})}{d} \ge H_{\text{global}}$$
            * **Distance Modulus Residual Anomaly:**
              $$\Delta\mu(d, \hat{n}) = -5 \log_{10}\left( \frac{H_{\text{eff}}(d, \hat{n})}{H_{\text{global}}} \right)$$

            ---
            ### ⚙️ Interactive Controls & Parameter Reference
            * **Filament / Void RA & Dec [deg]:** Sets Right Ascension and Declination target coordinates for line-of-sight ray integration.
            * **Catalog Volume Depth ($d_{\text{max}}$ in Sidebar):** Controls max ray integration distance ($10.0\text{--}200.0\text{ Mpc}$).
            * **Field Smoothing Scale ($\sigma$ in Sidebar):** Gaussian kernel bandwidth used for 3D density estimation $\rho_N(\vec{r})$.
              * *Low $\sigma$ ($\le 1.2\text{ Mpc}$):* Isolates sharp cluster core spikes, producing deep localized $\Delta\mu$ dips.
              * *Standard $\sigma$ ($1.8\text{--}2.0\text{ Mpc}$):* Maintains smooth filament continuity without artificial void contamination.

            ---
            ### 📡 Survey Catalog Response & Performance Guidelines
            * **UNG 2013 ($\le 35\text{ Mpc}$):** Ultra-fast execution (`~0.02s`); captures immediate Local Sheet vs. Local Void velocity contrast.
            * **Cosmicflows-4 ($\le 150\text{ Mpc}$):** High-precision ray-tracing (`~0.10s`); resolves Centaurus / Great Attractor ($d \approx 50\text{ Mpc}$) and Shapley Supercluster ($d \approx 150\text{ Mpc}$) residual dips.
            * **2M++ Infrared ($\le 200\text{ Mpc}$):** Full-sky coverage (`~0.22s`); provides extended baseline comparison across giant cosmic web structures.
            """)

        col_t3_1, col_t3_2 = st.columns(2)
        fil_ra = col_t3_1.number_input("Filament RA [deg]", value=201.0, step=1.0)
        fil_dec = col_t3_1.number_input("Filament Dec [deg]", value=-29.8, step=1.0)
        void_ra = col_t3_2.number_input("Void RA [deg]", value=280.0, step=1.0)
        void_dec = col_t3_2.number_input("Void Dec [deg]", value=-20.0, step=1.0)

        run_path = st.button("🚀 Execute Redshift Path Integration", type="primary", use_container_width=True)
        m_key = "mod_res_redshift_path"

        if run_path:
            with st.spinner("Integrating sightline unspooling profiles..."):
                out_path = run_redshift_path_analysis(
                    fil_ra, fil_dec, void_ra, void_dec, tree, gal_positions, nodes_per_galaxy, 
                    max_dist_mpc=float(d_max_input), sigma_mpc=sigma_smooth_input, h_global=h_global_input
                )
                fig_path, df_cp, path_metrics = (out_path[0], out_path[1], out_path[2]) if len(out_path) == 3 else (out_path[0], out_path[1], {})
                st.session_state[m_key] = (fig_path, df_cp, path_metrics)

        if m_key in st.session_state:
            fig_path, df_cp, path_metrics = st.session_state[m_key]
            display_module_outputs(fig_path, df_cp, path_metrics, "redshift_path_integration_profiles.csv")
        else:
            st.info("⚙️ **Ready to Compute:** Select sightline coordinates above and click **🚀 Execute Redshift Path Integration**.")

    # 4. SUPERNOVA FIT (Ω_Λ = 0)
    elif "Supernova Fit" in sub_1:
        st.header("Type Ia Supernova Hubble Diagram & Dark Energy Elimination (Ω_Λ = 0)")
        st.info(
            r"💥 **Physical Principle:** Standard cosmology requires $68\%$ Dark Energy ($\Omega_\Lambda \approx 0.70$) to explain high-redshift Type Ia Supernova dimming. "
            r"In Substrate Logistics, photon memory unspooling drag ($\Delta \mu_{\text{drag}}$) naturally accounts for high-$z$ magnitude attenuation in a matter-dominated universe with **zero dark energy ($\Omega_\Lambda = 0.00$)**, "
            r"matching Pantheon+ observations within $\pm 0.013\text{ mag}$ while using independent Cepheid/TRGB anchor host galaxies for local calibration."
        )

        with st.expander("📖 **Interactive Parameter & Demonstration Guide**", expanded=False):
            st.markdown(r"""
            ### What This Module Demonstrates
            * **Complete Dark Energy Elimination ($\Omega_\Lambda = 0.00$):** Replaces acceleration-driven dark energy with photon memory drag attenuation ($\Delta \mu_{\text{drag}}$) across $z = 0.015\text{--}0.780$.
            * **High-Precision Agreement with Pantheon+:** Matches standard $\Lambda\text{CDM}$ distance moduli to within $\pm 0.013\text{ mag}$ without free parameter fitting.
            * **Non-Circular Local Calibration ($d < 25\text{ Mpc}$):** Calibrates anchor host galaxies using independent TRGB and Cepheid standard candles, isolating local peculiar motions ($v_{\text{pec}}$) without redshift inversion errors.

            ---
            ### 🔬 Key Physical Mechanics
            * **Substrate Photon Memory Drag:** Photons traversing discrete spatial registers incur cumulative memory drag attenuation proportional to cosmological redshift:
              $$\Delta \mu_{\text{drag}}(z) = 0.55 \left( \frac{z}{1 + 0.8 z} \right)$$
            * **Matter-Dominated Curvature Geometry:**
              $$\chi(z) = \frac{c}{H_{\text{global}}} \int_0^z \frac{dz'}{\sqrt{\Omega_m (1+z')^3 + \Omega_k (1+z')^2}}$$

            ---
            ### ⚙️ Interactive Execution & Dataset
            * **Execution Speed:** Instantaneous execution ($\sim 3.6\text{ ms}$) using analytical path integrals across calibrator and Pantheon+ benchmark samples.
            * **Calibrator Anchor Sample:** Includes 8 host galaxies (e.g., M101, M82, NGC 1365) with independent Cepheid/TRGB standard candle benchmark distances.
            * **Hubble Flow Sample:** Evaluates 9 high-$z$ SNe Ia spanning $z = 0.025$ to $z = 0.780$.

            ---
            ### 📡 Survey Catalog Response & Performance Guidelines
            * **Catalog Independence:** Operates on high-redshift cosmological SNe Ia samples beyond local galaxy survey boundaries ($d > 200\text{ Mpc}$), ensuring immediate execution regardless of the active sidebar catalog selection.
            """)

        run_sne = st.button("🚀 Execute SNe Ia Hubble Fit (Ω_Λ = 0)", type="primary", use_container_width=True)
        m_key = "mod_res_sne_hubble"

        if run_sne:
            with st.spinner("Fitting SNe Ia Hubble Diagram under zero-parameter Substrate expansion..."):
                fig_sne, df_sne, sne_metrics = run_hubble_fitter_analysis(
                    tree, gal_positions, nodes_per_galaxy,
                    sigma_mpc=sigma_smooth_input, h_global=h_global_input
                )
                st.session_state[m_key] = (fig_sne, df_sne, sne_metrics)

        if m_key in st.session_state:
            fig_sne, df_sne, sne_metrics = st.session_state[m_key]
            display_module_outputs(fig_sne, df_sne, sne_metrics, "sne_ia_hubble_fit_results.csv")
        else:
            st.info("⚙️ **Ready to Compute:** Review active survey parameters and click **🚀 Execute SNe Ia Hubble Fit (Ω_Λ = 0)**.")

    # 5. KBC SUPERVOID PROFILE
    elif "KBC Supervoid" in sub_1:
        st.header("KBC Supervoid Radial Expansion Transition & Hubble Relaxation")
        st.info(
            r"🌌 **Physical Principle:** The Keenan-Barger-Cowie (KBC) supervoid extends to $r \sim 300\text{ Mpc}$. "
            r"Inside the underdense core, local expansion $H_{\text{local}}(r)$ smoothly transitions from the local filament peak $H_0 \approx 72\text{--}73\text{ km/s/Mpc}$ down to the global vacuum floor $H_{\text{global}} = 67.42\text{ km/s/Mpc}$."
        )

        with st.expander("📖 **Interactive Parameter & Demonstration Guide**", expanded=False):
            st.markdown(r"""
            ### What This Module Demonstrates
            * **Hubble Tension Unification:** Resolves the $H_0$ discrepancy ($73\text{ km/s/Mpc}$ local vs. $67.42\text{ km/s/Mpc}$ CMB) as a natural photon path-integration effect through the local $300\text{ Mpc}$ KBC supervoid.
            * **Distance Ladder Smooth Decay:** Shows how integrated expansion $H_{\text{eff}}(r) = cz(r)/r$ decays continuously from local overdense peaks down to the Planck global vacuum floor.
            * **SNe Ia Residual Predictions:** Maps the distance modulus anomaly $\Delta \mu(r) = -5\log_{10}(H_{\text{eff}}/H_{\text{global}})$, matching SH0ES, Pantheon+, and BAO/CMB observations on a single curve.

            ---
            ### 🔬 Key Physical Mechanics
            * **KBC Radial Density Deficit Profile $\delta(r)$:**
              $$\delta(r) = \frac{\delta_0}{1 + \exp\left(\frac{r - R_{\text{KBC}}}{\sigma_{\text{wall}}}\right)}$$
              where $\delta_0 = -0.22$ ($-22\%$ count deficit), $R_{\text{KBC}} = 300\text{ Mpc}$, and $\sigma_{\text{wall}} = 35\text{ Mpc}$.
            * **Path-Integrated Effective Expansion $H_{\text{eff}}(r)$:**
              $$cz(r) = \int_0^r H_{\text{local}}(s)\,ds \implies H_{\text{eff}}(r) = \frac{cz(r)}{r}$$
            * **Vacuum Floor Relaxation:** At $r > 300\text{ Mpc}$, path integration averages over background void volume, locking integrated expansion strictly onto $H_{\text{global}} = 67.42\text{ km/s/Mpc}$.

            ---
            ### ⚙️ Interactive Controls & Parameter Reference
            * **Central Density Deficit ($\delta_0$):** Controls the depth of the local KBC supervoid count deficit (default $-0.22$ or $-22\%$).
            * **Local Core Peak ($H_{\text{peak}}$):** Sets the foreground Local Sheet peak expansion rate (default $72.80\text{ km/s/Mpc}$).
            * **Global Vacuum Floor ($H_{\text{global}}$):** Fundamental hardware floor set by Planck 2018 ($67.42\text{ km/s/Mpc}$).

            ---
            ### 📡 Model Response & Horizon Relaxation Guidelines
            * **Foreground ($r \le 15\text{ Mpc}$):** Dominated by the Local Sheet overdensity; $H_{\text{eff}} \approx 70.0\text{--}72.8\text{ km/s/Mpc}$.
            * **Intermediate Volume ($50\text{--}150\text{ Mpc}$):** Path integration smooths expansion down to $67.6\text{--}68.1\text{ km/s/Mpc}$, matching Pantheon+ SNe Ia.
            * **Horizon Limit ($r \ge 300\text{ Mpc}$):** Complete relaxation to $67.42\text{ km/s/Mpc}$, agreeing with DESI / eBOSS BAO and Planck CMB baselines.
            """)

        col_kbc1, col_kbc2 = st.columns(2)
        delta_kbc_input = col_kbc1.slider("Central Density Deficit delta_0", -0.35, -0.10, -0.22, 0.01)
        h_peak_input = col_kbc2.slider("Local Filament Core Peak H_peak [km/s/Mpc]", 70.0, 76.0, 72.80, 0.2)
        
        run_kbc = st.button("🚀 Execute KBC Void Relaxation Analysis", type="primary", use_container_width=True)
        m_key = "mod_res_kbc_void"

        if run_kbc:
            with st.spinner("Calculating KBC Supervoid radial transition..."):
                out_kbc = run_kbc_void_analysis(r_max_mpc=500.0, delta_0_kbc=delta_kbc_input, h_local_peak=h_peak_input, h_global=h_global_input)
                fig_kbc, df_kbc, kbc_metrics = (out_kbc[0], out_kbc[1], out_kbc[2]) if len(out_kbc) == 3 else (out_kbc[0], None, out_kbc[1])
                st.session_state[m_key] = (fig_kbc, df_kbc, kbc_metrics)

        if m_key in st.session_state:
            fig_kbc, df_kbc, kbc_metrics = st.session_state[m_key]
            display_module_outputs(fig_kbc, df_kbc, kbc_metrics, "kbc_supervoid_radial_profile.csv")
        else:
            st.info("⚙️ **Ready to Compute:** Set density deficit and peak expansion above, then click **🚀 Execute KBC Void Relaxation Analysis**.")

# ------------------------------------------------------------------------------
# DOMAIN 2: COSMIC WEB & VELOCITY DYNAMICS
# ------------------------------------------------------------------------------
with domain_2:
    sub_2 = st.radio(
        "Select Module:",
        ["🌊 Bulk Flow Field", "🧱 Void Boundary Dynamics", "📐 Redshift Space Distortions (Kaiser & FoG)", "🌀 Galactic Rotation Curves & Kinematic Underflow", "🧊 Interactive 3D Cosmic Web Density Viewer"],
        horizontal=True
    )
    st.divider()

    # 6. BULK FLOW FIELD
    if "Bulk Flow" in sub_2:
        st.header("Local Sheet & Macro Bulk Flow Field")
        st.info(
            r"🌊 **Physical Principle:** In Substrate Logistics, peculiar velocity is driven by spatial gradients in the local unspooling rate: "
            r"$\vec{v}_{\text{inflow}} = \alpha \cdot \nabla H_{\text{local}}(\vec{r})$. "
            r"Galaxies drift toward regions of higher compiled density where unspooling rates are elevated."
        )

        with st.expander("📖 **Interactive Parameter & Demonstration Guide**", expanded=False):
            st.markdown(r"""
            ### What This Module Demonstrates
            * **Peculiar Velocity Generation:** Replaces abstract Newtonian dark matter pulling forces with a direct physical response to spatial metric gradients ($\nabla H_{\text{local}}$).
            * **Scale-Dependent Flows:** Demonstrates how local sharp filament infalls ($10\text{--}20\text{ Mpc}$) transition into macro bulk flow channels ($100\text{--}200\text{ Mpc}$) directed toward supercluster attractors (e.g., Laniakea).
            * **Void Outflows & Filament Inflows:** Proves that vector quivers diverge away from flat void centers ($\nabla H \to 0$) and converge onto dense wall corridors.

            ---
            ### 🔬 Key Physical Mechanics
            * **Expansion Gradient Drift:** Peculiar motion is modeled as coordinate conveyor drift directed toward regions of elevated compiled density:
              $$\vec{v}_{\text{inflow}}(\vec{r}) = +\alpha \nabla H_{\text{local}}(\vec{r})$$
            * **Void Floor Motionless Baseline:** Inside deep voids ($\rho_N \to 0$), expansion locks onto $H_{\text{global}} = 67.42\text{ km/s/Mpc}$. Because $\nabla H = 0$, peculiar velocity naturally drops to $0.00\text{ km/s}$.

            ---
            ### ⚙️ Interactive Controls & Parameter Reference
            * **Gradient Coupling ($\alpha$):** Linear coupling constant ($\text{Mpc}^2$). Controls how strongly an expansion gradient $\nabla H$ translates into physical peculiar velocity ($\text{km/s}$).
            * **Spatial Extent [Mpc]:** Sets the radius of the 2D evaluation plane. Zooming out to $100\text{--}200\text{ Mpc}$ averages over large cosmic voids, mapping macro-scale flow channels.
            * **Grid Resolution ($N \times N$):** Spatial sampling density. Higher grid resolutions capture fine-grained density spikes around compact cluster cores.
            * **Field Smoothing Scale ($\sigma$ in Sidebar):** Gaussian bandwidth kernel used to construct continuous density fields $\rho_N(\vec{r})$ from catalog point galaxies.

            ---
            ### 📡 Survey Catalog Response & Performance Guidelines
            * **UNG 2013 (Local Volume ≤ 35 Mpc | $\sigma = 1.8\text{ Mpc}$):** Resolves tight local sheet flows ($3.36\text{ km/s}$ mean drift) and sharp cluster infall peaks (up to $80.38\text{ km/s}$).
            * **Cosmicflows-4 / 2M++ (Extended ≤ 150–200 Mpc | $\sigma = 15.0\text{ Mpc}$):** Averages over macro cosmic web volumes ($V \sim 10^7\text{ Mpc}^3$), mapping large-scale flow channels with mean drift speeds of $\sim 3.85\text{ km/s}$ and peak filament infalls of $10.88\text{ km/s}$.
            """)

        col_bf1, col_bf2, col_bf3 = st.columns(3)
        grid_n = col_bf1.slider("Grid Resolution (N x N)", 20, 50, 30, 5, help="Sampling density. Higher values capture sharper filament gradients.")
        bounds_val = col_bf2.slider("Spatial Extent [Mpc]", 5.0, float(d_max_input), float(min(25.0, float(d_max_input))), 5.0, help="Physical size of the 2D slice plane box.")
        alpha_scale = col_bf3.slider("Gradient Coupling alpha", 50.0, 200.0, 120.0, 10.0, help="Coupling strength converting expansion gradients grad H into velocity (km/s).")

        run_bf = st.button("🚀 Execute Bulk Flow Vector Field", type="primary", use_container_width=True)
        m_key = "mod_res_bulk_flow"

        if run_bf:
            with st.spinner("Computing 2D spatial gradient velocity field..."):
                out_bf = run_bulk_flow_analysis(
                    tree, gal_positions, nodes_per_galaxy, X_gal, Y_gal, Z_gal, 
                    grid_size=grid_n, spatial_bounds=bounds_val, 
                    sigma_mpc=sigma_smooth_input, h_global=h_global_input, alpha_scale=alpha_scale
                )
                fig_bf, df_bf, bf_metrics = (out_bf[0], out_bf[1], out_bf[2]) if len(out_bf) == 3 else (out_bf[0], None, out_bf[1])
                st.session_state[m_key] = (fig_bf, df_bf, bf_metrics)

        if m_key in st.session_state:
            fig_bf, df_bf, bf_metrics = st.session_state[m_key]
            display_module_outputs(fig_bf, df_bf, bf_metrics, "bulk_flow_velocity_grid.csv")
        else:
            st.info("⚙️ **Ready to Compute:** Select grid resolution and coupling above, then click **🚀 Execute Bulk Flow Vector Field**.")

    # 7. VOID BOUNDARY DYNAMICS
    elif "Void Boundary" in sub_2:
        st.header("Cosmic Void Wall Boundary Dynamics & Velocity Shear")
        st.info(
            "🧱 **Physical Principle:** Deep inside cosmic voids ($\rho_N \to 0$), local expansion locks onto the global vacuum floor ($H_{\text{global}} = 67.42\\text{ km/s/Mpc}$). "
            "Crossing into overdense wall corridors increases compiled node density, driving radial expansion boosts and boundary velocity shear spikes."
        )

        with st.expander("📖 **Interactive Parameter & Demonstration Guide**", expanded=False):
            st.markdown(r"""
            ### What This Module Demonstrates
            * **Void Core Vacuum Floor Locking ($r \to 0$):** Deep inside empty cosmic voids where matter density vanishes ($\rho_{\text{macro}} \to 0$), expansion locks onto the fundamental hardware floor ($H_{\text{global}} = 67.42\text{ km/s/Mpc}$).
            * **Wall Unspooling Field Boost ($r \sim 10\text{--}16\text{ Mpc}$):** As light rays cross into surrounding galaxy wall structures, compiled matter density smoothly boosts local expansion up to $71\text{--}76\text{ km/s/Mpc}$.
            * **Boundary Velocity Shear Spike ($dv_{\text{rec}}/dr$):** Differential expansion across wall boundaries generates a velocity shear spike ($\sim 72\text{--}89\text{ km/s/Mpc}$) that dynamically evacuates voids without requiring dark energy.

            ---
            ### 🔬 Key Physical Mechanics
            * **Dual-Scale Density Kernel:** Uses internal fixed kernels ($\sigma_{\text{macro}} = 1.5\text{ Mpc}$ and $\sigma_{\text{filament}} = 0.5\text{ Mpc}$) to decouple macro intergalactic background volume from fine filament core density.
            * **Density-Driven Unspooling:**
              $$H_{\text{local}}(\vec{r}) = H_{\text{global}} + \frac{1}{\text{UNIT\_CONV}} \sqrt{\frac{8\pi \mathcal{K}_\Omega \rho_{\text{macro}}(\vec{r})}{3}}$$
            * **Boundary Velocity Shear:**
              $$\frac{dv_{\text{rec}}}{dr} = H_{\text{local}}(r) + r \frac{dH_{\text{local}}}{dr}$$

            ---
            ### ℹ️ Kernel Architecture Note
            * **Sidebar Smoothing Scale ($\sigma$):** The global sidebar $\sigma$ slider does **not** affect this module. Void boundary dynamics utilizes hardcoded dual-scale kernels ($\sigma_{\text{macro}} = 1.5\text{ Mpc}$, $\sigma_{\text{filament}} = 0.5\text{ Mpc}$) to maintain strict spatial isolation between filament cores and void boundaries across all catalogs.

            ---
            ### 📡 Survey Catalog Response
            * **UNG 2013 (Local Volume):** Captures the immediate Local Void boundary out to $4.5\text{ Mpc}$ ($H_{\text{peak}} \approx 68.35\text{ km/s/Mpc}$).
            * **Cosmicflows-4 / 2M++:** Traces deep wall corridors (e.g., Virgo Wall at $15.7\text{ Mpc}$), resolving velocity shear spikes up to $88.74\text{ km/s/Mpc}$.
            """)

        run_vd = st.button("🚀 Execute Void Boundary Shear Analysis", type="primary", use_container_width=True)
        m_key = "mod_res_void_boundary"

        if run_vd:
            with st.spinner("Evaluating void boundary wall velocity shear..."):
                out_vd = run_void_boundary_analysis(
                    tree, gal_positions, nodes_per_galaxy, h_global=h_global_input
                )
                fig_vd, df_vd, vd_metrics = (out_vd[0], out_vd[1], out_vd[2]) if len(out_vd) == 3 else (out_vd[0], None, out_vd[1])
                st.session_state[m_key] = (fig_vd, df_vd, vd_metrics)

        if m_key in st.session_state:
            fig_vd, df_vd, vd_metrics = st.session_state[m_key]
            display_module_outputs(fig_vd, df_vd, vd_metrics, "void_boundary_wall_profile.csv")
        else:
            st.info("⚙️ **Ready to Compute:** Review active survey parameters and click **🚀 Execute Void Boundary Shear Analysis**.")

    # 8. REDSHIFT SPACE DISTORTIONS
    elif "Redshift Space Distortions" in sub_2:
        st.header("2D Anisotropic Redshift Space Distortions (Kaiser & FoG)")
        st.info(
            r"📐 **Physical Principle:** In Substrate Logistics, redshift space distortions do not rely on separate empirical parameter fits ($f\sigma_8$ vs $\sigma_v$). "
            r"Large-scale coherent infall ($v_{\text{infall}} = \alpha \cdot \nabla H_{\text{local}}$) and cluster-core virial dispersion ($v_{\text{virial}} \propto \sqrt{\max(0, \rho_N/\bar{\rho}_{\text{baseline}} - 1)}$) "
            r"are evaluated simultaneously zero-parameter directly from the active catalog's 3D compiled node density field $\rho_N(\vec{r})$."
        )

        with st.expander("📖 **Interactive Parameter & Demonstration Guide**", expanded=False):
            st.markdown(r"""
            ### What This Module Demonstrates
            * **Unified Zero-Parameter Mapping:** Eliminates empirical $\Lambda\text{CDM}$ curve fitting ($f\sigma_8$, $\sigma_v$) by deriving both Kaiser squashing and Fingers-of-God elongation simultaneously from the active catalog's 3D compiled density field $\rho_N(\vec{r})$.
            * **Kaiser Coherent Squashing:** On large scales ($\sigma > 5.0\text{ Mpc}$), spatial unspooling gradients $\nabla H_{\text{local}}$ drive inward galaxy drifts ($v_{\text{infall}}$), squashing outer 2D correlation contours $\xi(\sigma, \pi)$ horizontally along the line-of-sight axis ($\pi$).
            * **Fingers-of-God Virial Elongation:** In overdense cluster cores ($\rho_N > \bar{\rho}_{\text{baseline}}$), high orbital dispersion velocities ($v_{\text{virial}}$) stretch inner correlation contours vertically into elongated "fingers".

            ---
            ### 🔬 Key Physical Mechanics
            * **Dynamic Catalog Baseline ($\bar{\rho}_{\text{baseline}}$):** Evaluates volumetric mean compiled node density dynamically across the survey sphere:
              $$\bar{\rho}_{\text{baseline}} = \frac{\sum N_i}{\frac{4}{3}\pi R_{\text{max}}^3} \quad [\text{nodes/m}^3]$$
            * **Redshift-Space Coordinate Transformation:**
              $$\vec{s} = \sigma \hat{e}_\perp + \left(\pi_{\text{real}} + \frac{v_{\text{infall},\parallel}(\vec{r}) + v_{\text{virial},\parallel}(\vec{r})}{H_{\text{global}}}\right)\hat{e}_\parallel$$

            ---
            ### ⚙️ Interactive Controls & Parameter Reference
            * **Infall Gradient Coupling ($\alpha$):** Linear coupling constant ($\text{Mpc}^2$) converting unspooling gradients $\nabla H_{\text{local}}$ into coherent line-of-sight infall velocity ($v_{\text{infall}}$ in $\text{km/s}$).
            * **Field Smoothing Scale ($\sigma$ in Sidebar):** Gaussian bandwidth kernel used to construct continuous 3D density field $\rho_N(\vec{r})$. Lower settings expose sharp cluster core FoG elongation, whereas higher settings emphasize macro Kaiser squashing.

            ---
            ### 📡 Survey Catalog Response & Performance Guidelines
            * **UNG 2013 (Local Volume):** Traces local sheet dynamics ($d \le 35\text{ Mpc}$); resolves local group virial dispersion.
            * **Cosmicflows-4 / 2M++:** Evaluates extended survey volumes ($150\text{--}200\text{ Mpc}$). At $\sigma = 5.0\text{ Mpc}$, samples ~4,200 overdense regions with mean 1-$\sigma$ virial dispersion $v_{\text{virial}} \approx 546.4\text{ km/s}$ and peak infall velocity $v_{\text{infall}} \approx 93.0\text{ km/s}$.
            """)

        col_rsd1, _ = st.columns([1, 1])
        alpha_rsd = col_rsd1.slider("Infall Gradient Coupling alpha", 50.0, 200.0, 120.0, 10.0, help="Linear coupling constant converting expansion gradients grad H into coherent infall velocity (km/s).")

        run_rsd = st.button("🚀 Execute 2D RSD Correlation Map", type="primary", use_container_width=True)
        m_key = "mod_res_rsd"

        if run_rsd:
            with st.spinner("Generating 2D Kaiser & Fingers-of-God correlation contours..."):
                out_rsd = run_redshift_space_distortions(
                    tree=tree, gal_positions=gal_positions, nodes_per_galaxy=nodes_per_galaxy,
                    alpha_inflow=alpha_rsd, sigma_mpc=sigma_smooth_input, h_global=h_global_input
                )
                fig_rsd, df_rsd, rsd_metrics = (out_rsd[0], out_rsd[1], out_rsd[2])
                st.session_state[m_key] = (fig_rsd, df_rsd, rsd_metrics)

        if m_key in st.session_state:
            fig_rsd, df_rsd, rsd_metrics = st.session_state[m_key]
            display_module_outputs(fig_rsd, df_rsd, rsd_metrics, "redshift_space_distortion_2d_correlation.csv")
        else:
            st.info("⚙️ **Ready to Compute:** Select gradient coupling above and click **🚀 Execute 2D RSD Correlation Map**.")

    # 12. GALACTIC ROTATION CURVES
    elif "Galactic Rotation" in sub_2:
        st.header("Galactic Rotation Curves & Kinematic Underflow Floor")
        st.info(
            r"🌀 **Physical Principle:** In Substrate Logistics, flat galactic rotation curves require zero Cold Dark Matter halos. "
            r"When localized baryonic acceleration drops below the dynamic Universal Baseline Acceleration ($g_{\text{baryon}} < a_\Omega = \frac{c H_{\text{local}}}{2\pi}$), "
            r"the grid encounters a Kinematic Underflow Error, clamping outer orbital speeds to "
            r"$v_{\text{flat}} = \sqrt[4]{N_{\text{total}} K_\Omega a_\Omega}$ (equivalent to $\sqrt[4]{G M_B a_\Omega}$ in classical limits)."
        )
        
        with st.expander("📖 **Interactive Parameter & Demonstration Guide**", expanded=False):
            st.markdown(r"""
            ### What This Module Demonstrates
            * **Zero Dark Matter Halos:** Replaces dark matter halo fitting with a single hardware limit—the Kinematic Underflow Floor ($a_\Omega$).
            * **Flat Outer Rotation Curves:** Demonstrates how outer orbital speeds lock onto a constant plateau $v_{\text{flat}} = \sqrt[4]{N_{\text{total}} \mathcal{K}_\Omega a_\Omega}$ as baryonic acceleration decays.
            * **Baryonic Tully-Fisher Relation:** Proves that $v_{\text{flat}} \propto M_B^{1/4}$ is a direct algebraic consequence of metric underflow clamping rather than an empirical mystery.
            * **Live 3D Density Coupling:** Shows how local cosmic environment ($H_{\text{local}}$) dynamically alters the underflow threshold $a_\Omega(\mathbf{x}) = \frac{c H_{\text{local}}(\mathbf{x})}{2\pi}$.

            ---
            ### 🔬 Key Physical Mechanics
            * **Kinematic Underflow Coupling:** When baryonic acceleration $g_{\text{baryon}}$ drops below $a_\Omega$, total effective acceleration scales geometrically:
              $$g_{\text{eff}}(r) = \sqrt{g_{\text{baryon}}(r)^2 + g_{\text{baryon}}(r) \cdot a_\Omega}$$
            * **Proton Node Gravitational Constant ($\mathcal{K}_\Omega$):** Acceleration is evaluated directly from proton counts $N_{\text{total}} = M_B / m_p$ without empirical $G$ fitting:
              $$g_{\text{baryon}}(r) = \frac{N_{\text{total}} \mathcal{K}_\Omega}{r^2} \quad \left(\equiv \frac{G M_B}{r^2}\right)$$
            * **Transition Radius ($R^*$):** The exact boundary where effective acceleration crosses the underflow floor ($g_{\text{eff}} = a_\Omega$), occurring at $g_{\text{baryon}}(R^*) \approx 0.618 \cdot a_\Omega$.

            ---
            ### ⚙️ Interactive Controls & Parameter Reference
            * **Rotation Curve Mode:**
              * *Benchmark Galaxy Presets:* Load standard targets (NGC 3198, Milky Way, DDO 154) to inspect classic rotation curves.
              * *Live Catalog Query:* Select any real galaxy from the active 3D survey catalog to extract its exact spatial position and evaluate $a_\Omega$ from local 3D node density.
            * **Baryonic Mass ($M_B$):** Total visible baryonic mass (stars + gas) in solar masses ($M_\odot$).
            * **Kernel Scale ($\sigma$ in Sidebar):** Gaussian smoothing window used for 3D density evaluation during catalog queries.

            ---
            ### 📡 Survey Catalog Response & Performance Guidelines
            * **Local Volume (UNG 2013):** Samples nearby galaxies within $35\text{ Mpc}$, isolating local group density fluctuations.
            * **Extended Catalogs (Cosmicflows-4 / 2M++):** Evaluates galaxies in deep cosmic web environments out to $200\text{ Mpc}$, showing how cluster environments boost $a_\Omega$.
            """)

        mode_rot = st.radio(
            "Select Rotation Curve Mode:",
            ["🎯 Benchmark Galaxy Presets (SPARC / MW / Dwarf)", "📡 Live Catalog Query (3D Density-Coupled a_Ω)"],
            horizontal=True
        )

        if "Benchmark" in mode_rot:
            col_r1, col_r2 = st.columns(2)
            target_preset_rot = col_r1.selectbox(
                "Select Galaxy Preset:",
                list(ROTATION_PRESETS.keys())
            )

            default_m_solar = ROTATION_PRESETS[target_preset_rot]["M_baryon_solar"]
            
            if default_m_solar >= 1e11:
                init_val, unit_lbl, scale_fac, step_v = default_m_solar / 1e12, "Baryonic Mass [10¹² M☉]:", 1e12, 0.05
            elif default_m_solar >= 1e8:
                init_val, unit_lbl, scale_fac, step_v = default_m_solar / 1e9, "Baryonic Mass [10⁹ M☉]:", 1e9, 1.0
            else:
                init_val, unit_lbl, scale_fac, step_v = default_m_solar / 1e6, "Baryonic Mass [10⁶ M☉]:", 1e6, 10.0

            custom_m_rot_scaled = col_r2.number_input(
                unit_lbl,
                value=float(init_val),
                step=float(step_v),
                help="Baryonic mass baseline. No dark matter halos!"
            )
            custom_m_rot_solar = custom_m_rot_scaled * scale_fac
            use_cat_rot = False
            target_pos_rot = None
            sigma_k_rot = 1.8

        else: # Live Catalog Query Mode
            st.markdown("### 📡 Live 3D Catalog Density Sampling")
            
            gal_options_rot = []
            for i in range(len(gal_names)):
                d_mpc = float(np.linalg.norm(gal_positions[i]))
                m_solar = float((nodes_per_galaxy[i] * M_PROTON) / M_SOLAR_KG)
                m_str = format_baryonic_mass(m_solar)
                gal_options_rot.append(f"#{i}: {gal_names[i]} (d={d_mpc:.1f} Mpc | Mass={m_str})")

            col_rc1, col_rc2, col_rc3 = st.columns([2.2, 1.0, 1.2])
            selected_opt_rot = col_rc1.selectbox(
                "Select Target Galaxy from Catalog:",
                options=gal_options_rot,
                index=min(5, len(gal_options_rot)-1)
            )

            gal_idx_rot = int(selected_opt_rot.split(":")[0].replace("#", ""))
            target_pos_rot = gal_positions[gal_idx_rot]
            target_name_rot = gal_names[gal_idx_rot]
            cat_m_solar_rot = float((nodes_per_galaxy[gal_idx_rot] * M_PROTON) / M_SOLAR_KG)

            col_rc1.caption(f"📍 **Position:** X={target_pos_rot[0]:.2f}, Y={target_pos_rot[1]:.2f}, Z={target_pos_rot[2]:.2f} Mpc | **Name:** {target_name_rot}")

            sigma_k_rot = col_rc2.slider(
                "Kernel Scale σ [Mpc]:",
                min_value=0.5, max_value=5.0, value=1.8, step=0.1
            )

            if cat_m_solar_rot >= 1e11:
                in_v, u_lbl, s_fac = cat_m_solar_rot / 1e12, "Baryonic Mass [10¹² M☉]:", 1e12
            elif cat_m_solar_rot >= 1e8:
                in_v, u_lbl, s_fac = cat_m_solar_rot / 1e9, "Baryonic Mass [10⁹ M☉]:", 1e9
            else:
                in_v, u_lbl, s_fac = cat_m_solar_rot / 1e6, "Baryonic Mass [10⁶ M☉]:", 1e6

            custom_m_rot_scaled = col_rc3.number_input(u_lbl, value=float(in_v), step=1.0)
            custom_m_rot_solar = custom_m_rot_scaled * s_fac
            use_cat_rot = True
            target_preset_rot = f"Catalog: {target_name_rot}"

        run_rot = st.button("🚀 Execute Galactic Rotation Analysis", type="primary", use_container_width=True)
        m_key = "mod_res_rotation"

        if run_rot:
            with st.spinner("Calculating rotation curve & Kinematic Underflow floor..."):
                fig_rot, df_rot, rot_metrics = run_rotation_curve_analysis(
                    galaxy_target_name=target_preset_rot,
                    custom_M_baryon=custom_m_rot_solar,
                    h_global=h_global_input,
                    use_catalog_query=use_cat_rot,
                    target_pos_mpc=target_pos_rot,
                    sigma_mpc=sigma_k_rot,
                    tree=tree,
                    gal_positions=gal_positions,
                    nodes_per_galaxy=nodes_per_galaxy
                )
                st.session_state[m_key] = (fig_rot, df_rot, rot_metrics)

        if m_key in st.session_state:
            fig_rot, df_rot, rot_metrics = st.session_state[m_key]
            display_module_outputs(fig_rot, df_rot, rot_metrics, "substrate_galactic_rotation_curve.csv")
        else:
            st.info("⚙️ **Ready to Compute:** Select galaxy above and click **🚀 Execute Galactic Rotation Analysis**.")

    # 13. INTERACTIVE 3D COSMIC WEB DENSITY VIEWER
    elif "Interactive 3D Cosmic Web" in sub_2:
        st.header("Interactive 3D Cosmic Web Density Viewer & Filament Isosurfaces")
        st.info(
            r"🧊 **Physical Principle:** Interactive 3D Cartesian mapping $(X, Y, Z)$ of active survey catalog galaxies. "
            r"Directly visualizes zero-parameter compiled baryonic mass distribution $M_B$, local metric unspooling rate field $H_{\text{local}}(\vec{r})$, "
            r"and volumetric intergalactic filament web isosurfaces across cosmic volumes up to $150\text{--}200\text{ Mpc}$."
        )

        with st.expander("📖 **Interactive Parameter & Demonstration Guide**", expanded=False):
            st.markdown(r"""
            ### What This Module Demonstrates
            * **3D Cosmic Web Rendering:** Interactively maps galaxy positions in ICRS Cartesian coordinates $(X, Y, Z)$ [Mpc] alongside 3D volumetric filament isosurface meshes generated by local density gradients.
            * **Zero-Parameter Mass & Metric Coupling:** Color-codes each galaxy node by its localized unspooling expansion rate $H_{\text{local}}(\vec{r})$, derived directly from compiled baryonic matter ($M_B = N \cdot m_p$) without dark matter halos.
            * **Volumetric Isosurface Web:** Visualizes 3D overdense filament corridors and cluster nodes where compiled matter accelerates local metric expansion above the global vacuum floor ($H_{\text{global}} = 67.42\text{ km/s/Mpc}$).

            ---
            ### 🔬 Key Physical Mechanics
            * **Local Metric Unspooling Field:** Evaluates local expansion $H_{\text{local}}(\vec{r})$ directly from the 3D compiled node density field $\rho_N(\vec{r})$:
              $$H_{\text{local}}(\vec{r}) = H_{\text{global}} + \sqrt{\frac{8\pi \mathcal{K}_\Omega \rho_N(\vec{r})}{3}}$$
            * **Baryonic Node Mass Scaling:** Galaxy marker sizes scale logarithmically with compiled baryonic mass $M_B = N \cdot m_p$, replacing invisible Dark Matter halos with direct photometric luminosity compilation.

            ---
            ### ⚙️ Interactive Controls & Parameter Reference
            * **Max Display Galaxies (`max_galaxies`):** Sets the point display cap (e.g., 500 to 5,000 nodes) to maintain smooth Plotly WebGL rendering frame rates across large survey volumes.
            * **3D Isosurface Grid Resolution (`grid_res`):** Sampling density ($N \times N \times N$) for evaluating the 3D scalar volumetric unspooling field $H_{\text{local}}(\vec{r})$.
            * **Field Smoothing Scale ($\sigma$ in Sidebar):** Gaussian kernel bandwidth ($\text{Mpc}$) used for 3D density evaluation.
              * *Astrophysical Macro Scale ($\sigma = 13.30\text{ Mpc}$ | $3\sigma = 39.9\text{ Mpc}$ search radius):* Samples the standard enclosed volume ($\sim 2.66 \times 10^5\text{ Mpc}^3$) where astrophysics evaluates macro Hubble flow expansion.
              * *Local Fine-Grained ($\sigma = 1.8\text{--}3.0\text{ Mpc}$):* Resolves sharp localized filament channels and individual cluster density spikes.

            ---
            ### 📡 Survey Catalog Response & Performance Guidelines
            * **Cosmicflows-4 ($150\text{ Mpc}$ Benchmark):** Evaluated at macro scale ($\sigma = 13.30\text{ Mpc}$), $H_{\text{local}}$ spans **$68.04\text{ to }73.80\text{ km/s/Mpc}$** (mean $70.22\text{ km/s/Mpc}$), smoothly bridging local filament peaks down to background void floors.
            * **UNG 2013 & 2M++ Infrared:** UNG 2013 isolates the Local Volume ($d \le 35\text{ Mpc}$), whereas 2M++ extends full-sky coverage out to $200\text{ Mpc}$.
            """)

        col_v1, col_v2, col_v3 = st.columns(3)
        max_gal_input = col_v1.slider("Max Display Galaxies", 500, 5000, 2000, 250, help="Node display cap for smooth rendering.")
        grid_res_input = col_v2.slider("3D Isosurface Grid Resolution", 12, 30, 18, 2, help="Volumetric sampling grid density.")
        show_iso = col_v3.checkbox("Render 3D Filament Web Mesh", value=True)
        iso_opacity = col_v3.slider("Isosurface Mesh Opacity", 0.05, 0.60, 0.20, 0.05) if show_iso else 0.25

        run_3d = st.button("🚀 Render Interactive 3D Density Field", type="primary", use_container_width=True)
        m_key = "mod_res_3d_density_viewer"

        if run_3d:
            with st.spinner("Building interactive 3D volumetric web figure..."):
                fig_3d, df_3d, metrics_3d = run_3d_density_viewer_analysis(
                    tree, gal_positions, nodes_per_galaxy, gal_names,
                    sigma_mpc=sigma_smooth_input,
                    h_global=h_global_input,
                    grid_res=grid_res_input,
                    max_galaxies=max_gal_input,
                    show_isosurface=show_iso,
                    isosurface_opacity=iso_opacity
                )
                st.session_state[m_key] = (fig_3d, df_3d, metrics_3d)

        if m_key in st.session_state:
            fig_3d, df_3d, metrics_3d = st.session_state[m_key]
            
            # Display Top Metric KPI Cards
            if metrics_3d:
                cols = st.columns(len(metrics_3d))
                for col, (k, v) in zip(cols, metrics_3d.items()):
                    col.metric(k, v)

            # Render Plotly Interactive Chart
            st.plotly_chart(fig_3d, use_container_width=True)

            # Display Download Button & Summary Dataframe
            if df_3d is not None and not df_3d.empty:
                st.subheader("🎯 Active 3D Rendered Galaxies")
                st.download_button(
                    label="📥 Export 3D Galaxy Dataset (CSV)",
                    data=df_3d.to_csv(index=False).encode('utf-8'),
                    file_name="3d_density_web_galaxies.csv",
                    mime="text/csv",
                    type="primary"
                )
                st.dataframe(df_3d, use_container_width=True)
        else:
            st.info("⚙️ **Ready to Render:** Configure display controls above and click **🚀 Render Interactive 3D Density Field**.")

# ------------------------------------------------------------------------------
# DOMAIN 3: ASTROPHYSICAL & HORIZON TESTS
# ------------------------------------------------------------------------------
with domain_3:
    sub_3 = st.radio(
        "Select Module:",
        ["📏 BAO Ruler Warping", "🛰️ CMB Dipole Anomaly", "🔍 Zero-Dark-Matter Gravitational Lensing"],
        horizontal=True
    )
    st.divider()

    # 9. BAO RULER WARPING
    if "BAO Ruler" in sub_3:
        st.header("Baryon Acoustic Oscillation Ruler Warping & AP Test")
        st.info(
            r"📏 **First-Principles Derivation:** The BAO sound horizon ($r_s = 149.21\text{ Mpc}$) is derived strictly "
            r"from structural invariants without thermodynamic tuning: the maximum causal horizon $R_{\text{max}} = c \cdot \mathcal{T}_\Omega = 3.9017 \times 10^{27}\text{ m}$ "
            r"is subjected to 3D topological drag penalty $R_{\text{restricted}} = R_{\text{max}} \cdot (3 \cdot \mathcal{C}_\delta) = 5.5248 \times 10^{24}\text{ m}$ "
            r"and divided by the dodecahedral compilation scalar ($1.2$) to yield $R_{\text{BAO}} = 4.6040 \times 10^{24}\text{ m} \equiv \mathbf{149.21\text{ Mpc}}$. "
            r"As photon rays traverse overdense filaments, higher compiled density boosts local expansion ($H_{\text{eff}} > H_{\text{global}}$), "
            r"deforming the acoustic sphere into an anisotropic ellipsoid: **radial squashing ($r_\parallel < r_s$)** and "
            r"**transverse dilation ($r_\perp > r_s$)**, quantified by $F_{\text{AP}} = \alpha_\perp / \alpha_\parallel = C^{3/2}$."
        )

        with st.expander("📖 **Interactive Parameter & Demonstration Guide**", expanded=False):
            st.markdown(r"""
            ### What This Module Demonstrates
            * **Anisotropic Geometric Warping:** Demonstrates how an ideal $149.21\text{ Mpc}$ spherical acoustic shell deforms into an elongated ellipsoid along overdense supercluster corridors ($F_{\text{AP}} \approx 1.050\text{--}1.127$) while remaining nearly isotropic across deep cosmic voids ($F_{\text{AP}} \approx 1.0055$, $+0.55\%$).
            * **Validation of Local BAO Structures (Ho'oleilana):** Provides a physical mechanism for empirical discoveries like *Ho'oleilana* in Cosmicflows-4—a giant local BAO structure dilated to $r_\perp \approx 152.5\text{--}155.3\text{ Mpc}$—without requiring dark energy fine-tuning or anomalous Hubble constants.
            * **Horizon Volume Relaxation:** Demonstrates that path integration out to the true sound horizon ($149.21\text{ Mpc}$) averages over local void and filament fluctuations, bringing void corridor predictions into strict agreement with DESI / eBOSS isotropy bounds ($\le 0.6\%$).

            ---
            ### 🔬 Key Physical Mechanics
            * **Bare-Metal Acoustic Horizon ($r_s = 149.21\text{ Mpc}$):**
              $$r_s = \frac{c \cdot \mathcal{T}_\Omega \cdot (3 \mathcal{C}_\delta)}{1.2} = \frac{(2.9979 \times 10^8) \cdot (1.3015 \times 10^{19}) \cdot (1.416 \times 10^{-3})}{1.2} \equiv 149.21\text{ Mpc}$$
            * **Alcock-Paczynski Distortion Parameter ($F_{\text{AP}}$):**
              $$C = \frac{H_{\text{eff}}}{H_{\text{global}}}, \quad \alpha_\parallel = \frac{1}{C}, \quad \alpha_\perp = \sqrt{C} \implies F_{\text{AP}} = \frac{\alpha_\perp}{\alpha_\parallel} = C^{3/2}$$

            ---
            ### ⚙️ Interactive Controls & Parameter Reference
            * **Ray Integration Depth [Mpc]:** Distance integrated along ray paths.
              * *Full Horizon ($149.21\text{ Mpc}$ - Recommended):* Integrates across the true physical acoustic sound horizon, correctly balancing foreground cluster nodes with background void volume.
              * *Local Foreground ($20.0\text{ Mpc}$):* Samples only immediate local galaxy neighborhoods, exaggerating local density spikes.
            * **AP Scan Declination [deg]:** Declination angle for the $360^\circ$ Right Ascension sky scan plane, mapping anisotropic dilation variations across 2D slices of the cosmic web.
            * **Field Smoothing Scale ($\sigma$ in Sidebar):** Gaussian kernel bandwidth used to evaluate continuous 3D node density fields $\rho_N(\vec{r})$.
              * *Coarse ($\sigma = 5.0\text{ Mpc}$):* Smooths discrete cluster core spikes (e.g., Coma $F_{\text{AP}} = 1.1272$, $+12.72\%$; Virgo $F_{\text{AP}} = 1.0677$, $+6.77\%$).
              * *Fine ($\sigma = 1.8\text{ Mpc}$):* Resolves sharp virial cluster boundaries, increasing peak node anisotropy while maintaining deep void floor convergence.

            ---
            ### 📡 Survey Catalog Response & Performance Guidelines
            * **Cosmicflows-4 ($\le 150\text{ Mpc}$):** Maps major cosmic web features (Coma, Virgo, Centaurus, Perseus-Pisces), predicting transverse dilations of $r_\perp = 151.6\text{--}155.3\text{ Mpc}$ along overdense corridors.
            * **2M++ Infrared ($\le 200\text{ Mpc}$):** Extends path integration through full supercluster complexes, smoothing void floor dilation down to $F_{\text{AP}} < 1.002$.
            """)

        col_bao1, col_bao2 = st.columns(2)
        depth_bao_input = col_bao1.slider("Ray Integration Depth [Mpc]", 20.0, float(d_max_input), 149.21, 10.0, help="Setting depth to 149.21 Mpc integrates across the true physical BAO acoustic scale.")
        dec_scan_input = col_bao2.slider("AP Scan Declination [deg]", -60.0, 60.0, 0.0, 5.0, help="Declination slice plane used for the 360° Right Ascension sky scan.")

        run_bao = st.button("🚀 Execute BAO Alcock-Paczynski Scan", type="primary", use_container_width=True)
        m_key = "mod_res_bao_warping"

        if run_bao:
            with st.spinner("Ray-tracing BAO sound horizon acoustic sphere warping..."):
                out_bao = run_bao_warping_analysis(
                    tree, gal_positions, nodes_per_galaxy,
                    path_depth_mpc=depth_bao_input,
                    sigma_mpc=sigma_smooth_input,
                    h_global=h_global_input,
                    dec_scan_deg=dec_scan_input
                )
                fig_bao, df_bao, bao_metrics = out_bao[0], out_bao[1], out_bao[2]
                st.session_state[m_key] = (fig_bao, df_bao, bao_metrics)

        if m_key in st.session_state:
            fig_bao, df_bao, bao_metrics = st.session_state[m_key]
            display_module_outputs(fig_bao, df_bao, bao_metrics, "bao_ruler_warping_alcock_paczynski.csv")
        else:
            st.info("⚙️ **Ready to Compute:** Select integration depth and scan declination, then click **🚀 Execute BAO Alcock-Paczynski Scan**.")


    # 10. CMB DIPOLE ANOMALY
    elif "CMB Dipole" in sub_3:
        st.header("Hybrid Kinematic & Unspooling CMB Dipole Anomaly")
        st.info(
            r"🛰️ **Physical Principle:** Combines Solar System Doppler motion ($369.0\text{ km/s}$) "
            r"with zero-parameter Substrate path-integrated metric unspooling anisotropy. "
            r"Integrated redshifting through 3D compiled matter density produces a zero-parameter "
            r"dipole axis tilt pointing toward the Centaurus / Great Attractor overdensity peak."
        )

        with st.expander("📖 **Interactive Parameter & Demonstration Guide**", expanded=False):
            st.markdown(r"""
            ### What This Module Demonstrates
            * **Hybrid CMB Dipole Signal:** Proves that the observed $3.35\text{ mK}$ CMB dipole is not a pure kinematic motion effect, but a hybrid composite combining Solar peculiar motion with directional metric unspooling along local matter fields.
            * **Dynamic Observable Axis Shift ($\sim 6.5^\circ\text{--}10^\circ$ Tilt):** Standard $\Lambda\text{CDM}$ assumes the CMB dipole is 100% kinematic Doppler motion pointing toward $(264.0^\circ, +48.0^\circ)$ (Leo/Crater). In Substrate Logistics, vector-adding the Solar Doppler dipole ($3.355\text{ mK}$) and the directional metric unspooling dipole ($0.45\text{--}0.80\text{ mK}$ integrated along 3D catalog density) shifts the net observable dipole axis to **$(273.2^\circ\text{--}276.0^\circ, +50.3^\circ)$**—inducing a dynamic angular tilt toward the Centaurus / Great Attractor overdensity peak.
            * **Resolving the NVSS / CatWISE Quasar Dipole Tension:** Explains why radio galaxy and quasar count surveys (e.g., Secrest et al.) measure dipole amplitudes over twice as large as expected from pure motion alone.

            ---
            ### 🔬 Key Physical Mechanics
            * **Thermodynamic Monopole Subtraction ($\bar{\Delta T}$):** Isolates pure directional sky fluctuations ($\ell \ge 1$) by subtracting the isotropic mean:
              $$\Delta T_{\text{sub}}(\hat{n}) = \Delta T_{\text{sub, raw}}(\hat{n}) - \bar{\Delta T}_{\text{sub, raw}}$$
            * **Adaptive Density Kernel:** Evaluates line-of-sight ray path integrals ($c \Delta z = \int \Delta H \, ds$) using adaptive spatial kernel smoothing ($\sigma(s) = \sigma_0 \sqrt{1 + s/10}$) to ensure continuous density integration without discrete galaxy sampling spikes.

            ---
            ### ⚙️ Interactive Controls & Parameter Reference
            * **HEALPix Resolution ($N_{\text{side}}$):** Spatial pixel density ($N_{\text{side}}=16 \to 3,072$ pixels; $N_{\text{side}}=32 \to 12,288$ pixels).
            * **Integration Depth [Mpc]:** Distance integrated along ray paths. Integrating out to $150\text{--}200\text{ Mpc}$ incorporates full macro-structures (Great Attractor, Shapley Concentration).
            * **Field Smoothing Scale ($\sigma$ in Sidebar):** Base Gaussian bandwidth kernel ($\sigma_0$). Values around $2.5\text{--}3.5\text{ Mpc}$ bridge sparse dust regions in 2M++ while maintaining filament density contrast.

            ---
            ### 📡 Survey Catalog Response & Performance Guidelines
            * **UNG 2013 (Local Volume $\le 35\text{ Mpc}$):** Evaluates local volume redshifting; renders in `< 0.8s`.
            * **2M++ Infrared ($\le 200\text{ Mpc}$, $\sigma = 3.5\text{ Mpc}$):** Full 3D ray-tracing across $12,288$ sky directions in `~3.85s`, yielding net composite dipole shifts of $\sim 6.5^\circ\text{--}10.2^\circ$.
            """)
            
        with st.expander(r"ℹ️ **Technical Note: Field Smoothing Scale ($\sigma$) & Performance**", expanded=False):
            st.markdown(r"""
            ### Impact of the Field Smoothing Scale ($\sigma$)

            The smoothing parameter $\sigma$ defines the 3D Gaussian kernel width used to construct continuous density fields $\rho_N(\vec{r})$ from discrete galaxy catalog positions. 

            #### 1. Computational Scaling ($O(\sigma^3)$ Spatial Volume)
            Neighbor galaxy queries search a spherical volume with radius $r_{\text{search}} = 3\sigma$. Because volume scales cubically ($V \propto \sigma^3$), increasing $\sigma$ drastically expands the number of catalog galaxies processed per ray-tracing point step:

            | Smoothing Scale ($\sigma$) | Search Radius ($3\sigma$) | Search Volume ($V$) | Relative Computational Work |
            | :--- | :--- | :--- | :--- |
            | **0.5 Mpc** | **1.5 Mpc** | **14.1 Mpc³** | **0.02×** (Ultra-fast) |
            | **1.2 Mpc** | **3.6 Mpc** | **195.4 Mpc³** | **0.30×** (Fast) |
            | **1.8 Mpc** *(Default)* | **5.4 Mpc** | **659.6 Mpc³** | **1.00×** (Baseline) |
            | **3.0 Mpc** | **9.0 Mpc** | **3,053.6 Mpc³** | **4.63×** (Moderate) |
            | **5.0 Mpc** | **15.0 Mpc** | **14,137.2 Mpc³** | **21.43×** (Slow) |
            | **10.0 Mpc** | **30.0 Mpc** | **113,097.3 Mpc³** | **171.47×** (Very Slow) |
            | **20.0 Mpc** | **60.0 Mpc** | **904,778.7 Mpc³** | **1,371.74×** (Execution Freeze) |

            > **Warning on High $\sigma$ Settings:** Settings of **$\sigma > 5.0\text{ Mpc}$** will perform hundreds of millions of floating-point operations across full-sky HEALPix rays ($12,288$ pixels $\times$ ray steps), causing execution times to scale into hours without throwing explicit errors.

            ---

            #### 2. Astrophysical Regime Recommendations

            * **$\sigma = 0.5\text{--}1.0\text{ Mpc}$ (Compact Structural Cores):**
            Use for high-resolution studies of compact cluster cores (e.g., Virgo, Fornax) or individual galactic rotation curves where sharp density gradients drive localized address garbage collection.
            * **$\sigma = 1.2\text{--}2.0\text{ Mpc}$ (Intergalactic Filaments & Default):**
            Optimal balance for Local Volume catalogs (UNG 2013, Cosmicflows-4). Maintains line-of-sight density continuity across intergalactic filaments while preventing artificial blurring into empty voids.
            * **$\sigma = 2.5\text{--}3.5\text{ Mpc}$ (Extended Surveys & Full-Sky Dipoles):**
            Recommended for deep infrared surveys (2M++) out to **200 Mpc** to bridge sparse catalog sampling gaps in galactic plane dust regions.
            """)

        col_cmb1, col_cmb2 = st.columns(2)
        nside_cmb = col_cmb1.selectbox("HEALPix Resolution (Nside)", [16, 32], index=1)
        depth_cmb = col_cmb2.slider("Integration Depth [Mpc]", 5.0, float(d_max_input), float(min(20.0, d_max_input)), 2.5)

        run_cmb = st.button("🚀 Execute CMB Dipole Anisotropy Shift", type="primary", use_container_width=True)
        m_key = "mod_res_cmb_dipole"

        if run_cmb:
            with st.spinner("Computing HEALPix hybrid kinematic + unspooling CMB dipole shift..."):
                out_cmb = run_cmb_dipole_analysis(
                    tree, gal_positions, nodes_per_galaxy, 
                    nside=nside_cmb, d_max_mpc=depth_cmb, sigma_mpc=sigma_smooth_input, h_global=h_global_input
                )
                fig_cmb, df_cmb, cmb_metrics = (out_cmb[0], out_cmb[1], out_cmb[2]) if len(out_cmb) == 3 else (out_cmb[0], None, out_cmb[1])
                st.session_state[m_key] = (fig_cmb, df_cmb, cmb_metrics)

        if m_key in st.session_state:
            fig_cmb, df_cmb, cmb_metrics = st.session_state[m_key]
            display_module_outputs(fig_cmb, df_cmb, cmb_metrics, "cmb_dipole_shift_lookup.csv")
        else:
            st.info("⚙️ **Ready to Compute:** Select HEALPix Nside and integration depth, then click **🚀 Execute CMB Dipole Anisotropy Shift**.")
            
    # 11. Gravitational Lensing
    elif "Gravitational Lensing" in sub_3:
        st.header("Zero-Dark-Matter Gravitational Lensing & NFW Halo Eradication")
        st.info(
            r"🔍 **Physical Principle:** In Substrate Logistics, optical light deflection does not rely on Cold Dark Matter halos. "
            r"Light deflection is decoupled into localized baryonic address deletion ($\theta_{\text{baryon}} \propto 1/b$) plus the linear metric unspooling tax ($\theta_{\text{scaffold}} = \frac{4 a_\Omega b}{c^2}$). "
            r"Standard cosmology mistakes the linear unspooling tax for an invisible quadratic dark matter halo ($M_{\text{DM}} \propto b^2$). "
            r"The module dynamically evaluates zero-parameter physical crossover horizons ($b_{\text{cross}}$ and $b_{1:1}$) where optical deflection parity and $1:1$ mass equivalence are reached."
        )
        
        with st.expander("📖 **Interactive Parameter & Demonstration Guide**", expanded=False):
            st.markdown(r"""
            ### What This Module Demonstrates
            * **Zero-Dark-Matter Optical Deflection:** Bends light rays using only visible baryonic mass ($M_{\text{baryon}}$) plus a linear metric unspooling tax ($\theta_{\text{scaffold}} \propto b$), eliminating particle dark matter halos.
            * **Eradication of NFW Dark Matter Halos:** Proves that standard cosmology misinterprets the linear unspooling tax as an artificial quadratic dark matter halo ($M_{\text{DM}} = \frac{a_\Omega b^2}{G} \propto b^2$).
            * **Zero-Parameter Physical Horizons:** Dynamically calculates $b_{\text{cross}}$ ($50\%$ deflection parity) and $b_{1:1}$ ($100\%$ total baryonic mass equivalence) for any target galaxy or cluster.

            ---
            ### 🔬 Key Physical Mechanics
            * **Dual-Component Deflection Angle:**
              $$\theta_{\text{total}}(b) = \frac{4 G M_{\text{enc}}(b)}{c^2 b} + \frac{4 a_\Omega b}{c^2}$$
            * **Projected Hernquist Mass Profile:** Models extended cluster/galaxy mass distributions without central singularities:
              $$M_{\text{enc}}(b) = M_{\text{baryon}} \left[ \frac{b^2}{(b + a_{\text{scale}})^2} \right]$$
            * **1:1 Enclosed Mass Horizon ($b_{\text{cross}}$):** Exact crossover radius where $\theta_{\text{baryon}} = \theta_{\text{scaffold}}$:
              $$b_{\text{cross}} = \sqrt{\frac{M_{\text{baryon}} \mathcal{K}_\Omega M_\odot}{a_\Omega m_p}} - a_{\text{scale}}$$

            ---
            ### ⚙️ Interactive Controls & Parameter Reference
            * **Lensing Mode Selection:** Toggle between preset benchmark targets (e.g., Abell 1689, Milky Way) and live 3D catalog query mode.
            * **Target Galaxy Dropdown (Catalog Mode):** Select any galaxy from the active survey catalog to automatically sample its 3D position $(X, Y, Z)$ and compiled baryonic mass.
            * **Baryonic Mass Baseline:** Adjust the visible baryonic mass (stars + gas) in Solar units ($M_\odot$).
            * **Kernel Scale ($\sigma$):** Adjust Gaussian smoothing bandwidth for local node density ($\rho_N$) and $H_{\text{local}}$ evaluation.

            ---
            ### 📡 Survey Catalog Response & Performance Guidelines
            * **UNG 2013 (Local Volume):** Ideal for analyzing nearby dwarf and spiral galaxy footprints ($d \le 35\text{ Mpc}$).
            * **Cosmicflows-4 & 2M++:** Deep survey catalogs tracing massive cluster centers out to $150\text{--}200\text{ Mpc}$, where local density boosts $H_{\text{local}}$ and dynamic baseline acceleration $a_\Omega$.
            """)

        # Mode Selection: Presets vs Live Catalog Query
        mode_lensing = st.radio(
            "Select Lensing Mode:",
            ["🎯 Benchmark Target Presets (Deep-Sky Lenses)", "📡 Live Catalog Query (3D Density-Coupled a_Ω)"],
            horizontal=True
        )

        if "Benchmark" in mode_lensing:
            col_l1, col_l2 = st.columns(2)
            target_preset = col_l1.selectbox(
                "Select Lens Target Preset:",
                list(LENS_PRESETS.keys())
            )
            
            default_m_solar_12 = LENS_PRESETS[target_preset]["M_baryon_solar"] / 1e12
            step_val = 0.05 if default_m_solar_12 < 1.0 else 0.5
            fmt_str = "%.2f" if default_m_solar_12 < 1.0 else "%.1f"

            custom_m_scaled = col_l2.number_input(
                "True Baryonic Mass Baseline [10¹² M☉]:",
                value=float(default_m_solar_12),
                step=step_val,
                format=fmt_str,
                help="Strictly the visible baryonic mass (stars + gas). No dark matter halos!"
            )
            custom_m_solar = custom_m_scaled * 1e12
            use_cat = False
            target_pos = None
            sigma_k = 1.8

        else: # Live Catalog Query Mode
            st.markdown("### 📡 Live 3D Catalog Density Sampling")
            st.info(
                "💡 **How Catalog Lensing Works:** This mode queries the active 3D `cKDTree` galaxy dataset loaded above. "
                "Selecting a galaxy automatically fetches its catalog name, spatial coordinates, and compiled baryonic mass from luminosity."
            )

            # Infobox for Kernel Scale
            with st.expander("ℹ️ Kernel Scale Smoothing (σ) Infobox & Guidelines", expanded=False):
                st.markdown(
                    "**What is the Kernel Scale (σ)?**\n"
                    "The Kernel Scale $\\sigma$ defines the 3D Gaussian smoothing window used to calculate spatial node density $\\rho_N(\\vec{r})$.\n\n"
                    "* **$\\sigma = 0.5 - 1.0 \\text{ Mpc}$ (Compact Group Scale):** Isolates tight, localized mass concentrations. High sensitivity to individual neighbor galaxies.\n"
                    "* **$\\sigma = 1.8 \\text{ Mpc}$ (Standard Cluster Default):** Balances point-like galaxy spikes while smoothing over cluster-scale matter distributions.\n"
                    "* **$\\sigma = 3.0 - 5.0 \\text{ Mpc}$ (Macro Filament / Supercluster Scale):** Smoothes density across large cosmic web structures, ideal for regional background $H_{\\text{local}}$ estimations.\n"
                )

            # Build adaptive dropdown labels for each galaxy in the catalog
            gal_options = []
            for i in range(len(gal_names)):
                d_mpc = float(np.linalg.norm(gal_positions[i]))
                m_solar = float((nodes_per_galaxy[i] * M_PROTON) / M_SOLAR_KG)
                
                # Format mass adaptively so it never displays as 0.000
                if m_solar >= 1e11:
                    m_label = f"{m_solar / 1e12:.3f} × 10¹² M☉"
                elif m_solar >= 1e8:
                    m_label = f"{m_solar / 1e9:.1f} × 10⁹ M☉"
                else:
                    m_label = f"{m_solar / 1e6:.1f} × 10⁶ M☉"
                    
                label = f"#{i}: {gal_names[i]} (d={d_mpc:.1f} Mpc | Mass={m_label})"
                gal_options.append(label)

            col_c1, col_c2, col_c3 = st.columns([2.2, 1.0, 1.2])
            
            # Galaxy Selectbox dropdown
            selected_option = col_c1.selectbox(
                "Select Target Galaxy from Catalog:",
                options=gal_options,
                index=min(10, len(gal_options)-1)
            )
            
            # Extract target index from option string
            gal_idx_target = int(selected_option.split(":")[0].replace("#", ""))
            target_pos = gal_positions[gal_idx_target]
            target_name = gal_names[gal_idx_target]
            
            # Derive target galaxy's baryonic mass in Solar Masses
            catalog_m_solar = float((nodes_per_galaxy[gal_idx_target] * M_PROTON) / M_SOLAR_KG)

            col_c1.caption(f"📍 **Target Position:** X={target_pos[0]:.2f}, Y={target_pos[1]:.2f}, Z={target_pos[2]:.2f} Mpc | **Name:** {target_name}")

            sigma_k = col_c2.slider(
                "Kernel Scale σ [Mpc]:",
                min_value=0.5,
                max_value=5.0,
                value=1.8,
                step=0.1,
                help="Gaussian smoothing window for local density estimation."
            )

            # Adaptive Central Mass Input Scale
            if catalog_m_solar >= 1e11:
                input_val = catalog_m_solar / 1e12
                unit_label = "Compiled Mass [10¹² M☉]:"
                fmt_str = "%.3f"
                step_val = 0.01
            elif catalog_m_solar >= 1e8:
                input_val = catalog_m_solar / 1e9
                unit_label = "Compiled Mass [10⁹ M☉]:"
                fmt_str = "%.1f"
                step_val = 1.0
            else:
                input_val = catalog_m_solar / 1e6
                unit_label = "Compiled Mass [10⁶ M☉]:"
                fmt_str = "%.1f"
                step_val = 1.0

            custom_m_scaled = col_c3.number_input(
                unit_label,
                value=input_val,
                step=step_val,
                format=fmt_str,
                help="Auto-populated with target galaxy catalog mass. You can override if desired."
            )

            # Convert user input back to Solar Masses for profile calculation
            if "10¹²" in unit_label:
                custom_m_solar = custom_m_scaled * 1e12
            elif "10⁹" in unit_label:
                custom_m_solar = custom_m_scaled * 1e9
            else:
                custom_m_solar = custom_m_scaled * 1e6

            use_cat = True
            target_preset = f"Catalog: {target_name}"

        # Execute Button
        run_lens = st.button("🚀 Execute Gravitational Lensing Analysis", type="primary", use_container_width=True)
        m_key = "mod_res_lensing"

        if run_lens:
            with st.spinner("Computing zero-G optical deflection profile..."):
                fig_lens, df_lens, lens_metrics = run_lensing_analysis(
                    lens_target_name=target_preset,
                    custom_M_baryon=custom_m_solar,
                    h_global=h_global_input,
                    use_catalog_query=use_cat,
                    target_pos_mpc=target_pos,
                    sigma_mpc=sigma_k,
                    tree=tree,
                    gal_positions=gal_positions,
                    nodes_per_galaxy=nodes_per_galaxy
                )
                st.session_state[m_key] = (fig_lens, df_lens, lens_metrics)

        if m_key in st.session_state:
            fig_lens, df_lens, lens_metrics = st.session_state[m_key]
            display_module_outputs(fig_lens, df_lens, lens_metrics, "substrate_gravitational_lensing.csv")
        else:
            st.info("⚙️ **Ready to Compute:** Select lensing mode above and click **🚀 Execute Gravitational Lensing Analysis**.")

# ------------------------------------------------------------------------------
# DOMAIN 4: MASTER REGULATORY PROTOCOL & BANDWIDTH ALLOCATION
# ------------------------------------------------------------------------------
with domain_4:
    sub_4 = st.radio(
        "Select Module:",
        [
            "🗑️ Spatial Address Garbage Collection Engine", 
            "⏱️ CPU Bandwidth Throttling & Clock Dilation"
        ],
        horizontal=True
    )
    st.divider()

    if "Garbage Collection" in sub_4:
        st.header("Spatial Address Garbage Collection Engine (E = 0 Memory Purge)")
        st.info(
            r"🖥️ **Physical Principle:** Gravity is reclassified as the zero-energy background address garbage collection "
            r"protocol ($E=0$) of the spatial grid. Local Informational Loads ($N = M/m_p$) execute spatial deletion "
            r"at rate $g_\Omega = N \mathcal{K}_\Omega / R^2$, pulling spatial coordinate addresses inward at $v_{\text{inflow}} = \sqrt{2g_\Omega R}$. "
            r"At the Event Horizon ($R_{\text{EH}} = 2N\bar{\lambda}_p / \mathcal{T}_\Omega^2$), inflow hits the absolute clock limit $1.0c$, triggering a Kernel Lock."
        )

        with st.expander("📖 **Interactive Parameter & Demonstration Guide**", expanded=False):
            st.markdown(r"""
            ### What This Module Demonstrates
            * **Zero-Energy Memory Purge ($E = 0$):** Demonstrates that gravity is zero-energy spatial address deletion. Empty coordinate pointers are purged from the metric index without emitting messenger particles, eliminating the graviton category error.
            * **Discrete Pixel Purge Rate ($\dot{N}_{\text{pixels}}$):** Evaluates exact sub-node voxel deletions ($V_{\text{node}} = \bar{\lambda}_p^3 = 9.301 \times 10^{-48}\text{ m}^3$) per second across boundary surfaces (e.g., $6.1328 \times 10^{65}\text{ px/s}$ for Earth).
            * **Coordinate Inflow & Escape Velocity Equivalence ($v_{\text{inflow}} = v_{\text{escape}}$):** Proves that spatial address inflow velocity $v_{\text{inflow}} = \sqrt{2 g_\Omega R}$ is physically identical to classical escape velocity, while low-orbit bypass speed evaluates to $v_{\text{orbit}} = \sqrt{g_\Omega R} = v_{\text{escape}} / \sqrt{2}$.
            * **Kernel Lock Safety Margin & $1.0c$ Horizon Limit:** Calculates remaining processing headroom $\max(0, 1.0 - v/c) \times 100\%$. Sub-light objects preserve wide safety margins (e.g., $99.9963\%$ for Earth at $0.0037\%$ capacity), whereas black hole horizons lock down at $100\%$ capacity ($0.0000\%$ margin) right at $1.0c$.

            ---
            ### 🔬 Key Physical Mechanics
            * **Spatial Address Deletion Rate ($g_\Omega$):**
              $$g_\Omega(r) = \frac{N \mathcal{K}_\Omega}{r^2}$$
            * **Coordinate Inflow Speed ($v_{\text{inflow}}$):**
              $$v_{\text{inflow}}(r) = \sqrt{2 g_\Omega(r) r} = \sqrt{\frac{2 N \mathcal{K}_\Omega}{r}}$$
            * **Discrete Voxel Purge Intensity ($\dot{N}_{\text{pixels}}$):**
              $$\dot{N}_{\text{pixels}}(r) = \frac{4 \pi r^2 v_{\text{inflow}}(r)}{\bar{\lambda}_p^3}$$

            ---
            ### ⚙️️ Interactive Controls & Parameter Reference
            * **Select Celestial Preset or Custom Load:** Choose pre-configured hardware targets (Earth, Moon, Sun, Cygnus X-1, Sagittarius A*, M87*) or configure custom mass and surface radius targets.
            * **Evaluation Radial Extent ($r / R_{\text{surface}}$):** Sets the outward multiplier for radial profile evaluation (default $10.0\times$).
            * **Custom Mass & Radius Inputs:** Enabled when selecting "Custom Mass & Radius Target". Allows arbitrary $M$ (kg) and $R$ (m) evaluation.

            ---
            ### 📡 Celestial Benchmark Presets & Behavioral Guidelines
            * **Terrestrial & Solar Targets (Earth, Moon, Sun):** Feature low capacity utilization ($< 0.01\%$) and high safety margins ($> 99.99\%$). Display inline safety badges.
            * **Black Hole Targets (Cygnus X-1, Sgr A*, M87*):** Reach $100\%$ capacity ($v_{\text{inflow}} = 1.0c$) at boundary radius $R_{\text{lock}}$, triggering Kernel Lock state visualization and horizon boundary markers.
            """)

        col_gc1, col_gc2 = st.columns(2)
        preset_choice = col_gc1.selectbox(
            "Select Celestial Preset or Custom Load:",
            list(CELESTIAL_PRESETS.keys()) + ["Custom Mass & Radius Target"]
        )
        r_max_ratio_input = col_gc2.slider(
            "Evaluation Radial Extent (r / R_surface)", 
            2.0, 50.0, 10.0, 1.0,
            help="Multiplier setting how far outward in units of surface radius to calculate the inflow field."
        )

        custom_m, custom_r = None, None
        if preset_choice == "Custom Mass & Radius Target":
            c_col1, c_col2 = st.columns(2)
            custom_m = c_col1.number_input("Custom Mass [kg]", value=5.9722e24, format="%.3e")
            custom_r = c_col2.number_input("Custom Surface Radius [m]", value=6.3712e6, format="%.3e")

        run_gc = st.button("🚀 Execute Garbage Collection Analysis", type="primary", use_container_width=True)
        m_key = "mod_res_garbage_collection"

        if run_gc:
            with st.spinner("Computing spatial address deletion rates and inflow velocity profiles..."):
                fig_gc, df_gc, gc_metrics = run_spatial_garbage_collection_analysis(
                    preset_key=preset_choice if preset_choice in CELESTIAL_PRESETS else "Custom",
                    custom_mass_kg=custom_m,
                    custom_radius_m=custom_r,
                    r_max_ratio=r_max_ratio_input
                )
                st.session_state[m_key] = (fig_gc, df_gc, gc_metrics)

        if m_key in st.session_state:
            fig_gc, df_gc, gc_metrics = st.session_state[m_key]
            display_module_outputs(
                fig_gc, df_gc, gc_metrics, 
                "spatial_garbage_collection_profile.csv"
            )
        else:
            st.info("⚙️ **Ready to Compute:** Select a celestial benchmark above and click **🚀 Execute Garbage Collection Analysis**.")

    elif "CPU Bandwidth" in sub_4:
        st.header("CPU Bandwidth Throttling & Gravitational Time Dilation")
        st.info(
            r"⏱️ **Physical Principle:** Time is localized metric processing latency. "
            r"Governed by the Pythagorean Bandwidth Allocation ($c^2 = v^2 + u^2$), external kinematic movement ($v_\Omega^2$) "
            r"and gravitational address deletion ($g_\Omega^2$) consume finite metric bandwidth. "
            r"The residual processing capacity available for atomic clock transitions evaluates strictly from dynamic loads to "
            r"$u_{\text{clock}} = \sqrt{1.0 - (g_\Omega^2 + v_\Omega^2)}$. "
            r"At the Event Horizon ($R_{\text{EH}}$), $g_\Omega^2 = 1.0$, dropping residual clock update capacity to $u_{\text{clock}} = 0.0000$ ($0\%$ CPU) "
            r"and triggering a 100% Kernel Lock."
        )

        with st.expander("📖 **Interactive Parameter & Demonstration Guide**", expanded=False):
            st.markdown(r"""
            ### What This Module Demonstrates
            * **CPU Throttling vs. Time Curvature:** Proves that atomic clock dephasing is literal hardware CPU throttling under heavy computational load ($g_\Omega^2 + v_\Omega^2$), eliminating continuous spacetime warping.
            * **Reference Frame Disambiguation:** Disambiguates absolute frequency shifts vs. unthrottled deep space ($\Delta f_{\text{space}} = -2.5040 \times 10^{-10}$) from relative frequency shifts observed by Earth ground observers ($\Delta f_{\text{Earth}} = +4.4661 \times 10^{-10}$).
            * **GPS Satellite Clock Calibration:** Demonstrates why GPS satellite atomic clocks run faster ($+38.59\,\mu\text{s/day}$ / $+4.4661 \times 10^{-10}$) relative to Earth ground clocks due to weaker gravitational deletion loads at $20,180\text{ km}$ altitude.
            * **ISS Orbit Clock Calibration:** Demonstrates why ISS atomic clocks run slower ($-24.62\,\mu\text{s/day}$) relative to Earth ground clocks due to high LEO orbital velocity ($7.67\text{ km/s}$).
            * **Event Horizon 100% CPU Freeze ($u_{\text{clock}} = 0$):** Demonstrates that at $R_{\text{EH}}$, $100\%$ of processing bandwidth is dedicated to spatial address deletion ($g_\Omega^2 = 1.0$), leaving zero residual capacity for atomic state updates ($u_{\text{clock}} = 0.0000$).

            ---
            ### 🔬 Key Physical Mechanics
            * **Pythagorean Bandwidth Allocation:**
              $$c^2 = v^2 + u^2 \implies u = c \sqrt{1 - \frac{v^2}{c^2}}$$
            * **Absolute Shift vs. Deep Space ($u = 1.0$):**
              $$\left(\frac{\Delta f}{f_0}\right)_{\text{space}} = u_{\text{clock}} - 1.0$$
            * **Relative Shift vs. Earth Ground ($u_{\text{Earth}}$):**
              $$\left(\frac{\Delta f}{f_0}\right)_{\text{Earth}} = \frac{u_{\text{target}} - u_{\text{Earth}}}{u_{\text{Earth}}}$$

            ---
            ### ⚙ Interactive Controls & Parameter Reference
            * **Select Celestial Benchmark or Custom Orbit:** Choose pre-configured hardware targets (Earth Surface, GPS Orbit, ISS Orbit, Solar Photosphere, Sirius B, Crab Pulsar, Sagittarius A*) or configure custom mass, radius, and velocity.
            * **Observer Kinematic Mode:** Toggle between `Static (v = 0)` (pure gravitational load) and `Orbital (v = v_orbit)` (includes circular orbital routing load).
            * **Internal Cell Tension ($\alpha_\Omega$):** Adjust static structural cell tension load (fine structure limit). Displays as internal cell lock bandwidth in Panel 1 without corrupting external clock drift metrics.

            ---
            ### 📡 Celestial Benchmark Presets & Hardware Guidelines
            * **Earth Surface & Satellite Constellations (GPS, ISS):** Captures high-precision microsecond daily clock drift differences ($+38.59\,\mu\text{s/day}$ GPS, $-24.62\,\mu\text{s/day}$ ISS vs ground).
            * **Relativistic Stellar Remnants (Sirius B, Crab Pulsar):** Demonstrates severe gravitational clock dephasing ($-22.39\text{ s/day}$ on Sirius B, $u_{\text{clock}} = 0.8097$ on Crab Pulsar).
            * **Event Horizon (Sagittarius A*):** Triggers $100\%$ CPU Kernel Lock ($u_{\text{clock}} = 0.0000$), freezing atomic state updates completely at $R_{\text{EH}}$.
            """)

        col_tb1, col_tb2, col_tb3 = st.columns(3)
        preset_tb = col_tb1.selectbox(
            "Select Celestial Benchmark or Custom Orbit:",
            list(THROTTLING_PRESETS.keys()) + ["Custom Mass, Radius & Velocity"]
        )
        obs_mode_tb = col_tb2.selectbox(
            "Observer Kinematic Mode:",
            ["Static (v = 0)", "Orbital (v = v_orbit)"],
            help="Static evaluates pure gravitational deletion dephasing; Orbital includes circular orbital speed routing load."
        )
        alpha_input = col_tb3.slider(
            "Internal Cell Tension alpha_Omega", 
            0.0, 0.1, 0.0, 0.005,
            help="Static internal structural cell tension load."
        )

        custom_m, custom_r, custom_v = None, None, None
        if preset_tb == "Custom Mass, Radius & Velocity":
            c_col1, c_col2, c_col3 = st.columns(3)
            custom_m = c_col1.number_input("Custom Mass [kg]", value=5.9722e24, format="%.3e")
            custom_r = c_col2.number_input("Custom Radius [m]", value=6.3712e6, format="%.3e")
            custom_v = c_col3.number_input("Custom Velocity [m/s]", value=0.0, step=100.0)

        run_tb = st.button("🚀 Execute CPU Bandwidth Throttling Analysis", type="primary", use_container_width=True)
        m_key = "mod_res_cpu_bandwidth"

        if run_tb:
            with st.spinner("Calculating metric processing capacity and clock frequency shift..."):
                fig_tb, df_tb, tb_metrics = run_cpu_bandwidth_throttling_analysis(
                    preset_key=preset_tb if preset_tb in THROTTLING_PRESETS else "Custom",
                    custom_mass_kg=custom_m,
                    custom_radius_m=custom_r,
                    custom_velocity_m_s=custom_v,
                    alpha_tension=alpha_input,
                    observer_mode=obs_mode_tb
                )
                st.session_state[m_key] = (fig_tb, df_tb, tb_metrics)

        if m_key in st.session_state:
            fig_tb, df_tb, tb_metrics = st.session_state[m_key]
            display_module_outputs(
                fig_tb, df_tb, tb_metrics, 
                "cpu_bandwidth_throttling_profile.csv"
            )
        else:
            st.info("⚙️ **Ready to Compute:** Select a benchmark above and click **🚀 Execute CPU Bandwidth Throttling Analysis**.")