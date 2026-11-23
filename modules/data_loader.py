import os
import numpy as np
import pandas as pd
import streamlit as st
from scipy.spatial import cKDTree
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.vizier import Vizier

# ==========================================
# 1. EXACT FIRST-PRINCIPLES HARDWARE INVARIANTS (ZERO CODATA G)
# ==========================================
C_LIGHT          = 299792458.0              # Exact SI speed of light [m/s]
K_BOLTZMANN      = 1.380649e-23             # Exact SI Boltzmann constant [J/K]
M_PROTON         = 1.67262192595e-27        # Proton mass anchor [kg] (CODATA 2022)
M_PROTON_MEV     = 938.27208816             # Proton rest energy [MeV/c^2] (CODATA 2022)
LAMBDA_BAR_P     = 2.10308910051e-16        # Reduced proton Compton wavelength [m] (CODATA 2022)
V_NODE           = LAMBDA_BAR_P ** 3        # Sub-node spatial voxel volume [m^3] (~9.301929e-48)
MEV_FM3_TO_PA    = 1.602176634e32           # 1 MeV/fm^3 in Pascals [J/m^3]

# 120-Cell Poincare Dodecahedral Geometry (kappa = |I*| * E_f / F = 120 * 5 / 12 = 50)
KAPPA            = 50.0
KAPPA_5          = KAPPA ** 5               # 3.125e8 (5D state-space partition)
A_OMEGA_TENSOR   = KAPPA ** 10              # 9.765625e16 (10D metric attenuation)
B_MAX            = 400.0 / 3.0              # Bare-metal buffer scalar (8/3 * kappa)
GEOM_SCALAR      = 1.2                      # 36 routing channels / 30 dodecahedral edges
DRAG_TAX_3D      = 0.017032                 # 3D volumetric 36-channel routing drag tax (1.7032%)

# Exact Algebraic Universal Topological Drag & Closed-Form T_Omega / K_Omega
C_DELTA                 = (np.sqrt(5.0) - 2.0) / 500.0        # 4.7213595500e-4
SQRT_C_DELTA            = np.sqrt(C_DELTA)                    # 0.0217286897
C_DELTA_FOURTH_ROOT     = np.sqrt(SQRT_C_DELTA)               # 0.1474065
C_DELTA_INV_FOURTH_ROOT = 1.0 / C_DELTA_FOURTH_ROOT           # 6.7839605
DISK_TRUNCATION_RATIO   = 0.5 * np.log(1.0 / C_DELTA)         # 3.829119 (2D exponential disk cliff)

T_OMEGA          = (B_MAX * A_OMEGA_TENSOR) / (1.0 + C_DELTA) # 1.3014688631e19
T_OMEGA_SQ       = T_OMEGA ** 2                               # 1.6938212016e38
K_OMEGA          = (LAMBDA_BAR_P * (C_LIGHT ** 2)) / T_OMEGA_SQ # 1.1159160239e-37 m^3/(node*s^2)

# Astronomical Unit Conversions
M_SOLAR_KG       = 1.98847e30
AU_TO_METER      = 1.495978707e11
MPC_TO_METER     = 3.0856775814913673e22
KPC_TO_METER     = MPC_TO_METER / 1000.0
UNIT_CONV        = 1000.0 / MPC_TO_METER                      # Exact s^-1 per (km/s/Mpc) (~3.24077929e-20)
RAD_TO_ARCSEC    = (180.0 * 3600.0) / np.pi

# Dual Energy-Density Ledgers & Topological Drag Pressure Floor (P_delta)
EPSILON_MAX_PA   = (M_PROTON * (C_LIGHT ** 2)) / V_NODE       # 1.61609e37 Pa (100.869 GeV/fm^3)
P_DELTA_MICRO_PA = C_DELTA * EPSILON_MAX_PA                   # 7.6302e33 Pa
B_QCD_MEV_FM3    = P_DELTA_MICRO_PA / MEV_FM3_TO_PA           # 47.624 MeV/fm^3
R_YUKAWA_FM      = (LAMBDA_BAR_P * 1.0e15) * C_DELTA_INV_FOURTH_ROOT # 1.42673 fm
M_PION_MEV       = M_PROTON_MEV * C_DELTA_FOURTH_ROOT         # 138.307 MeV/c^2

# Sub-Node Bohr / Electromagnetic Phase-Lock Constants
ALPHA_INV             = 137.035999084          # Inverse Fine-Structure / Golden Routing Angle
PROTON_ELECTRON_RATIO = 1836.15267343          # 1836-state saturated core mass ratio (m_p / m_e)
E_RYD_J               = (M_PROTON * (C_LIGHT ** 2)) / (2.0 * PROTON_ELECTRON_RATIO * (ALPHA_INV ** 2)) # 13.6057 eV
T_LOCK_RYD_K          = (2.0 * E_RYD_J) / (3.0 * K_BOLTZMANN)  # 105,258 K thermal unlock temperature

SIGMA_0_RAW      = (2.0 * M_PROTON * (C_LIGHT ** 2)) / (V_NODE * T_OMEGA_SQ * KAPPA_5)
SIGMA_0_LEDGER   = SIGMA_0_RAW * (1.0 - DRAG_TAX_3D)          # 6.002304e-10 J/m^3
P_DELTA_MACRO_PA = C_DELTA * SIGMA_0_LEDGER                   # 2.833903e-13 Pa (0.2834 pPa)
P_DELTA_K_CM3    = (P_DELTA_MACRO_PA / K_BOLTZMANN) * 1.0e-6  # 2.0526e4 K cm^-3

# Zero-Sum 120-Cell Matter Partitions, H_global, S8 & Cluster Mass Bias
OMEGA_M_VOID     = (12.0 / KAPPA) - SQRT_C_DELTA              # 0.218271
OMEGA_L_VOID     = 1.0 - OMEGA_M_VOID                         # 0.781729
OMEGA_M_GLOBAL   = (12.0 / KAPPA) + GEOM_SCALAR * SQRT_C_DELTA # 0.266074
OMEGA_L_GLOBAL   = 1.0 - OMEGA_M_GLOBAL                       # 0.733926
OMEGA_M_CORE     = 0.2840                                     # Compiled cluster core saturation
SIGMA_8_VOID     = np.sqrt(OMEGA_L_VOID)                      # 0.884154
S8_VOID          = SIGMA_8_VOID * np.sqrt(OMEGA_M_VOID / 0.300)   # 0.7542
S8_GLOBAL        = SIGMA_8_VOID * np.sqrt(OMEGA_M_GLOBAL / 0.300) # 0.8327
S8_CORE          = SIGMA_8_VOID * np.sqrt(OMEGA_M_CORE / 0.300)   # 0.8603
CLUSTER_BIAS_1MB = OMEGA_M_VOID / OMEGA_M_GLOBAL              # 0.8203 (17.97% lensing deficit)

# Derived Cosmological Expansion & BAO Horizon
H_LAMBDA_SI      = np.sqrt((8.0 * np.pi * LAMBDA_BAR_P * SIGMA_0_LEDGER) / (3.0 * M_PROTON * T_OMEGA_SQ))
H_LAMBDA         = float(H_LAMBDA_SI / UNIT_CONV)             # 59.6163 km/s/Mpc
H_GLOBAL         = float(H_LAMBDA / np.sqrt(OMEGA_L_VOID))    # 67.4275 km/s/Mpc
H_OMEGA_SI       = H_GLOBAL * UNIT_CONV
A_OMEGA          = (C_LIGHT * H_OMEGA_SI) / (2.0 * np.pi)     # 1.04262e-10 m/s^2
RHO_GC_EXACT     = (SIGMA_0_LEDGER / (C_LIGHT ** 2)) / M_PROTON # 3.9928 nodes/m^3
RHO_GC_THRESHOLD = 4.0                                        # Integer 4-channel He^4 ceiling
R_BAO_MPC        = (C_LIGHT * T_OMEGA * (3.0 * C_DELTA) / GEOM_SCALAR) / MPC_TO_METER # 149.2489 Mpc

# Photometric Mass-to-Light & Absolute Solar Magnitudes
M_B_SUN            = 5.44              # Solar absolute B-band magnitude (UNG 2013 / Cosmicflows)
M_K_SUN            = 3.27              # Solar absolute K-band magnitude (2M++ Infrared)
M_W1_SUN           = 3.24              # Solar absolute WISE W1 (3.4 um) Vega magnitude
M_L_B              = 0.44              # B-band mass-to-light ratio (Bell & de Jong 2001)
M_L_K              = 0.60              # K-band mass-to-light ratio (Near-Infrared)
M_L_W1             = 0.50              # WISE W1 stellar M/L (McGaugh & Schombert 2014 / Lelli+ 2016)
WISE_W1_AB_TO_VEGA = 2.699             # Jarrett et al. (2011) W1_Vega = W1_AB - 2.699

MASS_SOURCE_LABELS = {
    1: "K-band (McGaugh+14)",
    2: "WISE W1 (Lelli+16)",
    3: "B-V Color (Bell+03)",
    4: "B/Single-Band Fallback",
}

# Exact J2000 Galactic -> ICRS orthogonal rotation matrix (columns are Galactic x, y, z unit axes in ICRS)
R_GAL_TO_ICRS = SkyCoord(
    l=[0.0, 90.0, 0.0] * u.deg, b=[0.0, 0.0, 90.0] * u.deg, frame='galactic'
).icrs.cartesian.xyz.value

def galactic_to_icrs_vecs(vecs_gal):
    """Rotates (N, 3) Cartesian vectors from Galactic frame into the ICRS frame of gal_positions."""
    return np.asarray(vecs_gal, dtype=float) @ R_GAL_TO_ICRS.T

# ==========================================
# HELPER PARSERS FOR MULTI-CATALOG SCHEMAS
# ==========================================
def format_baryonic_mass(m_solar_val):
    """Formats baryonic mass adaptively based on astronomical scale."""
    if m_solar_val >= 1e11:
        return f"{m_solar_val / 1e12:.3f} × 10¹² M☉"
    elif m_solar_val >= 1e8:
        return f"{m_solar_val / 1e9:.1f} × 10⁹ M☉"
    else:
        return f"{m_solar_val / 1e6:.1f} × 10⁶ M☉"

def extract_numeric_column(cat_table, candidates):
    """Safely extracts a numeric numpy array from candidate column names with type-safe mask filling."""
    col_names = cat_table.colnames
    for cand in candidates:
        matched = next((c for c in col_names if c.lower() == cand.lower()), None)
        if matched is not None:
            col_data = cat_table[matched]
            if hasattr(col_data, 'filled'):
                try:
                    vals = col_data.filled(np.nan)
                except (TypeError, ValueError):
                    vals = col_data.astype(float).filled(np.nan)
            else:
                vals = col_data
            num_vals = pd.to_numeric(vals, errors='coerce')
            arr_vals = np.array(num_vals, dtype=float)
            if not np.all(np.isnan(arr_vals)):
                return arr_vals
    return None

def extract_string_column(cat_table, candidates):
    """Safely extracts galaxy names/identifiers, handling string/int masked arrays and decoding byte strings."""
    col_names = cat_table.colnames
    for cand in candidates:
        matched = next((c for c in col_names if c.lower() == cand.lower()), None)
        if matched is not None:
            col_data = cat_table[matched]
            if hasattr(col_data, 'filled'):
                try:
                    vals = col_data.filled('')
                except (TypeError, ValueError):
                    vals = col_data.astype(str).filled('')
            else:
                vals = col_data
            
            clean_str = []
            for i, item in enumerate(vals):
                if hasattr(item, 'decode'):
                    s = item.decode('utf-8', errors='ignore').strip()
                else:
                    s = str(item).strip()
                if s.startswith("b'") and s.endswith("'"):
                    s = s[2:-1].strip()
                if s.startswith('b"') and s.endswith('"'):
                    s = s[2:-1].strip()
                if s in ['', 'nan', 'None', '--', 'null', 'masked']:
                    s = f"Galaxy #{i}"
                clean_str.append(s)
                
            arr_str = np.array(clean_str, dtype=str)
            if not np.all([s.startswith("Galaxy #") for s in arr_str]):
                return arr_str
    return None

def parse_sky_coordinates(cat_table):
    """Robustly extracts RA/Dec or Galactic coordinates across UNG, Cosmicflows-4, and 2M++."""
    col_names = cat_table.colnames
    coord_pairs = [
        (['_RAJ2000', '_RA', 'RA2000', 'RAJ2000', 'RAdeg', 'RA_deg', 'RA', 'ra'],
         ['_DEJ2000', '_DE', 'DEC2000', 'DEJ2000', 'DEdeg', 'DEC_deg', 'DEC', 'dec', 'DE'], False),
        (['_l', 'GLON', 'glon', 'l_deg', 'l'],
         ['_b', 'GLAT', 'glat', 'b_deg', 'b'], True)
    ]

    for ra_cands, dec_cands, is_galactic in coord_pairs:
        ra_col_name = next((c for c in col_names if c.lower() in [x.lower() for x in ra_cands]), None)
        dec_col_name = next((c for c in col_names if c.lower() in [x.lower() for x in dec_cands]), None)

        if ra_col_name is not None and dec_col_name is not None:
            ra_raw = cat_table[ra_col_name]
            dec_raw = cat_table[dec_col_name]

            if hasattr(ra_raw, 'filled'):
                try:
                    ra_raw = ra_raw.filled(np.nan)
                except (TypeError, ValueError):
                    ra_raw = ra_raw.astype(float).filled(np.nan)
            if hasattr(dec_raw, 'filled'):
                try:
                    dec_raw = dec_raw.filled(np.nan)
                except (TypeError, ValueError):
                    dec_raw = dec_raw.astype(float).filled(np.nan)

            ra_num = pd.to_numeric(ra_raw, errors='coerce')
            dec_num = pd.to_numeric(dec_raw, errors='coerce')

            valid_count = np.sum(~np.isnan(ra_num) & ~np.isnan(dec_num))
            if valid_count > 0 and valid_count >= len(cat_table) * 0.3:
                ra_arr = np.array(ra_num, dtype=float)
                dec_arr = np.array(dec_num, dtype=float)
                if is_galactic:
                    coords = SkyCoord(l=ra_arr*u.deg, b=dec_arr*u.deg, frame='galactic').icrs
                    return coords.ra.deg, coords.dec.deg, 'icrs'
                else:
                    return ra_arr, dec_arr, 'icrs'

            ra_str = [str(x).strip() for x in ra_raw]
            dec_str = [str(x).strip() for x in dec_raw]

            for unit_pair in [(u.hourangle, u.deg), (u.deg, u.deg)]:
                try:
                    coords = SkyCoord(ra=ra_str, dec=dec_str, unit=unit_pair, frame='icrs')
                    return coords.ra.deg, coords.dec.deg, 'icrs'
                except Exception:
                    continue

    raise KeyError(f"Failed to locate valid sky coordinate columns in catalog. Available columns: {col_names}")

# ==========================================
# 2. MULTI-CATALOG CACHED PIPELINE
# ==========================================
def compute_multiband_stellar_nodes(
    dist_mpc, b_mag, v_mag, k_mag, w1_mag,
    default_m_l, m_sun_abs, fallback_mag=None
):
    """
    Vectorized multi-band & color-dependent stellar mass-to-light estimator:
      1) Near-IR K-band (McGaugh+14): M/L_K = 0.60, M_K_sun = 3.27
      2) Near-IR WISE W1 Vega (Lelli+16): M/L_W1 = 0.50, M_W1_sun = 3.24
      3) Optical (B - V) color (Bell+03 Table 7 - 0.15 dex Chabrier/Kroupa IMF)
      4) Legacy single-band fallback (default_m_l, m_sun_abs)
    Returns (nodes_multiband, source_codes).
    """
    dist = np.maximum(np.asarray(dist_mpc, dtype=float), 1.0e-6)
    n = dist.shape[0]

    def _to_float_arr(arr):
        if arr is None:
            return np.full(n, np.nan, dtype=float)
        return np.asarray(arr, dtype=float)

    b = _to_float_arr(b_mag)
    v = _to_float_arr(v_mag)
    k = _to_float_arr(k_mag)
    w1 = _to_float_arr(w1_mag)
    fb = b if fallback_mag is None else _to_float_arr(fallback_mag)

    dm = 5.0 * np.log10(dist * 1.0e5)

    # Priority 1: Near-IR K-band (-5 < k < 25)
    has_k = np.isfinite(k) & (k > -5.0) & (k < 25.0)
    k_safe = np.where(has_k, k, 12.0)
    m_k_abs = k_safe - dm
    mass_k = (10.0 ** (0.4 * (M_K_SUN - m_k_abs))) * M_L_K * M_SOLAR_KG

    # Priority 2: WISE W1 Vega (-5 < w1 < 25)
    has_w1 = (~has_k) & np.isfinite(w1) & (w1 > -5.0) & (w1 < 25.0)
    w1_safe = np.where(has_w1, w1, 12.0)
    m_w1_abs = w1_safe - dm
    mass_w1 = (10.0 ** (0.4 * (M_W1_SUN - m_w1_abs))) * M_L_W1 * M_SOLAR_KG

    # Priority 3: B - V optical color (Bell et al. 2003, Table 7 - 0.15 dex IMF offset)
    bv_pair_ok = np.isfinite(b) & np.isfinite(v)
    b_safe_bv = np.where(bv_pair_ok, b, 12.0)
    v_safe_bv = np.where(bv_pair_ok, v, 12.0)
    bv = b_safe_bv - v_safe_bv
    has_bv = (~has_k) & (~has_w1) & bv_pair_ok & (bv >= -0.2) & (bv <= 1.6)
    bv_clipped = np.clip(np.where(has_bv, bv, 0.5), 0.0, 1.3)
    log_ml_b = -0.942 + 1.737 * bv_clipped - 0.15
    ml_b = np.clip(10.0 ** log_ml_b, 0.10, 4.00)
    m_b_abs = b_safe_bv - dm
    mass_bv = (10.0 ** (0.4 * (M_B_SUN - m_b_abs))) * ml_b * M_SOLAR_KG

    # Priority 4: Legacy single-band fallback
    fb_ok = np.isfinite(fb) & (fb > -10.0) & (fb < 35.0)
    fb_safe = np.where(fb_ok, fb, 12.0)
    m_fb_abs = fb_safe - dm
    mass_fb = (10.0 ** (0.4 * (m_sun_abs - m_fb_abs))) * default_m_l * M_SOLAR_KG

    m_star_kg = np.select([has_k, has_w1, has_bv], [mass_k, mass_w1, mass_bv], default=mass_fb)
    source_codes = np.select([has_k, has_w1, has_bv], [1, 2, 3], default=4).astype(np.int8)

    return m_star_kg / M_PROTON, source_codes


def _get_catalog_config(catalog_name):
    """Returns schema candidates, cache filenames, and VizieR identifiers for a catalog."""
    if "Cosmicflows" in catalog_name:
        return {
            "cache_v2_filename": "cf4_catalog_cache_v2.npz",
            "cache_v1_filename": "cf4_catalog_cache.npz",
            "vizier_id": "J/ApJ/944/94",
            "default_m_l": M_L_B,
            "m_sun_abs": M_B_SUN,
            "dist_candidates": ['D', 'Dist', 'd', 'DMpc', 'D_Mpc', 'cz', 'Vcmb', 'DM', 'm-M'],
            "mag_candidates": ['Bmag', 'Kmag', 'i-mag', 'w1mag', 'm_b', 'bt', 'B', 'K', 'mag', 'm_k'],
            "name_candidates": ['Name', 'Name1', 'JNAME', 'galaxy', 'Galaxy', 'pgc', 'PGC', 'ID', 'Object', 'CF4'],
        }
    if "2M++" in catalog_name:
        return {
            "cache_v2_filename": "twompp_catalog_cache_v2.npz",
            "cache_v1_filename": "twompp_catalog_cache.npz",
            "vizier_id": "J/MNRAS/416/2840",
            "default_m_l": M_L_K,
            "m_sun_abs": M_K_SUN,
            "dist_candidates": ['Dist', 'cz', 'Vrec', 'V3k', 'Vcmb', 'cz2m', 'd', 'D', 'DM'],
            "mag_candidates": ['Ksmag', 'K20mag', 'K20', 'Kmag', 'Ks', 'K', 'mag', 'Ktot', 'm_k'],
            "name_candidates": ['Name', 'Name1', '2M++', '2MASX', 'ID', 'PGC', 'Object', 'galaxy'],
        }
    return {
        "cache_v2_filename": "ung_catalog_cache_v2.npz",
        "cache_v1_filename": "ung_catalog_cache.npz",
        "vizier_id": "J/AJ/145/101",
        "default_m_l": M_L_B,
        "m_sun_abs": M_B_SUN,
        "dist_candidates": ['Dist', 'd', 'D', 'cz', 'DM', 'm-M', 'DistMpc'],
        "mag_candidates": ['m_b', 'Bmag', 'B-mag', 'BMag', 'Kmag', 'K-mag', 'Ks', 'B', 'mag', 'bt'],
        "name_candidates": ['Name', 'Name1', 'JNAME', 'galaxy', 'Galaxy', 'PGC', 'UGCA', 'NGC', 'M', 'ID', 'Object'],
    }


def _filter_catalog_slice(gal_positions_all, nodes_all, X_all, Y_all, Z_all, gal_names_all, d_max_mpc):
    """Applies radial bounds (0.1, d_max_mpc] and builds a C++ cKDTree."""
    dists = np.linalg.norm(gal_positions_all, axis=1)
    mask = (dists > 0.1) & (dists <= d_max_mpc)
    gal_positions = gal_positions_all[mask]
    nodes_per_galaxy = nodes_all[mask]
    X, Y, Z = X_all[mask], Y_all[mask], Z_all[mask]
    gal_names = gal_names_all[mask]
    tree = cKDTree(gal_positions)
    return tree, gal_positions, nodes_per_galaxy, X, Y, Z, gal_names


@st.cache_resource(show_spinner=False)
def load_density_tree(catalog_name="UNG 2013 (Local Volume ≤ 35 Mpc)", d_max_mpc=35.0, mass_mode="multiband"):
    """
    Multi-catalog indexing engine with Streamlit resource caching and v2 multi-band persistence.
    Falls back gracefully to existing v1 cache if offline and v2 is not yet present.
    Returns: (tree, gal_positions_out, nodes_out, X_out, Y_out, Z_out, gal_names_out)
    """
    cfg = _get_catalog_config(catalog_name)
    mod_dir = os.path.dirname(__file__)
    cache_v2_path = os.path.abspath(os.path.join(mod_dir, cfg["cache_v2_filename"]))
    cache_v1_path = os.path.abspath(os.path.join(mod_dir, cfg["cache_v1_filename"]))
    vizier_id = cfg["vizier_id"]
    default_m_l = cfg["default_m_l"]
    m_sun_abs = cfg["m_sun_abs"]
    dist_candidates = cfg["dist_candidates"]
    mag_candidates = cfg["mag_candidates"]
    name_candidates = cfg["name_candidates"]

    # 1. Load from v2 Multi-Band Disk Cache if Exists
    if os.path.exists(cache_v2_path):
        try:
            data = np.load(cache_v2_path)
            gal_positions_all = data['gal_positions']
            if mass_mode == "multiband" and 'nodes_multiband' in data:
                nodes_all = data['nodes_multiband']
            elif mass_mode != "multiband" and 'nodes_legacy' in data:
                nodes_all = data['nodes_legacy']
            else:
                nodes_all = data['nodes_per_galaxy']
            X_all, Y_all, Z_all = data['X'], data['Y'], data['Z']
            if 'gal_names' in data:
                gal_names_all = data['gal_names']
            else:
                gal_names_all = np.array([f"Galaxy #{i}" for i in range(len(gal_positions_all))], dtype=str)
            return _filter_catalog_slice(
                gal_positions_all, nodes_all, X_all, Y_all, Z_all, gal_names_all, d_max_mpc
            )
        except Exception:
            pass

    # 2. Fetch Catalog via VizieR (with columns=["*"] and automatic mirror failover)
    vizier_mirrors = [
        "vizier.cfa.harvard.edu",
        "vizier.idia.ac.za",
        "vizier.nao.ac.jp",
        "vizier.china-vo.org",
        "vizier.cds.unistra.fr",
    ]

    cat_table = None
    for cat_id in [vizier_id, "J/AJ/127/2031"]:
        for mirror in vizier_mirrors:
            try:
                v = Vizier(columns=["*"], row_limit=-1, vizier_server=mirror, timeout=45)
                cats = v.get_catalogs(cat_id)
                if len(cats) > 0:
                    cat_table = cats[0]
                    break
            except Exception:
                continue
        if cat_table is not None:
            break

    # Graceful offline fallback to legacy v1 .npz cache if VizieR cannot be reached
    if cat_table is None:
        if os.path.exists(cache_v1_path):
            try:
                data = np.load(cache_v1_path)
                gal_positions_all = data['gal_positions']
                nodes_all = data['nodes_per_galaxy']
                X_all, Y_all, Z_all = data['X'], data['Y'], data['Z']
                gal_names_all = (
                    data['gal_names']
                    if 'gal_names' in data
                    else np.array([f"Galaxy #{i}" for i in range(len(gal_positions_all))], dtype=str)
                )
                return _filter_catalog_slice(
                    gal_positions_all, nodes_all, X_all, Y_all, Z_all, gal_names_all, d_max_mpc
                )
            except Exception:
                pass
        raise RuntimeError(
            f"Unable to retrieve catalog '{vizier_id}' from any VizieR mirror ({', '.join(vizier_mirrors)}) "
            f"and no local cache exists at '{cache_v1_path}'."
        )

    # Parse Coordinates, Distance, Legacy Magnitude, and Names
    ra_deg, dec_deg, coord_frame = parse_sky_coordinates(cat_table)
    dist_vals = extract_numeric_column(cat_table, dist_candidates)

    if dist_vals is not None and not np.all(np.isnan(dist_vals)):
        if np.nanmedian(dist_vals) > 500.0:
            dist_vals = dist_vals / H_GLOBAL
    else:
        dm_vals = extract_numeric_column(cat_table, ['DM', 'm-M', 'mod'])
        if dm_vals is not None and not np.all(np.isnan(dm_vals)):
            dist_vals = 10.0 ** ((dm_vals - 25.0) / 5.0)
        else:
            dist_vals = np.ones(len(ra_deg), dtype=float) * 10.0

    mag_vals = extract_numeric_column(cat_table, mag_candidates)
    if mag_vals is None or np.all(np.isnan(mag_vals)):
        mag_vals = np.full(len(ra_deg), 12.0, dtype=float)

    names_raw = extract_string_column(cat_table, name_candidates)
    if names_raw is None:
        names_raw = np.array([f"Galaxy #{i}" for i in range(len(ra_deg))], dtype=str)

    # Extract Multi-Band Photometry Independently
    def _extract_or_nan(candidates):
        arr = extract_numeric_column(cat_table, candidates)
        return np.full(len(ra_deg), np.nan, dtype=float) if arr is None else arr

    b_vals = _extract_or_nan(['Bmag', 'm_b', 'B-mag', 'BMag', 'bt', 'B'])
    v_vals = _extract_or_nan(['Vmag', 'V-mag', 'VMag', 'V'])
    k_vals = _extract_or_nan(['Ksmag', 'K20mag', 'K20', 'Kmag', 'K-mag', 'Ks', 'Ktot', 'm_k', 'K'])
    w1_vals = _extract_or_nan(['w1mag', 'W1mag', 'W1', 'w1'])

    # Best-effort PGC cross-match for Cosmicflows-4 WISE W1 photometry (Kourkchi+2020 J/ApJ/902/145)
    if "Cosmicflows" in catalog_name and np.all(np.isnan(w1_vals)):
        pgc_cf4 = extract_numeric_column(cat_table, ['PGC'])
        if pgc_cf4 is not None and not np.all(np.isnan(pgc_cf4)):
            for mirror in vizier_mirrors:
                try:
                    v_w1 = Vizier(columns=["PGC", "W1mag"], row_limit=-1, vizier_server=mirror, timeout=45)
                    w1_cats = v_w1.get_catalogs("J/ApJ/902/145")
                    if len(w1_cats) > 0:
                        t_w1 = w1_cats[0]
                        w1_pgc = extract_numeric_column(t_w1, ['PGC'])
                        w1_ab = extract_numeric_column(t_w1, ['W1mag', 'w1mag'])
                        if w1_pgc is not None and w1_ab is not None:
                            ok_w1 = np.isfinite(w1_pgc) & np.isfinite(w1_ab) & (w1_ab > -50.0) & (w1_ab < 30.0)
                            w1_vega_map = dict(
                                zip(w1_pgc[ok_w1].astype(np.int64), w1_ab[ok_w1] - WISE_W1_AB_TO_VEGA)
                            )
                            for idx_row in np.flatnonzero(np.isfinite(pgc_cf4)):
                                val_vega = w1_vega_map.get(int(pgc_cf4[idx_row]))
                                if val_vega is not None:
                                    w1_vals[idx_row] = val_vega
                        break
                except Exception:
                    continue

    # Filter Valid Mask
    has_any_mag = np.isfinite(mag_vals) | np.isfinite(k_vals) | np.isfinite(w1_vals) | np.isfinite(b_vals)
    valid_mask = (
        np.isfinite(ra_deg) & np.isfinite(dec_deg)
        & np.isfinite(dist_vals) & has_any_mag
        & (dist_vals > 0.1) & (dist_vals <= 300.0)
    )

    coords = SkyCoord(
        ra=ra_deg[valid_mask] * u.deg,
        dec=dec_deg[valid_mask] * u.deg,
        distance=dist_vals[valid_mask] * u.Mpc,
        frame=coord_frame,
    )
    mag_leg = np.where(np.isfinite(mag_vals[valid_mask]), mag_vals[valid_mask], 12.0)
    dist = dist_vals[valid_mask]
    b_sub = b_vals[valid_mask]
    v_sub = v_vals[valid_mask]
    k_sub = k_vals[valid_mask]
    w1_sub = w1_vals[valid_mask]
    names_clean = np.array([
        str(n).strip() if str(n).strip() not in ['', 'nan', 'None', '--', 'null', 'masked'] else f"Galaxy #{i}"
        for i, n in enumerate(names_raw[valid_mask])
    ], dtype=str)

    # Legacy Single-Band & Multi-Band Node Counts
    M_abs_leg = mag_leg - 5.0 * np.log10(dist * 1.0e5)
    L_lum_leg = 10.0 ** (0.4 * (m_sun_abs - M_abs_leg))
    nodes_legacy = (L_lum_leg * default_m_l * M_SOLAR_KG) / M_PROTON

    nodes_multiband, mass_source = compute_multiband_stellar_nodes(
        dist, b_sub, v_sub, k_sub, w1_sub, default_m_l, m_sun_abs, fallback_mag=mag_leg
    )

    cart = coords.cartesian
    X, Y, Z = cart.x.value, cart.y.value, cart.z.value
    gal_positions = np.vstack([X, Y, Z]).T

    # Save v2 Multi-Band Cache File
    try:
        np.savez_compressed(
            cache_v2_path,
            gal_positions=gal_positions,
            nodes_per_galaxy=nodes_multiband,
            nodes_legacy=nodes_legacy,
            nodes_multiband=nodes_multiband,
            mass_source=mass_source,
            b_mag=b_sub,
            v_mag=v_sub,
            k_mag=k_sub,
            w1_mag=w1_sub,
            X=X, Y=Y, Z=Z,
            gal_names=names_clean,
        )
    except Exception:
        pass

    nodes_active = nodes_multiband if mass_mode == "multiband" else nodes_legacy
    return _filter_catalog_slice(
        gal_positions, nodes_active, X, Y, Z, names_clean, d_max_mpc
    )


def get_catalog_photometry_stats(
    catalog_name="UNG 2013 (Local Volume ≤ 35 Mpc)",
    d_max_mpc=35.0,
    mass_mode="multiband",
    r_local_mpc=10.0
):
    """
    Returns photometric provenance breakdown (K-band, WISE W1, B-V Color, Fallback)
    plus the enclosed r_local_mpc (default 10 Mpc) sphere rho_N and H_local metrics.
    """
    _, gal_pos, nodes_active, _, _, _, _ = load_density_tree(
        catalog_name=catalog_name, d_max_mpc=d_max_mpc, mass_mode=mass_mode
    )
    n_tot = len(gal_pos)
    cfg = _get_catalog_config(catalog_name)
    mod_dir = os.path.dirname(__file__)
    cache_v2_path = os.path.abspath(os.path.join(mod_dir, cfg["cache_v2_filename"]))

    tier_keys = [(1, "K-band"), (2, "WISE W1"), (3, "B-V Color"), (4, "Fallback")]
    counts = {lbl: 0 for _, lbl in tier_keys}
    cache_version = "v1-fallback"

    if os.path.exists(cache_v2_path):
        cache_version = "v2"
        if mass_mode == "multiband":
            try:
                data = np.load(cache_v2_path)
                if 'mass_source' in data:
                    dists_all = np.linalg.norm(data['gal_positions'], axis=1)
                    mask = (dists_all > 0.1) & (dists_all <= d_max_mpc)
                    src_sub = data['mass_source'][mask]
                    for code, lbl in tier_keys:
                        counts[lbl] = int(np.sum(src_sub == code))
                else:
                    counts["Fallback"] = n_tot
            except Exception:
                counts["Fallback"] = n_tot
        else:
            counts["Fallback"] = n_tot
    else:
        counts["Fallback"] = n_tot

    percentages = {
        lbl: (100.0 * counts[lbl] / n_tot) if n_tot > 0 else 0.0
        for _, lbl in tier_keys
    }

    if cache_version == "v1-fallback":
        summary_str = "Photometry (v1 cache fallback): Single-Band / Constant M/L"
    elif mass_mode != "multiband":
        summary_str = "Photometry (v2): Legacy Constant M/L (Paper Baseline)"
    else:
        parts = []
        for _, lbl in tier_keys:
            cnt = counts[lbl]
            if cnt > 0:
                pct = percentages[lbl]
                pct_str = "<1%" if round(pct) == 0 else f"{round(pct)}%"
                parts.append(f"{pct_str} {lbl}")
        summary_str = "Photometry: " + (" | ".join(parts) if parts else "100% Fallback")

    dists_sub = np.linalg.norm(gal_pos, axis=1)
    m_loc = dists_sub <= float(r_local_mpc)
    n_loc = int(np.sum(m_loc))
    vol_m3 = (4.0 / 3.0) * np.pi * ((float(r_local_mpc) * MPC_TO_METER) ** 3)
    rho_N_10mpc = float(np.sum(nodes_active[m_loc]) / vol_m3) if n_loc > 0 else 0.0
    if rho_N_10mpc > 0.0:
        h_local_10mpc = float(H_GLOBAL + (np.sqrt((8.0 * np.pi * K_OMEGA * rho_N_10mpc) / 3.0) / UNIT_CONV))
    else:
        h_local_10mpc = float(H_GLOBAL)

    return {
        "catalog_name": catalog_name,
        "mass_mode": mass_mode,
        "cache_version": cache_version,
        "total_galaxies": n_tot,
        "counts": counts,
        "percentages": percentages,
        "summary_str": summary_str,
        "r_local_mpc": float(r_local_mpc),
        "n_local_10mpc": n_loc,
        "rho_N_10mpc": rho_N_10mpc,
        "h_local_10mpc": h_local_10mpc,
    }

# ==========================================
# 3. CORE DENSITY & EXPANSION EVALUATORS
# ==========================================
def compute_a_omega(h_local_kms_mpc):
    """Computes a_Omega from local expansion rate H_local [km/s/Mpc]."""
    h_si = h_local_kms_mpc * UNIT_CONV
    return (C_LIGHT * h_si) / (2.0 * np.pi)

def get_compiled_node_density(tree, gal_positions, nodes_per_galaxy, pos_mpc, sigma_mpc=1.8):
    """Evaluates 3D compiled node density rho_N(r) in nodes/m^3."""
    idx = tree.query_ball_point(pos_mpc, r=3.0 * sigma_mpc)
    if not idx:
        return 0.0
    dists_mpc = np.linalg.norm(gal_positions[idx] - pos_mpc, axis=1)
    kernel = np.exp(-0.5 * (dists_mpc / sigma_mpc)**2) / ((2.0 * np.pi * sigma_mpc**2)**1.5)
    rho_nodes_per_mpc3 = np.sum(nodes_per_galaxy[idx] * kernel)
    return rho_nodes_per_mpc3 / (MPC_TO_METER**3)

def get_hlocal(tree, gal_positions, nodes_per_galaxy, pos_mpc, sigma_mpc=1.8, h_global=H_GLOBAL):
    """Evaluates local unspooling rate H_local(r) = H_global + Delta_H(r)."""
    rho_N = get_compiled_node_density(tree, gal_positions, nodes_per_galaxy, pos_mpc, sigma_mpc)
    if rho_N <= 0.0:
        return h_global
    delta_H_rad = np.sqrt((8.0 * np.pi * K_OMEGA * rho_N) / 3.0)
    return h_global + (delta_H_rad / UNIT_CONV)

def get_a_omega_local(tree, gal_positions, nodes_per_galaxy, pos_mpc, sigma_mpc=1.8, h_global=H_GLOBAL):
    """Evaluates dynamic baseline acceleration a_Omega(r) in m/s^2."""
    h_loc = get_hlocal(tree, gal_positions, nodes_per_galaxy, pos_mpc, sigma_mpc=sigma_mpc, h_global=h_global)
    return compute_a_omega(h_loc)

def compute_kernel_volume_stats(sigma_mpc):
    """
    Calculates equivalent Gaussian kernel radii and volume metrics in Mpc and Mpc^3.
    - 1*sigma core sphere volume: V_1sigma = (4/3) * pi * sigma^3
    - 3*sigma enclosed radius: R_3sigma = 3 * sigma (encloses 99.7% Gaussian probability)
    - 3*sigma enclosed volume: V_3sigma = (4/3) * pi * (3*sigma)^3
    """
    r_3sigma = 3.0 * sigma_mpc
    v_1sigma = (4.0 / 3.0) * np.pi * (sigma_mpc ** 3)
    v_3sigma = (4.0 / 3.0) * np.pi * (r_3sigma ** 3)
    return {
        "r_3sigma_mpc": r_3sigma,
        "v_1sigma_mpc3": v_1sigma,
        "v_3sigma_mpc3": v_3sigma
    }
    
def get_compiled_node_density_batch(tree, gal_positions, nodes_per_galaxy, pos_array_mpc, sigma_mpc=1.8, workers=-1):
    """
    Evaluates 3D compiled node density rho_N(r) in [nodes/m^3] across an array of 3D
    positions (M, 3) using vectorized C++ cKDTree neighbor search and np.bincount.

    Note: The rho_N >= 4.0 nodes/m^3 compilation filter applies at the stellar/galactic
    phase level prior to Mpc-scale Gaussian smoothing.
    """
    pos_array = np.asarray(pos_array_mpc, dtype=float)
    if pos_array.size == 0:
        return np.array([], dtype=float)
    if pos_array.ndim == 1:
        pos_array = pos_array.reshape(1, 3)

    M = pos_array.shape[0]
    idx_lists = tree.query_ball_point(pos_array, r=3.0 * sigma_mpc, workers=workers)

    lengths = [len(x) for x in idx_lists]
    total_neighbors = sum(lengths)
    if total_neighbors == 0:
        return np.zeros(M, dtype=float)

    i_arr = np.repeat(np.arange(M, dtype=np.int64), lengths)
    j_arr = np.concatenate([np.asarray(x, dtype=np.int64) for x in idx_lists if len(x) > 0])

    diffs = gal_positions[j_arr] - pos_array[i_arr]
    dist_sq = diffs[:, 0]**2 + diffs[:, 1]**2 + diffs[:, 2]**2

    sigma_sq = sigma_mpc ** 2
    norm_factor = ((2.0 * np.pi * sigma_sq) ** 1.5) * (MPC_TO_METER ** 3)
    kernel_vals = np.exp(-0.5 * (dist_sq / sigma_sq)) / norm_factor
    weighted_nodes = nodes_per_galaxy[j_arr] * kernel_vals

    return np.bincount(i_arr, weights=weighted_nodes, minlength=M)


def get_hlocal_batch(tree, gal_positions, nodes_per_galaxy, pos_array_mpc, sigma_mpc=1.8, h_global=H_GLOBAL, workers=-1):
    """
    Evaluates local unspooling rate H_local across an array of 3D positions (M, 3)
    via get_compiled_node_density_batch using 100% vectorized NumPy indexing.
    """
    rho_N = get_compiled_node_density_batch(
        tree, gal_positions, nodes_per_galaxy, pos_array_mpc, sigma_mpc=sigma_mpc, workers=workers
    )
    if rho_N.size == 0:
        return rho_N

    coeff = np.sqrt((8.0 * np.pi * K_OMEGA) / 3.0) / UNIT_CONV
    delta_H = np.where(rho_N > 0.0, np.sqrt(np.maximum(rho_N, 0.0)) * coeff, 0.0)

    return h_global + delta_H


def evaluate_two_condition_compilation_gate(
    rho_N_m3, temp_K, radius_m, enclosed_mass_kg, mu_mol=0.60, v_bulk_ms=0.0
):
    """
    Zero-parameter Two-Condition Compilation Gate:
        Theta_compiled = Theta(rho_N - RHO_GC_THRESHOLD) * Theta(Lambda_lock - 1.0)
    where Lambda_lock = (g_river^2 + g_EM^2) / (v_th^2 + v_bulk^2) directly from the
    Absolute Bandwidth Invariant (Omega = g_river^2 + alpha_Omega^2 + v_Omega^2 <= 1.0).
    """
    rho_arr = np.asarray(rho_N_m3, dtype=float)
    T_arr   = np.maximum(np.asarray(temp_K, dtype=float), 1e-6)
    R_arr   = np.maximum(np.asarray(radius_m, dtype=float), 1e-18)
    N_enc   = np.asarray(enclosed_mass_kg, dtype=float) / M_PROTON
    mu_arr  = np.asarray(mu_mol, dtype=float)
    vb_arr  = np.asarray(v_bulk_ms, dtype=float)

    g_river_sq = (2.0 * N_enc * K_OMEGA) / (R_arr * (C_LIGHT ** 2))
    g_em_sq    = (2.0 * E_RYD_J) / (mu_arr * M_PROTON * (C_LIGHT ** 2))
    v_th_sq    = (3.0 * K_BOLTZMANN * T_arr) / (mu_arr * M_PROTON * (C_LIGHT ** 2))
    v_bulk_sq  = (vb_arr / C_LIGHT) ** 2
    v_disp_sq  = np.maximum(v_th_sq + v_bulk_sq, 1e-30)

    r_grav      = g_river_sq / v_disp_sq
    r_em        = g_em_sq / v_disp_sq
    lambda_lock = r_grav + r_em

    c1_density      = rho_arr >= RHO_GC_THRESHOLD
    c2_phase_locked = lambda_lock >= 1.0
    theta_compiled  = (c1_density & c2_phase_locked).astype(int)

    return {
        "g_river_sq": g_river_sq,
        "g_em_sq": g_em_sq,
        "v_disp_sq": v_disp_sq,
        "r_grav": r_grav,
        "r_em": r_em,
        "lambda_lock": lambda_lock,
        "c1_density": c1_density,
        "c2_phase_locked": c2_phase_locked,
        "theta_compiled": theta_compiled,
    }