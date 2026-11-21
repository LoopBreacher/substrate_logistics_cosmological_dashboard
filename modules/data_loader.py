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
M_B_SUN          = 5.44                # Solar absolute B-band magnitude (UNG 2013 / Cosmicflows)
M_K_SUN          = 3.27                # Solar absolute K-band magnitude (2M++ Infrared)
M_L_B            = 0.44                # B-band mass-to-light ratio (Bell & de Jong 2001)
M_L_K            = 0.60                # K-band mass-to-light ratio (Near-Infrared)

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
@st.cache_resource(show_spinner=False)
def load_density_tree(catalog_name="UNG 2013 (Local Volume ≤ 35 Mpc)", d_max_mpc=35.0):
    """
    Multi-catalog indexing engine with Streamlit resource caching and disk persistence.
    Returns: (tree, gal_positions_out, nodes_out, X_out, Y_out, Z_out, gal_names_out)
    """
    if "Cosmicflows" in catalog_name:
        cache_filename = "cf4_catalog_cache.npz"
        vizier_id = "J/ApJ/944/94"
        default_m_l = M_L_B
        m_sun_abs = M_B_SUN
        dist_candidates = ['D', 'Dist', 'd', 'DMpc', 'D_Mpc', 'cz', 'Vcmb', 'DM', 'm-M']
        mag_candidates = ['Bmag', 'Kmag', 'i-mag', 'w1mag', 'm_b', 'bt', 'B', 'K', 'mag', 'm_k']
        name_candidates = ['Name', 'Name1', 'JNAME', 'galaxy', 'Galaxy', 'pgc', 'PGC', 'ID', 'Object', 'CF4']
    elif "2M++" in catalog_name:
        cache_filename = "twompp_catalog_cache.npz"
        vizier_id = "J/MNRAS/416/2840"
        default_m_l = M_L_K
        m_sun_abs = M_K_SUN
        dist_candidates = ['Dist', 'cz', 'Vrec', 'V3k', 'Vcmb', 'cz2m', 'd', 'D', 'DM']
        mag_candidates = ['K20mag', 'K20', 'Kmag', 'Ks', 'K', 'mag', 'Ktot', 'm_k']
        name_candidates = ['Name', 'Name1', '2M++', '2MASX', 'ID', 'PGC', 'Object', 'galaxy']
    else:  # Default UNG 2013
        cache_filename = "ung_catalog_cache.npz"
        vizier_id = "J/AJ/145/101"
        default_m_l = M_L_B
        m_sun_abs = M_B_SUN
        dist_candidates = ['Dist', 'd', 'D', 'cz', 'DM', 'm-M', 'DistMpc']
        mag_candidates = ['m_b', 'Bmag', 'B-mag', 'BMag', 'Kmag', 'K-mag', 'Ks', 'B', 'mag', 'bt']
        name_candidates = ['Name', 'Name1', 'JNAME', 'galaxy', 'Galaxy', 'PGC', 'UGCA', 'NGC', 'M', 'ID', 'Object']

    cache_path = os.path.abspath(os.path.join(os.path.dirname(__file__), cache_filename))

    # 1. Load from Disk Cache if Exists
    if os.path.exists(cache_path):
        try:
            data = np.load(cache_path)
            gal_positions_all = data['gal_positions']
            nodes_all = data['nodes_per_galaxy']
            X_all, Y_all, Z_all = data['X'], data['Y'], data['Z']

            if 'gal_names' in data:
                gal_names_all = data['gal_names']
            else:
                gal_names_all = np.array([f"Galaxy #{i}" for i in range(len(gal_positions_all))], dtype=str)

            dists = np.linalg.norm(gal_positions_all, axis=1)
            mask = (dists > 0.1) & (dists <= d_max_mpc)

            gal_positions = gal_positions_all[mask]
            nodes_per_galaxy = nodes_all[mask]
            X, Y, Z = X_all[mask], Y_all[mask], Z_all[mask]
            gal_names = gal_names_all[mask]

            tree = cKDTree(gal_positions)
            return tree, gal_positions, nodes_per_galaxy, X, Y, Z, gal_names
        except Exception:
            pass

    # 2. Fetch Catalog via VizieR (with automatic mirror failover for CDS outage)
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
                v = Vizier(row_limit=-1, vizier_server=mirror, timeout=30)
                cats = v.get_catalogs(cat_id)
                if len(cats) > 0:
                    cat_table = cats[0]
                    break
            except Exception:
                continue
        if cat_table is not None:
            break

    if cat_table is None:
        raise RuntimeError(
            f"Unable to retrieve catalog '{vizier_id}' from any VizieR mirror ({', '.join(vizier_mirrors)})."
        )

    # Parse Coordinates, Distance, Magnitude, and Names
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
            dist_vals = np.ones(len(ra_deg)) * 10.0

    mag_vals = extract_numeric_column(cat_table, mag_candidates)
    if mag_vals is None or np.all(np.isnan(mag_vals)):
        mag_vals = np.full(len(ra_deg), 12.0)

    names_raw = extract_string_column(cat_table, name_candidates)
    if names_raw is None:
        names_raw = np.array([f"Galaxy #{i}" for i in range(len(ra_deg))], dtype=str)

    # Filter Valid Mask
    valid_mask = (
        (~np.isnan(ra_deg)) & (~np.isnan(dec_deg)) &
        (~np.isnan(dist_vals)) & (~np.isnan(mag_vals)) &
        (dist_vals > 0.1) & (dist_vals <= 300.0)
    )

    coords = SkyCoord(
        ra=ra_deg[valid_mask]*u.deg, 
        dec=dec_deg[valid_mask]*u.deg, 
        distance=dist_vals[valid_mask]*u.Mpc, 
        frame=coord_frame
    )
    mag = mag_vals[valid_mask]
    dist = dist_vals[valid_mask]
    names_clean = np.array([
        str(n).strip() if str(n).strip() not in ['', 'nan', 'None', '--', 'null', 'masked'] else f"Galaxy #{i}"
        for i, n in enumerate(names_raw[valid_mask])
    ], dtype=str)

    # Calculate Absolute Magnitude, Luminosity, Stellar Mass, and Node Count
    M_abs = mag - 5.0 * np.log10(dist * 1e5)
    L_lum = 10.0 ** (0.4 * (m_sun_abs - M_abs))
    M_compiled_kg = L_lum * default_m_l * M_SOLAR_KG

    cart = coords.cartesian
    X, Y, Z = cart.x.value, cart.y.value, cart.z.value
    gal_positions = np.vstack([X, Y, Z]).T
    nodes_per_galaxy = M_compiled_kg / M_PROTON

    # Save Cache File
    try:
        np.savez_compressed(
            cache_path,
            gal_positions=gal_positions,
            nodes_per_galaxy=nodes_per_galaxy,
            X=X, Y=Y, Z=Z,
            gal_names=names_clean
        )
    except Exception:
        pass

    # Filter to requested d_max_mpc
    dists_filtered = np.linalg.norm(gal_positions, axis=1)
    sub_mask = dists_filtered <= d_max_mpc

    gal_positions_out = gal_positions[sub_mask]
    nodes_out = nodes_per_galaxy[sub_mask]
    X_out, Y_out, Z_out = X[sub_mask], Y[sub_mask], Z[sub_mask]
    gal_names_out = names_clean[sub_mask]
    tree = cKDTree(gal_positions_out)

    return tree, gal_positions_out, nodes_out, X_out, Y_out, Z_out, gal_names_out

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