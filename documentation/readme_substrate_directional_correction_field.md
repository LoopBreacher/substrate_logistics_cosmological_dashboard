# Sky Anisotropy Maps: Full-Sky Expansion & Distance Residual Forecasts (`substrate_directional_correction_field.py`)

## 🌟 Executive Summary

The **Sky Anisotropy Maps** module ray-traces local cosmic unspooling rates across full-sky HEALPix pixels to map directional variations in the Hubble parameter ($`H_{\text{obs}}(\hat{n})`$) and predict supernova distance modulus residuals ($`\Delta\mu(\hat{n})`$) for Euclid and Rubin Observatory (LSST). By rotating Galactic HEALPix sightlines into the ICRS catalog frame (`R_GAL_TO_ICRS`) and batch-evaluating compiled proton node densities along each line of sight, the module replaces single-scalar isotropic FLRW assumptions with zero-parameter directional maps driven directly by 3D survey geometry.

---

## 🔬 Physical Foundation & First Principles

### 1. Path-Integrated Observed Expansion Field ($`H_{\text{obs}}(\hat{n})`$)

Along each full-sky Galactic unit vector $\hat{n}_G$, rays are transformed into the catalog's ICRS reference frame via the exact J2000 orthogonal rotation matrix (`R_GAL_TO_ICRS`):

$$\hat{n}_{\text{ICRS}} = \mathbf{R}_{\text{Gal}\to\text{ICRS}} \, \hat{n}_G$$

Integrating the localized unspooling rate $H_{\text{local}}(\vec{r})$ from $s = 0.1\text{ Mpc}$ out to integration depth $D_{\text{max}}$ gives the directional expansion rate:

$$H_{\text{obs}}(\hat{n}) = \frac{1}{D_{\text{max}}} \int_{0.1}^{D_{\text{max}}} H_{\text{local}}(s \, \hat{n}_{\text{ICRS}}) \, ds$$

where:

$$H_{\text{local}}(\vec{r}) = H_{\text{global}} + \sqrt{\frac{8\pi \mathcal{K}_\Omega \rho_N(\vec{r})}{3}}$$

with closed-form vacuum floor $H_{\text{global}} = 67.4275\text{ km/s/Mpc}$ and anchor constant $\mathcal{K}_\Omega = 1.11592 \times 10^{-37}\text{ m}^3/(\text{node}\cdot\text{s}^2)$.

### 2. Directional Multiplier & SNe Ia Residual Anomaly

* **Directional Expansion Multiplier ($C$):**

  $$C(\hat{n}) = \frac{H_{\text{obs}}(\hat{n})}{H_{\text{global}}}$$

* **Euclid / Rubin (LSST) Distance Modulus Anomaly ($\Delta\mu$):**

  $$\Delta\mu(\hat{n}) = -5 \log_{10}\left( C(\hat{n}) \right)$$

Sightlines piercing dense supercluster corridors yield $C(\hat{n}) > 1$ and $\Delta\mu < 0$ (supernovae appear systematically brighter than an isotropic baseline), whereas deep void sightlines converge toward $C(\hat{n}) \to 1.0$ and $\Delta\mu \to 0.0\text{ mag}$.

---

## ⚙️ Architectural & Technical Implementation

* **Primary Driver Function:** `run_anisotropic_sky_forecast(...)`
* **Cached Ray-Tracing Engine:** `_compute_sky_forecast_cached(...)`
* **Core Dependencies:** `healpy` (or `modules/healpy_compat.py` via `astropy_healpix`), `modules.data_loader` (`get_hlocal_batch`, `galactic_to_icrs_vecs`, `H_GLOBAL`)

### Input Parameters

* `tree`: `scipy.spatial.cKDTree` spatial index built over active ICRS galaxy coordinates.
* `gal_positions`: `np.ndarray` of shape $(N, 3)$ in Mpc (ICRS Cartesian).
* `nodes_per_galaxy`: `np.ndarray` of proton node counts $N = M / m_p$.
* `nside`: HEALPix resolution parameter ($8, 16, 32$; default `32`, corresponding to $768$, $3,072$, or $12,288$ sky pixels).
* `d_max_mpc`: Maximum ray-tracing integration depth in Mpc.
* `sigma_smooth_mpc`: Gaussian kernel bandwidth $\sigma$ in Mpc.
* `h_global`: Global vacuum floor expansion rate ($67.4275\text{ km/s/Mpc}$).

### Cache Safety & Coordinate Alignment

1. **J2000 Galactic-to-ICRS Rotation:** Pixel directions generated on the Galactic HEALPix sphere (`coord=["G"]`) are rotated into the ICRS catalog frame via `galactic_to_icrs_vecs(pix_vecs)` prior to spatial tree queries.
2. **Catalog-Aware Caching (`catalog_sig`):** `_compute_sky_forecast_cached` keys Streamlit's `@st.cache_data` store with `catalog_sig = (len(gal_positions), round(float(np.sum(nodes_per_galaxy)), 6))`, guaranteeing immediate cache invalidation whenever the user switches survey catalogs or volume depths.

---

## 📊 Visual Diagnostic Output & Graphics

The module outputs a 3-panel vertical stack of full-sky Mollweide projections rendered in Galactic coordinates (`coord=["G"]`, astronomical orientation `flip="astro"`):

1. **Panel 1 (`plasma`):** Observed directional expansion field $H_{\text{obs}}(\theta, \phi)$ $[\text{km/s/Mpc}]$.
2. **Panel 2 (`magma`):** Directional correction multiplier $C(\theta, \phi) = H_{\text{obs}} / H_{\text{global}}$.
3. **Panel 3 (`viridis_r`):** Predicted Euclid & Rubin (LSST) Type Ia supernova distance modulus anomaly $\Delta\mu(\theta, \phi)$ $[\text{mag}]$.

---

## 🎯 KPI Metrics & Export Deliverables

| Metric Key | Unit / Format | Physical Meaning |
| :--- | :--- | :--- |
| `Mean Sky Expansion` | `km/s/Mpc` | Full-sky average of path-integrated expansion $H_{\text{obs}}(\hat{n})$. |
| `Correction Range` | `x` | Minimum to maximum directional multiplier $C(\hat{n})$ across the sky. |
| `LSST/Euclid Anomaly Range` | `mag` | Predicted SNe Ia distance modulus residual span $\Delta\mu(\hat{n})$. |

### CSV Export Schema (`healpix_sky_anisotropy_lookup.csv`)

* `pixel_id`: HEALPix ring-ordered pixel index.
* `glon_deg`: Galactic longitude $l$ $[0^\circ \dots 360^\circ]$.
* `glat_deg`: Galactic latitude $b$ $[-90^\circ \dots +90^\circ]$.
* `H_obs_kmsMpc`: Path-integrated directional expansion rate $[\text{km/s/Mpc}]$.
* `Correction_Factor_C`: Directional expansion factor $C(\hat{n})$.
* `Delta_mu_mag`: Predicted distance modulus anomaly $\Delta\mu(\hat{n})$ $[\text{mag}]$.

---

## ⚡ Performance Benchmarks & Catalog Behavioral Notes

* **Cosmicflows-4 Benchmark ($N = 17,509$, $D_{\text{max}} = 150\text{ Mpc}$, $\sigma = 2.0\text{ Mpc}$, $N_{\text{side}} = 16$):**
  * **Mean Sky Expansion:** `68.38 km/s/Mpc`
  * **Directional Multiplier Range ($C$):** `1.0003x` to `1.1194x`
  * **Distance Modulus Anomaly Range ($\Delta\mu$):** `-0.2449 mag` to `-0.0007 mag`
* **Headless Execution Time ($N_{\text{side}} = 16$, $d_{\text{max}} = 15\text{ Mpc}$, full Matplotlib render):** `~1.65 s` via batched `cKDTree` evaluation (`get_hlocal_batch`).
