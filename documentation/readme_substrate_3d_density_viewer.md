# 🧊 Interactive 3D Cosmic Web Density Viewer (`substrate_3d_density_viewer`)

## 🌟 Executive Summary
The `substrate_3d_density_viewer` module provides an interactive 3D WebGL visualization engine powered by Plotly and Streamlit. It maps galaxy survey catalog geometry in 3D ICRS Cartesian coordinates $(X, Y, Z)$ [Mpc], inspects localized metric unspooling expansion rates $H_{\text{local}}(\vec{r})$, and renders volumetric intergalactic filament isosurface meshes directly from active survey catalogs without requiring Cold Dark Matter halos or empirical cosmological fitting parameters.

---

## 🔬 Physical Foundation & First Principles
Instead of modeling the universe as an expanding smooth background populated by dark matter halos, Substrate Logistics derives localized metric unspooling rates directly from the 3D compiled baryonic node density field $\rho_N(\vec{r})$:

* **Primary Unspooling Equation:**
  $$H_{\text{local}}(\vec{r}) = H_{\text{global}} + \sqrt{\frac{8\pi \mathcal{K}_\Omega \rho_N(\vec{r})}{3}}$$


* **Compiled Node Density Evaluation:**
  $$\rho_N(\vec{r}) = \sum_{i} \frac{N_i}{(2\pi \sigma^2)^{3/2}} \exp\left(-\frac{\Vert{}\vec{r}_i - \vec{r}\Vert{}^2}{2\sigma^2}\right) \quad \left[\text{nodes}/\text{m}^3\right]$$

* **Substrate Mechanics:**
  * **Baryonic Mass Node Conversion:** Each catalog galaxy contains integer proton nodes $N = M_B / m_p$ derived zero-parameter from photometric luminosity. Marker sizes scale logarithmically with visible baryonic mass $M_B$.
  * **Color-Coded Expansion Field:** Galaxy nodes and volumetric web meshes are color-coded by $H_{\text{local}}(\vec{r})$, highlighting overdense cluster corridors ($H_{\text{local}} \approx 72\text{--}74\text{ km/s/Mpc}$) versus underdense cosmic void interiors ($H_{\text{local}} \to H_{\text{global}} = 67.42\text{ km/s/Mpc}$).

---

## ⚙️ Architectural & Technical Implementation

* **Primary Function:** `run_3d_density_viewer_analysis(...)`
* **Input Parameters:**
  * `tree`: `scipy.spatial.cKDTree` index built over active survey galaxy coordinates.
  * `gal_positions`: `np.ndarray` of shape $(N, 3)$ representing ICRS Cartesian coordinates $(X, Y, Z)$ in Mpc.
  * `nodes_per_galaxy`: `np.ndarray` of integer proton node counts $N = M_B / m_p$.
  * `gal_names`: `list` or `np.ndarray` of target galaxy identifiers.
  * `sigma_mpc`: Gaussian field smoothing scale bandwidth $\sigma$ [Mpc].
  * `h_global`: Global vacuum floor expansion rate $H_{\text{global}} = 67.42\text{ km/s/Mpc}$.
  * `grid_res`: Volumetric sampling grid resolution $N \times N \times N$ for 3D isosurface mesh generation.
  * `max_galaxies`: Point display cap (e.g., 2,500) to optimize WebGL frame rates.
  * `show_isosurface`: `bool` toggle to enable/disable 3D filament web meshes.
  * `isosurface_opacity`: Transparency level ($0.05\text{--}0.60$) for inner void vs. outer shell visibility.

* **Algorithmic Pipeline:**
  1. **Point Selection & Mass Compilation:** Downsamples catalog nodes if $N_{\text{total}} > \text{max\_galaxies}$ and evaluates $M_B = (N \cdot m_p) / M_\odot$.
  2. **Field Evaluation:** Computes $H_{\text{local}}(\vec{r})$ for each displayed galaxy node using `get_hlocal(...)`.
  3. **3D Volumetric Mesh:** Builds a 3D grid across the spatial extent, evaluates $H_{\text{local}}$ at grid intersections, and generates `go.Isosurface` web traces.
  4. **Plotly Figure Assembly:** Renders `go.Scatter3d` nodes and `go.Isosurface` meshes formatted under a dark cosmological UI theme.

---

## 📊 Visual Diagnostic Output & Graphics

* **Interactive 3D WebGL Canvas:**
  * **X / Y / Z Axes:** Physical ICRS Cartesian coordinates in Mpc centered on Earth $(0, 0, 0)$.
  * **Galaxy Node Scatter (`go.Scatter3d`):** Color-coded by $H_{\text{local}}$ (Plasma color ramp) and scaled by baryonic mass $M_B$.
  * **3D Volumetric Isosurface (`go.Isosurface`):** Semi-transparent volumetric web shells (Viridis color ramp) highlighting high-density cosmic web corridors.
  * **Hover Tooltips:** Details galaxy name, 3D position vector $(X, Y, Z)$, radial distance [Mpc], compiled mass $[M_\odot]$, and $H_{\text{local}}$.

---

## 🎯 KPI Metrics & Export Deliverables

| Metric Key | Typical Value / Format | Physical Meaning |
| :--- | :--- | :--- |
| `Rendered Galaxies` | `500` to `5,000` | Total active galaxy nodes displayed on the 3D WebGL canvas. |
| `Peak Local Unspooling` | `73.80 km/s/Mpc` | Maximum $H_{\text{local}}$ recorded in dense cluster corridors. |
| `Void Floor Unspooling` | `68.04 km/s/Mpc` | Minimum $H_{\text{local}}$ in underdense void interiors. |
| `Mean Field Expansion` | `70.22 km/s/Mpc` | Volume-averaged expansion rate across the active survey domain. |

* **CSV Export Schema (`3d_density_web_galaxies.csv`):**
  * `galaxy_name`: Catalog target identifier.
  * `x_mpc`, `y_mpc`, `z_mpc`: ICRS Cartesian coordinates [Mpc].
  * `distance_mpc`: Radial distance $d = \sqrt{X^2 + Y^2 + Z^2}$ [Mpc].
  * `baryonic_mass_solar`: Photometrically compiled baryonic mass $[M_\odot]$.
  * `H_local_kmsMpc`: Local metric unspooling expansion rate $H_{\text{local}}$ $[\text{km/s/Mpc}]$.

---

## ⚡ Performance Benchmarks & Catalog Behavioral Notes

* **Cosmicflows-4 Macro Benchmark ($N=500$ galaxies | $d \le 150\text{ Mpc}$ | $\sigma = 13.30\text{ Mpc}$):**
  * **Enclosed Search Radius ($3\sigma$):** $39.9\text{ Mpc}$ (Search volume $V_{3\sigma} \approx 2.66 \times 10^5\text{ Mpc}^3$).
  * **Distance Range:** $1.39\text{ Mpc}$ to $149.90\text{ Mpc}$ (Mean: $91.03\text{ Mpc}$).
  * **Baryonic Mass Range:** $2.03 \times 10^7 M_\odot$ to $2.35 \times 10^{11} M_\odot$ (Mean: $1.03 \times 10^{11} M_\odot$).
  * **$H_{\text{local}}$ Field Range:** $68.04\text{ km/s/Mpc}$ (deep voids) to $73.80\text{ km/s/Mpc}$ (dense cluster cores), with a mean of $70.22\text{ km/s/Mpc}$.
* **Smoothing Scale ($\sigma$) Behavioral Impact:**
  * **Macro Scale ($\sigma = 13.30\text{ Mpc}$):** Evaluates $H_{\text{local}}$ over standard macro volumes where observational astrophysics assesses large-scale Hubble flow expansion.
  * **Fine Scale ($\sigma = 1.8\text{--}3.0\text{ Mpc}$):** Isolates individual cluster cores (e.g., Virgo, Coma), revealing localized peak spikes ($H_{\text{local}} > 76\text{ km/s/Mpc}$).