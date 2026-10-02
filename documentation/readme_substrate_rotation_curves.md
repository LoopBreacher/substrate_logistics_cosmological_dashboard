# Galactic Rotation Curves Engine (`substrate_rotation_curves.py`)

## 🌟 Executive Summary
The Galactic Rotation Curves module evaluates orbital rotation profiles $v(r)$ and radial acceleration profiles $g(r)$ across spiral and dwarf galaxies without requiring Cold Dark Matter (CDM) halos. Instead of fitting arbitrary halo profiles (e.g., NFW or isothermal spheres), orbital speeds are clamped at low accelerations by a Kinematic Underflow Floor ($a_\Omega = \frac{c H_{\text{local}}}{2\pi}$). The module supports preset benchmark galaxies (e.g., NGC 3198, Milky Way, DDO 154) and live 3D spatial catalog queries that calculate $a_\Omega(\mathbf{x})$ dynamically from the surrounding cosmic matter density.

---

## 🔬 Physical Foundation & First Principles

### 1. Dynamic Universal Baseline Acceleration ($a_\Omega$)
In Substrate Logistics, baseline acceleration is derived from local unspooling expansion rate $H_{\text{local}}(\mathbf{x})$:
$$a_\Omega(\mathbf{x}) = \frac{c \cdot H_{\text{local}}(\mathbf{x})}{2\pi}$$

In void regions or overdense filaments, $H_{\text{local}}(\mathbf{x})$ varies relative to the global floor ($H_{\text{global}} = 67.42 \text{ km/s/Mpc}$), dynamically adjusting $a_\Omega$.

### 2. Geometric Underflow Coupling
When localized baryonic acceleration drops near $a_\Omega$, effective acceleration couples geometrically:
$$g_{\text{eff}}(r) = \sqrt{g_{\text{baryon}}(r)^2 + g_{\text{baryon}}(r) \cdot a_\Omega}$$

Baryonic acceleration is computed directly from proton node counts $N_{\text{total}} = M_B / m_p$ and Substrate constant $\mathcal{K}_\Omega$:
$$g_{\text{baryon}}(r) = \frac{N_{\text{total}} \mathcal{K}_\Omega}{r^2} \quad \left(\equiv \frac{G M_B}{r^2}\right)$$

### 3. Extended Mass Profile (Exponential Disk)
To prevent core point-mass singularities ($r \to 0$), the enclosed mass profile $M_{\text{enc}}(r)$ models an exponential disk with scale length $R_d$:
$$M_{\text{enc}}(r) = M_B \left[ 1 - \left(1 + \frac{r}{R_d}\right) e^{-r/R_d} \right]$$

* **Exact Central Acceleration Limit:** $g_{\text{baryon}}(0) = \lim_{r \to 0} \frac{N_{\text{total}} \mathcal{K}_\Omega \left[1 - (1 + r/R_d) e^{-r/R_d}\right]}{r^2} = \mathbf{\frac{N_{\text{total}} \mathcal{K}_\Omega}{2 R_d^2}}$
* **Coordinate Origin Boundary:** $v(0) = 0.0 \text{ km/s}$

### 4. Asymptotic Velocity Floor (Baryonic Tully-Fisher Relation)
At large radii ($r \gg R^*$), $g_{\text{baryon}} \ll a_\Omega$, simplifying effective acceleration to $g_{\text{eff}} \to \sqrt{g_{\text{baryon}} a_\Omega} = \frac{\sqrt{N_{\text{total}} \mathcal{K}_\Omega a_\Omega}}{r}$. Substituting this into circular orbital velocity ($v = \sqrt{r \cdot g_{\text{eff}}}$) algebraically cancels $r$:
$$v_{\text{flat}} = \sqrt[4]{N_{\text{total}} \mathcal{K}_\Omega a_\Omega} \quad \left(\equiv \sqrt[4]{G M_B a_\Omega}\right)$$

This provides a zero-parameter derivation of the empirical Baryonic Tully-Fisher Relation ($M_B \propto v_{\text{flat}}^4$).

### 5. Transition Radius ($R^*$)
The transition boundary $R^*$ marks where effective acceleration crosses the underflow floor ($g_{\text{eff}}(R^*) = a_\Omega$), which occurs when baryonic acceleration decays to:
$$g_{\text{baryon}}(R^*) = \left(\frac{\sqrt{5}-1}{2}\right) a_\Omega \approx 0.618034 \cdot a_\Omega$$

---

## ⚙️ Architectural & Technical Implementation

* **Primary Driver Function:** `run_rotation_curve_analysis(...)`
* **Math Profile Engine:** `calculate_rotation_profile(...)`
* **Input Parameters:**
  * `galaxy_target_name` (*str*): Name of target preset or catalog galaxy.
  * `custom_M_baryon` (*float*): Total visible baryonic mass in solar units ($M_\odot$).
  * `custom_R_d` (*float*): Exponential disk scale length in kpc.
  * `h_global` (*float*): Global vacuum floor $H_{\text{global}}$ ($67.42 \text{ km/s/Mpc}$).
  * `use_catalog_query` (*bool*): Enables 3D catalog spatial density sampling.
  * `target_pos_mpc` (*np.ndarray*): 3D cartesian coordinates $(X, Y, Z)$ in Mpc.
  * `sigma_mpc` (*float*): Gaussian kernel scale for local density evaluation.
  * `tree`, `gal_positions`, `nodes_per_galaxy`: Indexed catalog data structures.

---

## 📊 Visual Diagnostic Output & Graphics

* **Panel 1 (Galactic Rotation Curve $v(r)$):** Plots Substrate rotation curve $v_{\text{rot}}(r)$ alongside classical Keplerian decay $v_{\text{Kepler}}(r) \propto 1/\sqrt{r}$, marking the flat asymptotic floor $v_{\text{flat}}$ and vertical transition radius $R^*$.
* **Panel 2 (Acceleration Profile $g(r)$ vs $a_\Omega$):** Logarithmic plot comparing $g_{\text{eff}}(r)$ and $g_{\text{baryon}}(r)$ against the Kinematic Underflow Floor $a_\Omega$.

---

## 🎯 KPI Metrics & Export Deliverables

| Metric Key | Format / Example | Physical Meaning |
| :--- | :--- | :--- |
| **True Baryonic Mass** | `88.0 × 10⁹ M☉` | Total visible stellar + gas mass (no DM halos). |
| **Disk Scale Length ($R_d$)** | `2.30 kpc` | Exponential surface density scale length. |
| **Kinematic Underflow Floor ($a_\Omega$)** | `1.129e-10 m/s²` | Dynamic minimum acceleration threshold. |
| **Transition Radius ($R^*$)** | `13.65 kpc` | Radius where $g_{\text{eff}} = a_\Omega$. |
| **Flat Asymptotic Velocity** | `186.8 km/s` | Outer plateau speed $v_{\text{flat}} = \sqrt[4]{N_{\text{total}}\mathcal{K}_\Omega a_\Omega}$. |

* **CSV Export Schema (`substrate_galactic_rotation_curve.csv`):**
  * `radius_kpc`: Radial distance $r$ from core in kiloparsecs.
  * `v_substrate_kms`: Substrate orbital velocity in km/s.
  * `v_kepler_kms`: Classical Keplerian velocity in km/s.
  * `g_eff_ms2`: Effective acceleration in $\text{m/s}^2$ (full floating-point precision).
  * `g_baryon_ms2`: Pure baryonic acceleration in $\text{m/s}^2$.
  * `a_omega_floor_ms2`: Dynamic underflow acceleration floor $a_\Omega$ in $\text{m/s}^2$.

---

## ⚡ Performance Benchmarks & Catalog Behavioral Notes

* **Mathematical Profile Engine (`calculate_rotation_profile`):** `~0.33 ms` (500 grid points).
* **Preset Rendering Execution Time (includes Matplotlib rendering):**
  * **NGC 3198:** `~1,200 ms`
  * **Milky Way:** `~440 ms`
  * **DDO 154:** `~530 ms`
* **Live Catalog Query Execution Time (3D `cKDTree` density evaluation):**
  * **UNG 2013 (≤ 35 Mpc, 3k galaxies):** `~480 ms`
  * **Cosmicflows-4 (≤ 150 Mpc, 35k galaxies):** `~520 ms`
  * **2M++ Infrared (≤ 200 Mpc, 70k galaxies):** `~440 ms`
* **Smoothing Scale ($\sigma$) Sensitivity:**
  * **Fine ($\sigma = 0.5\text{--}1.0 \text{ Mpc}$):** Captures compact galaxy group overdensities, boosting $H_{\text{local}}$ and elevating $a_\Omega$.
  * **Standard ($\sigma = 1.8\text{--}2.5 \text{ Mpc}$):** Recommended default; preserves filament continuity.
  * **Coarse ($\sigma \ge 5.0 \text{ Mpc}$):** Smooths over cosmic voids, relaxing $a_\Omega$ to the global floor $a_{\Omega,\text{global}} \approx 1.043 \times 10^{-10} \text{ m/s}^2$.