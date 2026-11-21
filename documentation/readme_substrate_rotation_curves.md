# Galactic Rotation Curves Engine (`substrate_rotation_curves.py`)

## 🌟 Executive Summary
The Galactic Rotation Curves module evaluates orbital rotation profiles $v(r)$ and radial acceleration profiles $g(r)$ across spiral and dwarf galaxies without requiring Cold Dark Matter (CDM) halos. Instead of fitting arbitrary halo profiles (e.g., NFW or isothermal spheres), orbital speeds are clamped at low accelerations by a Kinematic Underflow Floor ($a_\Omega = \frac{c H_{\text{local}}}{2\pi}$). The module supports preset benchmark galaxies (e.g., NGC 3198, Milky Way, DDO 154) and live 3D spatial catalog queries that calculate $a_\Omega(\mathbf{x})$ dynamically from the surrounding cosmic matter density.

---

## 🔬 Physical Foundation & First Principles

### 1. Dynamic Universal Baseline Acceleration ($a_\Omega$)
In Substrate Logistics, baseline acceleration is derived from local unspooling expansion rate $H_{\text{local}}(\mathbf{x})$:
$$a_\Omega(\mathbf{x}) = \frac{c \cdot H_{\text{local}}(\mathbf{x})}{2\pi}$$

In void regions or overdense filaments, $H_{\text{local}}(\mathbf{x})$ varies relative to the closed-form global floor ($H_{\text{global}} = 67.4275\text{ km/s/Mpc}$), dynamically adjusting $a_\Omega$.

### 2. Geometric Underflow Coupling
When localized baryonic acceleration drops near $a_\Omega$, effective acceleration couples geometrically:

$$g_{\text{eff}}(r) = \sqrt{g_{\text{baryon}}(r)^2 + g_{\text{baryon}}(r) \cdot a_\Omega}$$

Baryonic acceleration is computed directly from proton node counts $N_{\text{total}} = M_B / m_p$ and the derived Substrate anchor constant $\mathcal{K}_\Omega = 1.11592 \times 10^{-37}\text{ m}^3/(\text{node}\cdot\text{s}^2)$ (zero empirical $G$):

$$g_{\text{baryon}}(r) = \frac{N_{\text{total}} \mathcal{K}_\Omega}{r^2} \quad \left(\equiv \frac{G M_B}{r^2}\right)$$

### 3. Radial Baryonic Mass Geometry (Freeman Thin Disk vs. Spherical Enclosed Mass)
The solver supports two physical disk geometries via `disk_geometry`:

* **Freeman (1970) Thin Exponential Disk (Default):** Evaluates the exact 2D thin exponential disk radial acceleration using modified Bessel functions (`scipy.special.i0, i1, k0, k1`) at dimensionless radius $x = r / R_d$ ($y = x / 2$), safely clamping $x < 10^{-4}$ to $0$ and switching to the exact asymptotic $1 / x^2$ falloff for $y > 300$:

  $$g_{\text{baryon}}(r) = \frac{N_{\text{total}} \mathcal{K}_\Omega}{R_d^2} \cdot \frac{x}{2} \left[ I_0(y) K_0(y) - I_1(y) K_1(y) \right]$$

* **Spherical Enclosed Mass:** Models the integrated spherical exponential mass profile with scale length $R_d$:

  $$M_{\text{enc}}(r) = M_B \left[ 1 - \left(1 + \frac{r}{R_d}\right) e^{-r/R_d} \right]$$

  * **Exact Central Acceleration Limit:** 

    $$g_{\text{baryon}}(0) = \lim_{r \to 0} \frac{N_{\text{total}} \mathcal{K}_\Omega \left[1 - (1 + r/R_d) e^{-r/R_d}\right]}{r^2} = \frac{N_{\text{total}} \mathcal{K}_\Omega}{2 R_d^2}$$

* **Coordinate Origin Boundary:** $v(0) = 0.0\text{ km/s}$ in both geometries.

### 4. Asymptotic Velocity Floor (Baryonic Tully-Fisher Relation)
At large radii ($r \gg R^*$), $g_{\text{baryon}} \ll a_\Omega$, simplifying effective acceleration to $g_{\text{eff}} \to \sqrt{g_{\text{baryon}} a_\Omega} = \frac{\sqrt{N_{\text{total}} \mathcal{K}_\Omega a_\Omega}}{r}$. Substituting this into circular orbital velocity ($v = \sqrt{r \cdot g_{\text{eff}}}$) algebraically cancels $r$:

$$v_{\text{flat}} = \sqrt[4]{N_{\text{total}} \mathcal{K}_\Omega a_\Omega} \quad \left(\equiv \sqrt[4]{G M_B a_\Omega}\right)$$

This provides a zero-parameter derivation of the empirical Baryonic Tully-Fisher Relation ($M_B \propto v_{\text{flat}}^4$).

### 5. Transition Radius ($R^*$) & 2D Exponential Disk Truncation Cliff
* **Outer Transition Radius ($R^*$):** Identifies the outer crossing beyond the acceleration peak where $g_{\text{eff}}(R^*) = a_\Omega$ (where $g_{\text{baryon}}(R^*) \approx 0.618034 \cdot a_\Omega$):

  $$g_{\text{baryon}}(R^*) = \left(\frac{\sqrt{5}-1}{2}\right) a_\Omega \approx 0.618034 \cdot a_\Omega$$

  For ultra-diffuse or low-mass dwarf galaxies where $\max g_{\text{eff}}(r) < a_\Omega$ everywhere (e.g., DDO 154), the curve resides entirely below the acceleration floor and reports `Deep Underflow (< a_Ω)`.

* **2D Disk Truncation Cliff ($R_{\text{trunc}}$):** Derives the topological drag exponential disk cutoff boundary from `DISK_TRUNCATION_RATIO`:

  $$R_{\text{trunc}} = \frac{1}{2} \ln\left(\frac{1}{\mathcal{C}_\delta}\right) R_d \approx 3.8291 \, R_d$$

---

## ⚙️ Architectural & Technical Implementation

* **Primary Driver Function:** `run_rotation_curve_analysis(...)`
* **Math Profile Engine:** `calculate_rotation_profile(...)`
* **Input Parameters:**
  * `galaxy_target_name` (*str*): Name of target preset or catalog galaxy.
  * `custom_M_baryon` (*float*): Total visible baryonic mass in solar units ($M_\odot$).
  * `custom_R_d` (*float*): Exponential disk scale length in kpc.
  * `h_global` (*float*): Global vacuum floor $H_{\text{global}}$ ($67.4275\text{ km/s/Mpc}$).
  * `use_catalog_query` (*bool*): Enables 3D catalog spatial density sampling.
  * `target_pos_mpc` (*np.ndarray*): 3D cartesian coordinates $(X, Y, Z)$ in Mpc.
  * `sigma_mpc` (*float*): Gaussian kernel scale for local density evaluation.
  * `tree`, `gal_positions`, `nodes_per_galaxy`: Indexed catalog data structures.
  * `disk_geometry` (*str*): `"Freeman Thin Disk"` (default) or `"Spherical Enclosed Mass"`.

---

## 📊 Visual Diagnostic Output & Graphics

* **Panel 1 (Galactic Rotation Curve $v(r)$):** Plots Substrate rotation curve $v_{\text{rot}}(r)$ for the selected geometry alongside classical Keplerian decay $v_{\text{Kepler}}(r) \propto 1/\sqrt{r}$, marking the flat asymptotic floor $v_{\text{flat}}$, the outer transition radius $R^*$ (when crossed), and the vertical dashed gold 2D Disk Truncation Cliff ($R_{\text{trunc}} = 3.8291\,R_d$).
* **Panel 2 (Acceleration Profile $g(r)$ vs $a_\Omega$):** Logarithmic plot comparing $g_{\text{eff}}(r)$ and $g_{\text{baryon}}(r)$ against the Kinematic Underflow Floor $a_\Omega$.

---

## 🎯 KPI Metrics & Export Deliverables

| Metric Key | Format / Example | Physical Meaning |
| :--- | :--- | :--- |
| **True Baryonic Mass** | `32.0 × 10⁹ M☉` | Total visible stellar + gas mass (no DM halos). |
| **Disk Scale Length ($R_d$)** | `2.72 kpc` | Exponential surface density scale length. |
| **Kinematic Underflow Floor ($a_\Omega$)** | `1.043e-10 m/s²` | Dynamic minimum acceleration threshold. |
| **Transition Radius ($R^*$)** | `8.77 kpc` / `Deep Underflow (< a_Ω)` | Outer crossing where $`g_{\text{eff}} = a_\Omega`$ (if $`\max g_{\text{eff}} \ge a_\Omega`$). |
| **Flat Asymptotic Velocity** | `145.0 km/s` | Outer plateau speed $`v_{\text{flat}} = \sqrt[4]{N_{\text{total}}\mathcal{K}_\Omega a_\Omega}`$. |
| **Velocity @ Focus ($r_{\text{focus}}$)** | `153.4 km/s (@ 20.0 kpc)` | Substrate orbital velocity at the preset benchmark radius. |
| **Disk Truncation (`3.829 R_d`)** | `10.42 kpc` | Topological drag 2D exponential disk cutoff radius. |

* **CSV Export Schema (`substrate_galactic_rotation_curve.csv`):**
  * `radius_kpc`: Radial distance $r$ from core in kiloparsecs.
  * `v_substrate_kms`: Substrate orbital velocity in km/s.
  * `v_kepler_kms`: Classical Keplerian velocity in km/s.
  * `g_eff_ms2`: Effective acceleration in $\text{m/s}^2$ (full floating-point precision).
  * `g_baryon_ms2`: Pure baryonic acceleration in $\text{m/s}^2$.
  * `a_omega_floor_ms2`: Dynamic underflow acceleration floor $a_\Omega$ in $\text{m/s}^2$.

---

## ⚡ Performance Benchmarks & Catalog Behavioral Notes

* **Preset Benchmarks ($H_{\text{global}} = 67.4275\text{ km/s/Mpc}$, $a_\Omega = 1.043 \times 10^{-10}\text{ m/s}^2$):**
  * **NGC 3198 ($M_B = 3.2 \times 10^{10} M_\odot$, $R_d = 2.72\text{ kpc}$):** `v_flat = 145.0 km/s`, $R_{\text{trunc}} = 10.42\text{ kpc}$; **Freeman** $R^* = 8.77\text{ kpc}$ (`153.4 km/s` at $20\text{ kpc}$) vs. **Spherical** $R^* = 7.15\text{ kpc}$ (`148.6 km/s` at $20\text{ kpc}$).
  * **Milky Way ($M_B = 8.8 \times 10^{10} M_\odot$, $R_d = 2.30\text{ kpc}$):** `v_flat = 186.8 km/s`, $R_{\text{trunc}} = 8.81\text{ kpc}$; **Freeman** $R^* = 14.79\text{ kpc}$ (`256.1 km/s` at $8.1\text{ kpc}$) vs. **Spherical** $R^* = 13.65\text{ kpc}$ (`227.7 km/s` at $8.1\text{ kpc}$).
  * **DDO 154 ($M_B = 3.8 \times 10^8 M_\odot$, $R_d = 1.20\text{ kpc}$):** `v_flat = 47.9 km/s`, $R_{\text{trunc}} = 4.59\text{ kpc}$; **Deep Underflow (`< a_Ω`)** across all radii in both geometries (`50.4 km/s` Freeman vs `47.1 km/s` Spherical at $5.0\text{ kpc}$).
* **Execution Speed (Full Matplotlib Figure Render):** `~194 ms` headlessly across 500 radial points.
* **Smoothing Scale ($\sigma$) Sensitivity (Live Catalog Query Mode):**
  * **Fine ($\sigma = 0.5\text{--}1.0 \text{ Mpc}$):** Captures compact galaxy group overdensities, boosting $H_{\text{local}}$ and elevating $a_\Omega$.
  * **Standard ($\sigma = 1.8\text{--}2.5 \text{ Mpc}$):** Recommended default; preserves filament continuity.
  * **Coarse ($\sigma \ge 5.0 \text{ Mpc}$):** Smooths over cosmic voids, relaxing $a_\Omega$ to the global floor $a_{\Omega,\text{global}} \approx 1.0426 \times 10^{-10} \text{ m/s}^2$.