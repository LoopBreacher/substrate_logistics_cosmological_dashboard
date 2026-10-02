# Substrate Galactic Lensing: Zero-Dark-Matter Optical Deflection & NFW Halo Eradication

## 🌟 Executive Summary

The `substrate_galactic_lensing` module models optical light deflection and dimensionless convergence profiles ($\kappa$) across galaxy and cluster scales without particle dark matter halos. By decoupling localized baryonic spatial address deletion ($\theta_{\text{baryon}} \propto 1/b$) from the linear metric unspooling tax ($\theta_{\text{scaffold}} \propto b$), the engine demonstrates that standard cosmology mistakes background expansion drag for quadratic Dark Matter halos ($M_{\text{DM}} \propto b^2$). The module dynamically calculates physical crossover horizons ($b_{\text{cross}}$ and $b_{1:1}$), establishing zero-parameter predictions for optical deflection parity and halo mass illusions across both preset targets and live 3D survey catalogs.

---

## 🔬 Physical Foundation & First Principles

### 1. Dual-Component Optical Deflection Angle ($\theta$)

Light traversing past a baryonic mass concentration experiences optical bending from two distinct physical mechanisms:

$$\theta_{\text{total}}(b) = \theta_{\text{baryon}}(b) + \theta_{\text{scaffold}}(b)$$

* **Localized Baryonic Spatial Address Deletion ($\theta_{\text{baryon}}$):**
Deflection caused by visible baryonic mass compiled in spatial nodes ($N_B = M_B / m_p$):

$$\theta_{\text{baryon}}(b) = \frac{4 N_B(b) \mathcal{K}_\Omega}{c^2 b} = \frac{4 G M_B(b)}{c^2 b}$$


* **Metric Scaffolding Unspooling Tax ($\theta_{\text{scaffold}}$):**
A linear metric expansion tax incurred as photons traverse expanding background space governed by the dynamic baseline acceleration $a_\Omega = \frac{c H_{\text{local}}}{2\pi}$:

$$\theta_{\text{scaffold}}(b) = \frac{4 a_\Omega b}{c^2}$$



### 2. Extended Mass Profile (Hernquist Projected Distribution)

To eliminate central point-mass mathematical singularities ($b \to 0$) in cluster cores, enclosed baryonic mass is evaluated via a projected Hernquist profile with scale radius $a_{\text{scale}}$:

$$M_{\text{enc}}(b) = M_B \left[ \frac{b^2}{(b + a_{\text{scale}})^2} \right]$$

* **Core Behavior ($b \ll a_{\text{scale}}$):** $M_{\text{enc}}(b) \propto b^2$, holding core deflection finite at $b = 0$.
* **Asymptotic Boundary ($b \gg a_{\text{scale}}$):** Enclosed mass converges smoothly to $100\%$ of total visible baryonic mass ($M_B$).

### 3. Inferred Dark Matter Halo Illusion ($M_{\text{DM}} \propto b^2$)

When standard general relativity fits assume optical deflection is produced entirely by localized mass rather than background metric unspooling, equating $\theta_{\text{scaffold}}$ to a mass term generates an artificial quadratic dark matter halo:

$$N_{\text{DM}}(b) = \frac{a_\Omega b^2}{\mathcal{K}_\Omega} \implies M_{\text{DM}}(b) = \frac{a_\Omega b^2}{G}$$

This reproduces the $M_{\text{DM}} \propto b^2$ quadratic mass growth characteristic of empirical Navarro-Frenk-White (NFW) cluster halo profiles.

### 4. Zero-Parameter Physical Horizons ($b_{\text{cross}}$ and $b_{1:1}$)

The module analytically derives two physical transition horizons:

* **1:1 Enclosed Mass Radius ($b_{\text{cross}}$):**
The impact parameter where optical deflection parity occurs ($\theta_{\text{baryon}} = \theta_{\text{scaffold}}$) and inferred dark matter equals enclosed baryonic mass ($M_{\text{DM,inferred}} = M_{\text{baryon,enclosed}}$):

$$b_{\text{cross}} = \sqrt{\frac{M_{\text{baryon}} \mathcal{K}_\Omega M_\odot}{a_\Omega m_p}} - a_{\text{scale}}$$


* **1:1 Total Mass Horizon ($b_{1:1}$):**
The impact parameter where inferred dark matter equals $100\%$ of the target's total visible baryonic mass ($M_{\text{DM,inferred}} = M_{\text{baryon,total}}$):

$$b_{1:1,\text{total}} = \sqrt{\frac{M_{\text{baryon}} \mathcal{K}_\Omega M_\odot}{a_\Omega m_p}}$$



---

## ⚙️ Architectural & Technical Implementation

* **Primary Driver Function:** `run_lensing_analysis(...)`
* **Core Profile Evaluator:** `calculate_lensing_profile(...)`
* **Zero-Parameter Bounds Calculator:** `_compute_zero_parameter_lensing_bounds(...)`

### Input Parameters

* `lens_target_name`: Name string for benchmark presets or active catalog targets.
* `custom_M_baryon`: Total visible baryonic mass $M_B$ in solar masses ($M_\odot$).
* `custom_scale_radius_kpc`: Hernquist scale radius $a_{\text{scale}}$ in kpc.
* `h_global`: Global vacuum expansion floor $H_{\text{global}}$ ($\text{km/s/Mpc}$).
* `use_log_scale`: Boolean flag toggling logarithmic plot scaling.
* `use_catalog_query`: Boolean flag enabling live 3D KD-Tree catalog queries.
* `target_pos_mpc`: 3D Cartesian coordinates $[X, Y, Z]$ of catalog target in Mpc.
* `sigma_mpc`: Gaussian smoothing scale $\sigma$ for density evaluation.
* `tree`, `gal_positions`, `nodes_per_galaxy`: Shared spatial catalog structures.

### Algorithmic Execution Pipeline

1. **Density Sampling:** If `use_catalog_query=True`, queries the active `cKDTree` dataset to derive localized compiled node density $\rho_N(\vec{r})$, local expansion rate $H_{\text{local}}(\vec{r})$, and dynamic baseline acceleration $a_\Omega(\vec{r})$.
2. **Horizon & Grid Calculation:** Computes zero-parameter physical bounds $b_{\text{cross}}$ and $b_{1:1}$, constructing a 300-point linear or logarithmic impact parameter grid $b \in [b_{\text{min}}, b_{\text{max}}]$.
3. **Profile Synthesis:** Evaluates $\theta_{\text{total}}$, $\theta_{\text{baryon}}$, $\theta_{\text{scaffold}}$, $M_{\text{DM,inferred}}(b)$, and dimensionless convergence $\kappa(b)$.
4. **Adaptive Graphic & Table Generation:** Renders dual-panel Matplotlib figures with dynamic mass units ($M_\odot$, $10^6 M_\odot$, $10^9 M_\odot$, $10^{12} M_\odot$) and formats output DataFrames with scientific precision for sub-arcsecond deflection values.

---

## 📊 Visual Diagnostic Output & Graphics

* **Panel 1 (Light Deflection Profile):** Plots total deflection angle $\theta_{\text{total}}(b)$ [red], local baryonic component $\theta_{\text{baryon}}(b)$ [dashed blue], and scaffolding tax $\theta_{\text{scaffold}}(b)$ [dotted purple]. A vertical green guide line highlights the crossover focus $b_{\text{cross}}$ where deflection parity ($50\%$) occurs.
* **Panel 2 (Dark Matter Halo Illusion):** Plots inferred halo mass $M_{\text{DM}}(b) \propto b^2$ [solid orange] against the true baryonic baseline $M_{\text{baryon}}$ [dashed blue]. A vertical cyan guide line marks the $1:1$ Total Mass Horizon $b_{1:1}$ where $M_{\text{DM,inferred}} = M_{\text{baryon,total}}$.

---

## 🎯 KPI Metrics & Export Deliverables

| Metric Key | Display Format / Unit | Physical Meaning |
| --- | --- | --- |
| `Target Baseline Acceleration (a_Ω)` | `m/s²` (e.g., `1.043e-10 m/s²`) | Dynamic acceleration baseline derived from local expansion $H_{\text{local}}$. |
| `Local Expansion Rate (H_local)` | `km/s/Mpc` | Local unspooling rate at target 3D spatial position. |
| `1:1 Enclosed Mass Radius (b_cross)` | `kpc` | Impact parameter where $M_{\text{DM,inferred}} = M_{\text{baryon,enclosed}}$ ($\theta_{\text{baryon}} = \theta_{\text{scaffold}}$). |
| `1:1 Total Mass Horizon (b_1:1)` | `kpc` | Impact parameter where $M_{\text{DM,inferred}} = M_{\text{baryon,total}}$. |
| `Deflection @ Crossover` | `arcsec` | Optical light bending angle evaluated at $b_{\text{cross}}$. |
| `Total Dark/Light @ b_cross` | Ratio (e.g., `0.43 : 1`) | Inferred dark mass relative to total visible baryonic mass at $b_{\text{cross}}$. |

### CSV Export Schema (`substrate_gravitational_lensing.csv`)

* `impact_parameter_kpc`: Radial distance $b$ from lens center in kpc.
* `theta_total_arcsec`: Total optical deflection angle $\theta_{\text{total}}$ in arcseconds.
* `theta_baryon_arcsec`: Baryonic deflection component $\theta_{\text{baryon}}$ in arcseconds.
* `theta_scaffold_arcsec`: Metric scaffolding tax $\theta_{\text{scaffold}}$ in arcseconds.
* `convergence_kappa`: Dimensionless angular convergence $\kappa$.
* `inferred_M_DM_solar`: Inferred dark matter halo mass illusion $M_{\text{DM}}$ in $M_\odot$.
* `dark_to_light_ratio`: Ratio of $M_{\text{DM,inferred}}(b) / M_{\text{baryon,total}}$.

---

## ⚡ Performance Benchmarks & Catalog Behavioral Notes

### Execution Benchmarks

* **Preset Lens Benchmark Mode:** `< 0.01s` (Instantaneous analytical profile calculation)
* **UNG 2013 (≤ 35 Mpc,):** `~0.02s` per spatial query
* **Cosmicflows-4 (≤ 150 Mpc,):** `~0.05s` per spatial query
* **2M++ Infrared (≤ 200 Mpc,):** `~0.08s` per spatial query

### Dynamic Horizon Response Across Mass Scales

* **Cluster Core Scale ($M_B \ge 10^{13} M_\odot$):** Crossover radius expands to $b_{\text{cross}} \sim 70\text{--}150\text{ kpc}$, reproducing arc distortion radii in rich clusters like Abell 1689.
* **Milky Way Scale ($M_B \sim 1.5 \times 10^{11} M_\odot$):** Crossover occurs at $b_{\text{cross}} \approx 6.16\text{ kpc}$ and $b_{1:1} \approx 14.16\text{ kpc}$, explaining flat galactic deflection profiles without cold dark matter halos.
* **Dwarf Galaxy Scale ($M_B \sim 10^6\text{--}10^8 M_\odot$):** Crossover contracts to sub-kiloparsec scales ($b_{\text{cross}} < 1.0\text{ kpc}$), where the scaffolding tax rapidly dominates outer light deflection profiles.

---