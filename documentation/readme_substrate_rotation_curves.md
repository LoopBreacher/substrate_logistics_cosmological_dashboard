# Galactic Rotation Curves Engine (`substrate_rotation_curves.py`)

## 🌟 Executive Summary
The Galactic Rotation Curves module evaluates orbital rotation profiles $v(r)$ and radial acceleration profiles $g(r)$ across spiral, giant LSB, and dwarf galaxies without particle Cold Dark Matter (CDM) or empirical $G$. It supports two unified comparison modes:

1. **Topological Cover Pocket (Linear $\mathcal{K}_\Omega$ + $3P_\delta$ Burkert, Default):** Combines pure linear $\mathcal{K}_\Omega$ baryons (**Freeman thin exponential disk** or **Spherical Enclosed Mass**) with a momentum-conserving cored Burkert Cover Pocket anchored to the 3D Topological Drag Pressure Floor $3 P_\delta$, naturally recovering the Baryonic Tully-Fisher Acceleration $a_{\text{BTFR}} = 8.9844 \times 10^{-11}\ \text{m/s}^2$.
2. **Legacy Kinematic Underflow ($a_\Omega$ Floor):** Clamps low-acceleration regions via $g_{\text{eff}} = \sqrt{g_{\text{baryon}}^2 + g_{\text{baryon}}\,a_\Omega}$ using the global or local catalog expansion floor $a_\Omega(\mathbf{x}) = \frac{c H_{\text{local}}(\mathbf{x})}{2\pi}$.

In both modes, the engine overlays the inactive architecture as a dashed reference curve so outer Burkert roll-off beyond $6\text{--}8\,r_0$ vs. flat $a_\Omega$ plateaus can be inspected side-by-side.

---

## 🔬 Physical Foundation & First Principles

### 1. Zero-Parameter $3P_\delta$ Topological Cover Pocket (Recommended)
All Cover Pocket invariants are derived algebraically in [`modules/data_loader.py`](../modules/data_loader.py) from `K_OMEGA`, `M_PROTON`, `P_DELTA_MACRO_PA`, `OMEGA_M_VOID`, and `OMEGA_M_GLOBAL`:

- **Effective Coupling:** $G_{\text{eff}} = \mathcal{K}_\Omega / m_p = 6.671657 \times 10^{-11}\ \text{m}^3/(\text{kg}\cdot\text{s}^2)$
- **3D Drag Pressure Floor:** $P_{\delta,\text{3D}} = 3 P_{\delta,\text{macro}} = 8.501710 \times 10^{-13}\ \text{Pa}$
- **Peak Acceleration Scale:** $g_{3P_\delta} = \sqrt{8\pi G_{\text{eff}}(3P_\delta)} = 3.775634 \times 10^{-11}\ \text{m/s}^2$
- **Universal Burkert Surface Density:**
  $$\Sigma_p = \frac{g_{3P_\delta}}{f_{B,\max}\,\pi\,G_{\text{eff}}} = 169.47\ M_\odot/\text{pc}^2 \quad (f_{B,\max} = 0.508967\ \text{at}\ x_{g,\max} = 0.96341)$$
- **Mass-Locked Pocket Core Radius ($r_0$):** Locked via universal pocket-to-baryon ratio $\eta_{\text{pocket}} = \Omega_{m,\text{void}} / (\Omega_{m,\text{global}} - \Omega_{m,\text{void}}) = 4.566048$:
  $$r_0 = \sqrt{\frac{\eta_{\text{pocket}} M_B}{\pi\,\Sigma_p\,m_{\text{fac}}(x_{v,\max})}} \quad (m_{\text{fac}} = 2.792378\ \text{at}\ x_{v,\max} = 3.24461)$$
- **Linear Superposition & Circular Velocity:**
  $$g_{\text{tot}}(r) = g_{\text{baryon}}(r) + \pi G_{\text{eff}} \Sigma_p \frac{\ln(1+x^2) + 2\ln(1+x) - 2\arctan x}{x^2}, \quad x = \frac{r}{r_0}, \quad v_{\text{pocket}}(r) = \sqrt{r\,g_{\text{tot}}(r)}$$
- **Emergent Baryonic Tully-Fisher Acceleration:**
  $$a_{\text{BTFR}} = \frac{v_{B,\max}^4}{f_{B,\max}\,m_{\text{fac}}(x_{v,\max})}\,\eta_{\text{pocket}}\,g_{3P_\delta} = 8.9844 \times 10^{-11}\ \text{m/s}^2$$

### 2. Baryonic Mass Geometry (Freeman Thin Disk vs. Spherical Enclosed Mass)
* **Freeman (1970) Thin Exponential Disk (Default):** Evaluates the exact 2D thin exponential disk radial acceleration via modified Bessel functions (`scipy.special.i0, i1, k0, k1`) at $x = r / R_d$ ($y = x / 2$):
  $$g_{\text{baryon}}(r) = \frac{N_{\text{total}} \mathcal{K}_\Omega}{R_d^2} \cdot \frac{x}{2} \left[ I_0(y) K_0(y) - I_1(y) K_1(y) \right]$$
* **Spherical Enclosed Mass:** Integrates the spherical exponential mass profile $M_{\text{enc}}(r) = M_B [1 - (1 + r/R_d)e^{-r/R_d}]$.
* **2D Disk Truncation Cliff:** $R_{\text{trunc}} = \frac{1}{2}\ln(1/\mathcal{C}_\delta)\,R_d \approx 3.8291\,R_d$.

### 3. Legacy Kinematic Underflow Overlay ($a_\Omega$)
When running in Legacy mode (or as the dashed comparison overlay), radial acceleration couples via:
$$g_{\text{eff}}(r) = \sqrt{g_{\text{baryon}}(r)^2 + g_{\text{baryon}}(r)\cdot a_\Omega}, \qquad a_\Omega(\mathbf{x}) = \frac{c H_{\text{local}}(\mathbf{x})}{2\pi}$$
yielding flat asymptotic speed $v_{\text{flat}} = \sqrt[4]{N_{\text{total}}\mathcal{K}_\Omega a_\Omega}$ and outer crossing radius $R^*$ where $g_{\text{eff}}(R^*) = a_\Omega$ (or `Deep Underflow (< a_Ω)` when $\max g_{\text{eff}} < a_\Omega$).

---

## ⚙️ Architectural & Technical Implementation

* **Primary UI Driver:** `run_rotation_curve_analysis(..., gravity_model="Cover Pocket (Linear K_Ω + 3 P_δ Burkert)")`
* **Comparative Math Evaluator:** `calculate_cover_pocket_rotation_profile(...)`
* **Backward-Compatible Profile API:** `calculate_rotation_profile(M_baryon_solar, r_kpc_array, a_omega_active=A_OMEGA, R_d_kpc=None, disk_geometry="Freeman Thin Disk", gravity_model="Legacy Underflow (a_Ω)")`
* **Presets (`ROTATION_PRESETS`):** `NGC 3198 (Classic Spiral Benchmark)`, `Milky Way (Stellar Disk + Bulge + Gas)`, `DDO 154 (Gas-Rich Dwarf Galaxy)`, and `UGC 2885 (Giant Spiral Benchmark)`.

---

## 📊 Visual Diagnostic Output & Export Schema

* **Panel 1 (Rotation Curve $v(r)$):** Active primary curve (`#ef4444` solid) overlaid against the alternative mode (`#a855f7` dashed), `Cover Pocket Alone` (`#10b981` dash-dot), `Baryonic Keplerian` (`#38bdf8` dotted), locked core $r_0$, and truncation cliff $3.8291\,R_d$.
* **Panel 2 (Radial Acceleration Log Profile):** Compares $g_{\text{tot}}(r)$, alternative $g_{\text{alt}}(r)$, and $g_{\text{baryon}}(r)$ against $g_{3P_\delta} = 3.776 \times 10^{-11}\ \text{m/s}^2$, $a_{\text{BTFR}} = 8.984 \times 10^{-11}\ \text{m/s}^2$, and $a_\Omega$.
* **CSV Export Schema (`substrate_galactic_rotation_curve.csv`):** `radius_kpc`, `v_pocket_kms`, `v_pocket_only_kms`, `v_legacy_aOmega_kms`, `v_kepler_kms`, `v_substrate_kms`, `g_pocket_tot_ms2`, `g_legacy_eff_ms2`, `g_baryon_ms2`, `g_pocket_ms2`, `a_omega_floor_ms2`.

### KPI Metrics (`metrics` dictionary)
| Gravity Mode | KPI Keys Returned |
| :--- | :--- |
| **Cover Pocket (Default)** | `True Baryonic Mass`, `Locked Pocket Core (r_0)`, `Pocket Surface Density (Σ_p)`, `Emergent a_BTFR`, `Velocity @ Focus (r_focus)`, `Disk Truncation (3.829 R_d)` |
| **Legacy Underflow ($a_\Omega$)** | `True Baryonic Mass`, `Disk Scale Length (R_d)`, `Kinematic Underflow Floor (a_Ω)`, `Transition Radius (R*)`, `Flat Asymptotic Velocity`, `Velocity @ Focus (r_focus)`, `Disk Truncation (3.829 R_d)` |
