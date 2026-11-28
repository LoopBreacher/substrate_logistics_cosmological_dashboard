# Substrate Galactic Lensing: Unified Two-Ceiling Cover Pocket & Strong/Weak Lensing Engine

## 🌟 Executive Summary

The `substrate_galactic_lensing` module models 2D optical light deflection $\theta(b)$, surface convergence $\kappa(b)$, enclosed mean convergence $\bar{\kappa}(<b)$, strong-lensing Einstein radii $R_E$, and projected cylinder masses $M_{2D}(<250\ \text{kpc})$ without particle dark matter. It provides two selectable physical architectures:

1. **Topological Cover Pocket Lensing (Linear $\mathcal{K}_\Omega$ + Burkert, Default):** Evaluates exact closed-form 2D projected Burkert surface density $\Sigma_{2D}(X)$ across both the **Galaxy Pressure Floor ($3 P_\delta \implies \Sigma_p = 169.47\ M_\odot/\mathrm{pc}^2$)** and the **Isotropic Cluster Ceiling ($\Sigma_0 / 3 \implies \Sigma_p = 2,599.82\ M_\odot/\mathrm{pc}^2$)** linked by the exact topological bridge $1 / (9\mathcal{C}_\delta) = 235.34$.
2. **Legacy Linear Scaffold ($a_\Omega$):** Preserves the $1/b$ Hernquist baryonic term plus $\theta_{\text{scaffold}} = 4 a_\Omega b / c^2$ and $M_{\text{DM,inferred}} \propto b^2$ diagnostics.

Both modes report the zero-parameter **SZ-to-Lensing Hydrostatic Cluster Mass Bias** $(1 - b) = \Omega_{m,\text{void}} / \Omega_{m,\text{global}} = 0.8203$ (`-17.97%`).

---

## 🔬 Physical Foundation & First Principles

### 1. Unified Two-Ceiling Cover Pocket Architecture
* **Galaxy Floor ($M_B < 5 \times 10^{12}\ M_\odot$):**
  $$P_{\text{cap}} = 3 P_\delta = 8.5017 \times 10^{-13}\ \text{Pa} \implies g_{3P_\delta} = 3.7756 \times 10^{-11}\ \text{m/s}^2 \implies \Sigma_{p,\text{gal}} = 169.47\ M_\odot/\text{pc}^2$$
* **Isotropic Cluster Ceiling ($M_B \ge 5 \times 10^{12}\ M_\odot$):**
  $$P_{\text{cap}} = \frac{\Sigma_0}{3} = 2.0008 \times 10^{-10}\ \text{Pa} \implies g_{\text{iso}} = 5.7921 \times 10^{-10}\ \text{m/s}^2 \implies \Sigma_{p,\text{cluster}} = 2,599.82\ M_\odot/\text{pc}^2$$
* **Exact Cluster-to-Galaxy Hardware Bridge:**
  $$\frac{P_{\text{iso,cluster}}}{3 P_\delta} = \frac{\Sigma_0 / 3}{3\,\mathcal{C}_\delta\,\Sigma_0} = \frac{1}{9\mathcal{C}_\delta} = 235.337$$
* **Mass-Locked Core Radius:**
  $$r_0 = \sqrt{\frac{\eta_{\text{pocket}} M_B}{\pi\,\Sigma_p\,m_{\text{fac}}(x_{v,\max})}}, \qquad \eta_{\text{pocket}} = \frac{\Omega_{m,\text{void}}}{\Omega_{m,\text{global}} - \Omega_{m,\text{void}}} = 4.566048$$

### 2. Fast Analytic 2D Burkert Lensing Projection
Instead of slow nested `scipy.integrate.quad` loops, [`modules/data_loader.py`](../modules/data_loader.py) evaluates the exact closed-form projected factor $I(X) = \Sigma_{2D}(X) / (\rho_0 r_0)$ vectorized across $X = b / r_0$ (`< 1e-11` agreement with adaptive quadrature across `[1e-3, 10.0]`):
* `X < 1e-5`: $\pi / 2$
* `1e-5 <= X < 1 - 1e-4`: $\frac{\pi}{2s} + \frac{\ln(s+1)}{s} - \frac{\ln(1+q)}{q} - \ln X\left(\frac{1}{s} - \frac{1}{q}\right)$, with $s = \sqrt{1+X^2}$, $q = \sqrt{1-X^2}$
* `|X - 1| <= 1e-4`: smooth linear bridge around $I(1) = \frac{\pi}{2\sqrt{2}} + \frac{\ln(\sqrt{2}+1)}{\sqrt{2}} - 1$ with exact slope $-2/3$
* `X > 1 + 1e-4`: $\frac{\pi}{2s} + \frac{\operatorname{arctanh}(1/s)}{s} - \frac{\arccos(1/X)}{\sqrt{X^2-1}}$

Integrating $M_{2D,\text{pocket}}(<b) = 2\pi \int_0^b b'\Sigma_{2D}(b')\,db'$ via `cumulative_trapezoid` on a 2,000-point log grid and adding the compiled Hernquist baryonic contribution ($f_{\text{comp}} = 0.15$ for passive hot ICM clusters under the Two-Condition Compilation Gate; $f_{\text{comp}} = 1.0$ for galaxies) yields exact deflection angles, convergence profiles, and strong-lensing Einstein radii ($\bar{\kappa}(<R_E) = 1.0$) in `< 50 ms`.

### 3. Legacy Linear Scaffold Overlay ($a_\Omega$)
When running in Legacy mode, deflection follows:
$$\theta_{\text{total}}(b) = \frac{4 N_B(b)\mathcal{K}_\Omega}{c^2 b} + \frac{4 a_\Omega b}{c^2}, \qquad M_{\text{DM,inferred}}(b) = \frac{a_\Omega b^2}{\mathcal{K}_\Omega} m_p$$

---

## ⚙️ Architectural & Technical Implementation

* **Primary UI Driver:** `run_lensing_analysis(..., lensing_model="Cover Pocket (Linear K_Ω + Burkert)", ceiling_mode="Auto (3 P_δ Galaxy / Σ_0/3 Cluster)")`
* **Cover Pocket Lensing Engine:** `calculate_cover_pocket_lensing_profile(...)`
* **Presets (`LENS_PRESETS`):** `Bullet Cluster Subcluster (Bullet)` ($3.8 \times 10^{13}\ M_\odot$, $z_l = 0.296$), `Bullet Cluster Main (1E 0657-56)` ($1.1 \times 10^{14}\ M_\odot$, $z_l = 0.296$), `Abell 1689 Core (Cluster Lens)` ($1.5 \times 10^{14}\ M_\odot$, $z_l = 0.183$), and `Milky Way (Total Baryonic Footprint)` ($8.8 \times 10^{10}\ M_\odot$).

---

## 📊 Visual Diagnostic Output & Export Schema

* **Cover Pocket Mode Plots:**
  * **Panel 1 (Deflection Angle $\theta(b)$):** Plots Total, Burkert Cover Pocket, Compiled Baryons, and Legacy $a_\Omega$ Scaffold overlay with vertical $r_0$ and $R_E$ markers.
  * **Panel 2 (2D Convergence):** Plots local $\kappa(b)$ and enclosed mean $\bar{\kappa}(<b)$ against the $\kappa = 1.0$ Strong-Lensing Critical Line, highlighting $R_E$, $M_{2D}(<250\ \text{kpc})$, $(1-b) = 0.8203$, and the $1/(9\mathcal{C}_\delta)$ bridge.
* **CSV Export Schema (Cover Pocket Mode):** `impact_parameter_kpc`, `theta_total_arcsec`, `theta_baryon_arcsec`, `theta_pocket_arcsec`, `theta_legacy_scaffold_arcsec`, `kappa_local`, `kappa_mean_enclosed`, `M_2D_pocket_solar`, `M_2D_total_solar`.

---

## ⚡ Verified Benchmark Presets ($z_s = 1.0$, $H_{\text{global}} = 67.4275\ \text{km/s/Mpc}$)

| Lens Preset | Active Ceiling ($\Sigma_p$) | Locked $r_0$ [kpc] | Peak $\kappa(0)$ | Einstein Radius $R_E$ [kpc] | $M_{2D}(<250\ \text{kpc})$ [$M_\odot$] |
| :--- | :--- | ---: | ---: | ---: | ---: |
| **Bullet Cluster Subcluster (`Bullet`)** | $\Sigma_0/3$ (`2,599.82 M☉/pc²`) | `87.22` | `2.032` | `89.0` | `2.144e+14` |
| **Bullet Cluster Main (`1E 0657-56`)** | $\Sigma_0/3$ (`2,599.82 M☉/pc²`) | `148.40` | `2.259` | `152.6` | `3.695e+14` |
| **Abell 1689 Core** | $\Sigma_0/3$ (`2,599.82 M☉/pc²`) | `173.28` | `1.858` | `128.5` | `4.243e+14` |
| **Milky Way** | $3P_\delta$ (`169.47 M☉/pc²`) | `16.44` | `0.257` | `Sub-critical (< 1.0)` | `1.402e+12` |
