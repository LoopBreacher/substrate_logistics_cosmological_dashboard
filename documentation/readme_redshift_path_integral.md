# Redshift Path Integration: Line-of-Sight Velocity Accumulation & Supernova Residual Profiles

## 🌟 Executive Summary
The `redshift_path_integral` module computes first-principles line-of-sight path integrals of the local space-unspooling field $H_{\text{local}}(\vec{r})$ across multi-catalog volume horizons ($d \le 35\text{--}200\text{ Mpc}$). By ray-tracing photon paths through 3D compiled matter distributions, it calculates accumulated recession velocity $cz(d)$ and predicts distance modulus anomalies $\Delta\mu(d)$ along overdense filament corridors versus underdense void paths. This deterministic metric unspooling drag replaces empirical Dark Energy ($\Omega_\Lambda$) parameters with direct, zero-parameter observational forecasts for Rubin Observatory (LSST) and Euclid supernova surveys.

---

## 🔬 Physical Foundation & First Principles

* **Local Space-Unspooling Rate:**
  $$H_{\text{local}}(\vec{r}) = H_{\text{global}} + \sqrt{\frac{8\pi \mathcal{K}_\Omega \rho_N(\vec{r})}{3}}$$
  where $H_{\text{global}} = 67.42\text{ km/s/Mpc}$ represents the fundamental Planck vacuum floor, $\rho_N(\vec{r})$ is the 3D compiled proton node density ($\text{nodes/m}^3$), and $\mathcal{K}_\Omega = 1.11587 \times 10^{-37}\text{ m}^3/(\text{node}\cdot\text{s}^2)$ is the Substrate Anchor Constant.

* **Line-of-Sight Redshift Accumulation Integral:**
  $$cz(d, \hat{n}) = \int_0^d H_{\text{local}}(s \hat{n}) \, ds$$
  accumulates total expansion velocity ($cz$) along celestial unit ray vector $\hat{n} = (\text{RA}, \text{Dec})$ out to distance $d$.

* **Path-Integrated Effective Expansion Rate:**
  $$H_{\text{eff}}(d, \hat{n}) = \frac{cz(d, \hat{n})}{d}$$
  Because matter density is everywhere non-negative ($\rho_N \ge 0$), $H_{\text{eff}}$ provides a strict physical lower bound ($H_{\text{eff}} \ge H_{\text{global}}$).

* **Supernova Distance Modulus Residual Anomaly ($\Delta\mu$):**
  $$\Delta\mu(d, \hat{n}) = \mu_{\text{true}} - \mu_{\text{naive}} = -5 \log_{10}\left( \frac{H_{\text{eff}}(d, \hat{n})}{H_{\text{global}}} \right)$$
  Quantifies apparent Type Ia supernova magnitude brightness variations relative to an ideal homogeneous vacuum universe ($H_{\text{global}}$). Sightlines traversing dense filaments yield $H_{\text{eff}} > H_{\text{global}}$, generating negative distance modulus residuals ($\Delta\mu < 0$, appearing brighter than naive $cz \to d$ inversion implies).

---

## ⚙️ Architectural & Technical Implementation

* **Primary Function:** `run_redshift_path_analysis(...)`
* **Input Parameters:**
  * `fil_ra`, `fil_dec` (*float*): Filament target direction in degrees (Default: $\text{RA} = 201.0^\circ, \text{Dec} = -29.8^\circ$, Centaurus / Local Sheet corridor).
  * `void_ra`, `void_dec` (*float*): Void target direction in degrees (Default: $\text{RA} = 280.0^\circ, \text{Dec} = -20.0^\circ$, Local Void corridor).
  * `tree` (*scipy.spatial.cKDTree*): Indexed 3D galaxy catalog spatial tree.
  * `gal_positions` (*np.ndarray*): Cartesian position matrix $(N, 3)$ in Mpc.
  * `nodes_per_galaxy` (*np.ndarray*): Integer proton node counts ($N = M / m_p$).
  * `max_dist_mpc` (*float*): Maximum ray-tracing depth limit in Mpc ($10.0\text{--}200.0\text{ Mpc}$).
  * `sigma_mpc` (*float*): Gaussian kernel field smoothing scale $\sigma$ ($0.5\text{--}5.0\text{ Mpc}$).
  * `h_global` (*float*): Vacuum floor expansion rate ($67.42\text{ km/s/Mpc}$).
* **Algorithmic Pipeline:**
  1. **Ray Vector Synthesis:** Converts input celestial coordinates $(\text{RA}, \text{Dec})$ into 3D unit direction vectors $\hat{n}_{\text{fil}}$ and $\hat{n}_{\text{void}}$.
  2. **Field Sampling:** Evaluates $H_{\text{local}}(s \hat{n})$ at $N_{\text{steps}} = 250$ discrete points along each path via 3D KD-tree density queries.
  3. **Numerical Integration:** Performs cumulative Simpson rule path integration $\int_0^s H_{\text{local}}(s') ds'$ to calculate $cz(s)$ profiles.
  4. **Residual Derivation:** Computes true distance modulus $\mu_{\text{true}} = 5 \log_{10}(d \cdot 10^5)$ versus naive inversion distance modulus $\mu_{\text{naive}} = 5 \log_{10}((cz/H_{\text{global}}) \cdot 10^5)$ to extract $\Delta\mu(d)$.
  5. **Checkpoint Interpolation:** Extracts exact distance evaluation checkpoints ($5.0$, $10.0$, $15.0$, $20.0\text{ Mpc}$, etc.) for tabular output export.

---

## 📊 Visual Diagnostic Output & Graphics

* **Panel 1 (Redshift Accumulation Path Integral):**
  Plots accumulated expansion velocity $cz(d)$ along Filament (solid red line) and Void (dashed blue line) corridors against homogeneous Planck $H_0 = 67.42\text{ km/s/Mpc}$ (dotted gray) and SH0ES $H_0 = 73.04\text{ km/s/Mpc}$ (dash-dotted green) baselines. Exposes how velocity builds up faster along overdense structures.
* **Panel 2 (Euclid / Rubin LSST Predicted Residual Anomaly):**
  Plots distance modulus residual $\Delta\mu(d)$ along both sightlines relative to the $\Delta\mu = 0.0$ baseline. Demonstrates how void corridors remain flat near $\Delta\mu \approx 0.0\text{ mag}$ while filament corridors dive into negative residual troughs ($\Delta\mu \approx -0.05\text{ to }-0.16\text{ mag}$), directly predicting directional supernova brightness anomalies.

---

## 🎯 KPI Metrics & Export Deliverables

| Metric Key | Unit / Format | Physical Meaning |
| :--- | :--- | :--- |
| `Filament Peak H_eff` | `km/s/Mpc` | Maximum cumulative effective expansion rate along filament path. |
| `Void Floor H_eff` | `km/s/Mpc` | Minimum cumulative effective expansion rate along void path ($H_{\text{eff}} \to H_{\text{global}}$). |
| `Max Filament Δμ` | `mag` | Maximum negative distance modulus residual trough along filament path. |
| `Max Path Velocity Δcz` | `km/s` | Peak velocity deficit accumulated between filament and void sightlines ($(cz_{\text{fil}} - cz_{\text{void}})_{\max}$). |

* **CSV Export Schema (`redshift_path_integration_profiles.csv`):**
  * `Distance [Mpc]`: Checkpoint distance $d$ along ray path.
  * `Filament cz [km/s]`: Accumulated velocity along filament corridor.
  * `Filament H_eff [km/s/Mpc]`: Effective path expansion rate $cz_{\text{fil}} / d$.
  * `Filament Δμ [mag]`: Supernova distance modulus anomaly $\mu_{\text{true}} - \mu_{\text{naive}}$.
  * `Void cz [km/s]`: Accumulated velocity along void corridor.
  * `Void H_eff [km/s/Mpc]`: Effective path expansion rate $cz_{\text{void}} / d$.
  * `Void Δμ [mag]`: Supernova distance modulus anomaly along void path.
  * `Velocity Contrast [km/s]`: Accumulated differential velocity $cz_{\text{fil}} - cz_{\text{void}}$.

---

## ⚡ Performance Benchmarks & Catalog Behavioral Notes

* **Execution Benchmarks (Dual 250-Step Ray Paths):**
  * **UNG 2013 ($\le 35\text{ Mpc}$, ~3k galaxies):** `~0.02s`
  * **Cosmicflows-4 ($\le 150\text{ Mpc}$, ~35k galaxies):** `~0.10s`
  * **2M++ Infrared ($\le 200\text{ Mpc}$, ~70k galaxies):** `~0.22s`
* **Field Smoothing Scale ($\sigma$) Sensitivity:**
  * **Fine Resolution ($\sigma \le 1.2\text{ Mpc}$):** Resolves sharp, localized cluster core density spikes, creating steep localized dips in $\Delta\mu$ ($\le -0.15\text{ mag}$).
  * **Standard Baseline ($\sigma = 1.8\text{--}2.0\text{ Mpc}$):** Smoothly bridges intergalactic filaments while preserving void-wall density contrast.
  * **Coarse Resolution ($\sigma \ge 5.0\text{ Mpc}$):** Smooths cluster density into background voids, dampening peak velocity contrast between sightlines.