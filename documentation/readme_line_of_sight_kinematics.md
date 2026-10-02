# Line-of-Sight Kinematics & Peculiar Velocity Solver

## 🌟 Executive Summary
This module performs line-of-sight ray integration across the 3D cosmic web to decouple true cosmological expansion ($c z_{\text{exp}}$) from local peculiar velocity flows ($v_{\text{pec}}$). Because local distances ($d < 30\text{ Mpc}$) are determined via gold-standard standard candles (TRGB or Cepheids), attempting to infer distance directly from redshift is fundamentally flawed. This solver uses benchmark physical distances to accurately extract peculiar velocities across cluster infalls and void outflows.

---

## 🔬 Physical Foundation & First Principles

* **Line-of-Sight Expansion Path Integral:**
  $$c z_{\text{exp}} = \int_0^{d_{\text{benchmark}}} H_{\text{local}}(s \hat{n}) \, ds$$

* **Derived Line-of-Sight Peculiar Velocity ($v_{\text{pec}}$):**
  $$v_{\text{pec}} = c z_{\text{obs}} - c z_{\text{exp}}$$

* **Substrate Mechanics:**
  * **Density-Dependent Unspooling:** Photon paths traversing overdense cosmic web filaments encounter elevated compiled node density $\rho_N(s)$, increasing $H_{\text{local}}(s)$ along the ray path.
  * **Kinematic Decomposition:** Rather than assuming a uniform Hubble parameter ($c z = H_0 \cdot d$), the line-of-sight integral accounts for density fluctuations along the line of sight, allowing pure expansion ($c z_{\text{exp}}$) to be isolated from peculiar dynamics ($v_{\text{pec}}$).

---

## ⚙️ Architectural & Technical Implementation

* **Primary Function:** `run_line_of_sight_kinematics_analysis(...)`
* **Core Integrator:** `solve_line_of_sight_kinematics(...)`
* **Input Parameters:**
  * `ra_deg`, `dec_deg`: Target sky coordinates in Right Ascension and Declination (ICRS degrees).
  * `cz_obs`: Observed recessional velocity ($c z_{\text{obs}}$ in $\text{km/s}$).
  * `d_known`: Physical benchmark distance from TRGB or Cepheid measurements ($\text{Mpc}$).
  * `tree`: `scipy.spatial.cKDTree` catalog index.
  * `sigma_mpc`: Gaussian field smoothing scale ($\sigma = 1.8\text{ Mpc}$ default).
  * `h_global`: Vacuum floor Hubble constant ($H_{\text{global}} = 67.42\text{ km/s/Mpc}$).

---

## 📊 Visual Diagnostic Output & Graphics

* **Panel 1 (Line-of-Sight Expansion Profile):** Maps the local expansion rate $H_{\text{local}}(s)$ along the sightline from observer ($s=0$) out to the integration boundary, highlighting TRGB/Cepheid benchmark distances.
* **Panel 2 (Kinematic Decomposition Profile):** Compares accumulated Substrate unspooling $c z_{\text{exp}}(s)$ against the linear baseline ($c z = H_{\text{global}} \cdot d$) and observed velocity $c z_{\text{obs}}$, explicitly marking the derived peculiar velocity vector $v_{\text{pec}}$.

---

## 🎯 KPI Metrics & Export Deliverables

| Metric Key | Unit / Format | Physical Meaning |
| :--- | :--- | :--- |
| `Substrate Expansion (cz_exp)` | `km/s` | Path-integrated cosmological expansion velocity along sightline. |
| `Derived Peculiar Velocity (v_pec)` | `km/s` | Net line-of-sight peculiar motion ($c z_{\text{obs}} - c z_{\text{exp}}$). |
| `TRGB/Cepheid Benchmark` | `Mpc` | Primary physical distance baseline. |
| `Linear Hubble Distance` | `Mpc` | Uncorrected naive distance estimate ($c z_{\text{obs}} / H_{\text{global}}$). |

* **CSV Export Schema (`line_of_sight_kinematics_benchmark_comparison.csv`):**
  * `Galaxy Target`: Target name or celestial coordinate label.
  * `Environment`: Local cosmic web environment (e.g., Dense Cluster Core, Deep Local Void).
  * `Observed cz [km/s]`: Total spectroscopic recessional velocity.
  * `TRGB/Cepheid Benchmark [Mpc]`: Standard candle distance.
  * `Substrate Expansion cz_exp [km/s]`: Pure unspooling expansion contribution.
  * `Derived v_pec [km/s]`: Inward infall ($<0$) or outward recession ($>0$).
  * `Linear Hubble Distance [Mpc]`: Simple $c z / H_{\text{global}}$ distance.

---

## ⚡ Performance Benchmarks & Catalog Behavioral Notes

* **Execution Speed:** `~0.05s` per line-of-sight ray calculation across Local Volume catalogs ($N=3,000$ to $70,000$ galaxies).
* **Environment Sensitivity:**
  * **Inward Cluster Infall ($v_{\text{pec}} < 0$):** In dense environments like NGC 6946 ($v_{\text{pec}} = -220.4\text{ km/s}$), line-of-sight unspooling exceeds observed velocity due to local gravitational deceleration.
  * **Outward Recessional Flow ($v_{\text{pec}} > 0$):** Target galaxies embedded in cluster outward streams (e.g., NGC 1365 in Fornax, $v_{\text{pec}} = +478.5\text{ km/s}$) exhibit positive peculiar flow vectors.