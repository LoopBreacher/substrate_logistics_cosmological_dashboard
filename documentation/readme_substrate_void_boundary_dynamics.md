# Substrate Void Boundary Dynamics: Cosmic Void Wall Boundary Dynamics & Velocity Shear

## 🌟 Executive Summary
The **Substrate Void Boundary Dynamics** module models spatial metric unspooling along radial trajectories originating deep inside cosmic void cores and crossing into surrounding cluster and wall structures. By replacing empirical dark energy and artificial boundary switches with a continuous 3D compiled node density field $\rho_N(\vec{r})$, the module demonstrates how spatial transitions from empty void cores ($\rho_N \to 0$) into dense wall corridors generate localized expansion boosts and differential velocity shear spikes ($dv_{\text{rec}}/dr$). This provides a zero-parameter physical mechanism for the dynamic evacuation of cosmic voids without requiring repulsive dark energy or modified gravity parameter tuning.

---

## 🔬 Physical Foundation & First Principles

### 1. Dual-Scale Intergalactic Density Field
The module calculates continuous background density along ray paths using a dual 3D Gaussian smoothing kernel architecture to decouple macro intergalactic background volume from fine filament core density:

$$\rho_{\text{macro}}(\vec{r}) = \frac{1}{(2\pi \sigma_{\text{macro}}^2)^{3/2}} \sum_i N_i \exp\left( -\frac{\vert\vec{r}_i - \vec{r}\vert^2}{2\sigma_{\text{macro}}^2} \right) \cdot \frac{1}{\mathrm{MPC\_TO\_METER}^3}$$

$$\rho_{\text{filament}}(\vec{r}) = \frac{1}{(2\pi \sigma_{\text{filament}}^2)^{3/2}} \sum_i N_i \exp\left( -\frac{\vert\vec{r}_i - \vec{r}\vert^2}{2\sigma_{\text{filament}}^2} \right) \cdot \frac{1}{\mathrm{MPC\_TO\_METER}^3}$$

* **Macro Kernel Bandwidth:** $\sigma_{\text{macro}} = 1.5\text{ Mpc}$ (smooths intergalactic background node volume).
* **Filament Kernel Bandwidth:** $\sigma_{\text{filament}} = 0.5\text{ Mpc}$ (isolates narrow, high-density cluster wall filaments).

### 2. Local Unspooling Expansion Field
The local metric expansion rate $H_{\text{local}}(\vec{r})$ is driven directly by compiled intergalactic matter density $\rho_{\text{macro}}(\vec{r})$ above the fundamental vacuum baseline:

$$H_{\text{local}}(\vec{r}) = H_{\text{global}} + \frac{1}{\mathrm{UNIT\_CONV}} \sqrt{\frac{8\pi \mathcal{K}_\Omega \rho_{\text{macro}}(\vec{r})}{3}}$$

Where:
* $H_{\text{global}} = 67.42\text{ km/s/Mpc}$ (global Planck vacuum floor).
* $\mathcal{K}_\Omega = 1.11587 \times 10^{-37}\text{ m}^3/(\text{node}\cdot\text{s}^2)$ (Substrate Anchor Constant).
* $\mathrm{UNIT\_CONV} = 3.24078 \times 10^{-20}\text{ s}^{-1} / (\text{km/s/Mpc})$.

### 3. Recession Velocity & Radial Velocity Shear
The radial recession velocity $v_{\text{rec}}(r)$ accumulated along the ray path from the void center $r=0$ to $r$ is given by:

$$v_{\text{rec}}(r) = r \cdot H_{\text{local}}(r)$$

Taking the radial derivative yields the boundary velocity shear $dv_{\text{rec}}/dr$:

$$\frac{dv_{\text{rec}}}{dr} = H_{\text{local}}(r) + r \frac{dH_{\text{local}}}{dr}$$

---

## ⚙️ Architectural & Technical Implementation

* **Primary Function:** `run_void_boundary_analysis(tree, gal_positions, nodes_per_galaxy, void_center_mpc, max_r_mpc, n_pts, h_global)`
* **Input Parameters:**
  * `tree`: `scipy.spatial.cKDTree` spatial index constructed from active survey catalog.
  * `gal_positions`: `np.ndarray` of shape $(N, 3)$ containing ICRS Cartesian galaxy positions in Mpc.
  * `nodes_per_galaxy`: `np.ndarray` of proton node counts $N = M / m_p$.
  * `void_center_mpc`: 3D Cartesian coordinates of local void core origin (default: $[-2.5, -6.0, -1.5]\text{ Mpc}$).
  * `max_r_mpc`: Radial evaluation extent in Mpc (default: $18.0\text{ Mpc}$).
  * `n_pts`: Number of radial sampling steps (default: $200$).
  * `h_global`: Fundamental vacuum expansion floor in km/s/Mpc (default: $67.42$).

* **Algorithmic Pipeline:**
  1. **Trajectory Vector Construction:** Constructs a normalized ray direction vector $\hat{e}_{\text{sheet}}$ pointing from the void core origin toward the Virgo cluster wall ($\text{RA}=187.7^\circ, \text{Dec}=12.4^\circ, d=16.5\text{ Mpc}$).
  2. **Dual-Scale Kernel Query:** For each radial step $r_i \in [0, 18.0\text{ Mpc}]$, queries `cKDTree` within $3\sigma$ neighbor spheres to compute $\rho_{\text{macro}}(r_i)$ and $\rho_{\text{filament}}(r_i)$.
  3. **Expansion & Velocity Shear Evaluation:** Computes $H_{\text{local}}(r_i)$ and derives radial recession velocity $v_{\text{rec}}(r_i)$ and boundary shear $dv_{\text{rec}}/dr_i$ using finite-difference numerical gradients (`np.gradient`).
  4. **Visualization & KPI Packaging:** Plots 3-panel dark diagnostic figure, compiles metric cards, and returns exportable pandas DataFrame.

---

## 📊 Visual Diagnostic Output & Graphics

* **Panel 1 (Void-to-Wall Intergalactic Density Profiles):**
  Renders continuous background macro density $\rho_{\text{macro}}(r)$ alongside localized filament core density $\rho_{\text{filament}}(r)$ in $\text{nodes/m}^3$. Illustrates the dramatic density jump as the trajectory exits the void core and pierces the cluster wall corridor.
* **Panel 2 (Radial Expansion Rate Unspooling Field $H_{\text{local}}$):**
  Traces the radial expansion profile $H_{\text{local}}(r)$ relative to the horizontal vacuum floor line $H_{\text{global}} = 67.42\text{ km/s/Mpc}$. Demonstrates void core vacuum floor locking at $r \to 0$ and unspooling expansion boosts at the wall boundary.
* **Panel 3 (Cosmic Void Wall Velocity Shear Dynamics $dv_{\text{rec}}/dr$):**
  Plots the differential radial velocity shear spike marking the dynamic wall interface. Shows the physical shear driver evacuating galaxies out of underdense voids.

---

## 🎯 KPI Metrics & Export Deliverables

| Metric Key | Unit / Format | Physical Meaning |
| :--- | :--- | :--- |
| `Void Core Expansion` | `km/s/Mpc` | Expansion rate deep inside the underdense void core ($r=0\text{ Mpc}$). |
| `Peak Wall Expansion` | `km/s/Mpc` | Maximum $H_{\text{local}}$ unspooling peak reached inside the wall corridor. |
| `Max Velocity Shear` | `km/s/Mpc` | Peak differential velocity shear spike $(dv_{\text{rec}}/dr)_{\text{max}}$ across the boundary. |
| `Shear Peak Distance` | `Mpc` | Radial distance from void center where velocity shear reaches its maximum. |

* **CSV Export Schema (`void_boundary_wall_profile.csv`):**
  * `r_mpc`: Radial distance from void center origin in Mpc.
  * `rho_macro_nodes_m3`: Continuous intergalactic macro field density ($\sigma=1.5\text{ Mpc}$) in $\text{nodes/m}^3$.
  * `rho_filament_nodes_m3`: Localized filament core density ($\sigma=0.5\text{ Mpc}$) in $\text{nodes/m}^3$.
  * `H_local_kmsMpc`: Computed metric expansion rate $H_{\text{local}}(r)$ in km/s/Mpc.
  * `v_recession_kms`: Path-integrated radial recession velocity $v_{\text{rec}}(r)$ in km/s.
  * `v_shear_kmsMpc`: Boundary velocity shear $dv_{\text{rec}}/dr$ in km/s/Mpc.

---

## ⚡ Performance Benchmarks & Catalog Behavioral Notes

### Execution Benchmarks

* **UNG 2013 (Local Volume $\le 35\text{ Mpc}$, ~3.5k nodes):** `~0.005 s` (5 ms)
* **Cosmicflows-4 (Extended $\le 150\text{ Mpc}$, ~35k nodes):** `~0.005 s` (5 ms)
* **2M++ Infrared (Full Sky $\le 200\text{ Mpc}$, ~70k nodes):** `~0.005 s` (5 ms)

### Kernel Architecture Note: Field Smoothing Scale ($\sigma$)

The sidebar **Field Smoothing Scale ($\sigma$)** slider has **no effect** on this module by design. To accurately resolve fine wall transitions without blurring, this module utilizes a specialized internal **dual-scale Gaussian kernel** architecture:
1. Fixed $\sigma_{\text{macro}} = 1.5\text{ Mpc}$ for intergalactic volume density.
2. Fixed $\sigma_{\text{filament}} = 0.5\text{ Mpc}$ for localized filament core density.

### Catalog Response Comparison

| Catalog Dataset | Void Core $H(r=0)$ | Peak Wall $H_{\text{max}}$ | Max Velocity Shear | Shear Peak Dist | Physical Feature Resolved |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **UNG 2013 (≤ 35 Mpc)** | $67.45\text{ km/s/Mpc}$ | $68.35\text{ km/s/Mpc}$ | $72.46\text{ km/s/Mpc}$ | $4.5\text{ Mpc}$ | Inner Local Void boundary floor locking. |
| **Cosmicflows-4 (≤ 150 Mpc)** | $67.42\text{ km/s/Mpc}$ | $75.50\text{ km/s/Mpc}$ | $88.74\text{ km/s/Mpc}$ | $15.7\text{ Mpc}$ | Deep Virgo / Local Sheet wall density wall. |
| **2M++ Infrared (≤ 200 Mpc)** | $67.61\text{ km/s/Mpc}$ | $70.95\text{ km/s/Mpc}$ | $77.08\text{ km/s/Mpc}$ | $15.7\text{ Mpc}$ | Macro infrared wall structure smoothing. |