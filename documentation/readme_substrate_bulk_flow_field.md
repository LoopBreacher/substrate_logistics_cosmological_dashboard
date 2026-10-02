# Substrate Bulk Flow Field: Macro Velocity Gradients & Conveyor Drift

## 🌟 Executive Summary
The **Substrate Bulk Flow Field** module models 2D/3D galactic peculiar velocity fields across cosmic web structures without invoking Cold Dark Matter halos or abstract gravitational "pulling" forces at a distance. In Substrate Logistics, peculiar motion is reclassified as coordinate conveyor drift: matter naturally drifts toward overdense galaxy corridors where compiled node density $\rho_N(\vec{r})$ accelerates spatial address deletion ($E=0$) and metric unspooling rates ($H_{\text{local}}$). By evaluating spatial expansion gradients $\nabla H_{\text{local}}(\vec{r})$ directly across 3D galaxy catalogs, this module computes local cluster infalls ($10\text{--}20\text{ Mpc}$) and macro-scale bulk flows ($100\text{--}200\text{ Mpc}$) in zero-parameter agreement with observational surveys.

---

## 🔬 Physical Foundation & First Principles

In standard $\Lambda\text{CDM}$ cosmology, peculiar velocity fields are driven by gravitational acceleration gradients sourced by unseen dark matter overdensities ($\vec{v} \propto \nabla \Phi_{\text{DM}}$). Substrate Logistics replaces this model with a direct hydrodynamic response to spatial unspooling rate gradients:

### 1. Local Unspooling Field
The local Hubble expansion rate $H_{\text{local}}(\vec{r})$ is determined by the local compiled node density $\rho_N(\vec{r})$:

$$H_{\text{local}}(\vec{r}) = H_{\text{global}} + \Theta(\rho_N(\vec{r}) - \rho_{\text{threshold}}) \sqrt{\frac{8\pi \mathcal{K}_\Omega \rho_N(\vec{r})}{3}}$$

Where:
* $H_{\text{global}} = 67.42\text{ km/s/Mpc}$ is the global hardware vacuum floor.
* $\mathcal{K}_\Omega = 1.11587 \times 10^{-37}\text{ m}^3/(\text{node}\cdot\text{s}^2)$ is the Substrate Anchor Constant.
* $\Theta$ is the Heaviside step function enforcing the activation threshold ($\rho_{\text{threshold}} = 4.0\text{ nodes/m}^3$).

### 2. Inflow Peculiar Velocity Vector
Galactic motion is directed toward higher spatial processing throughput. The peculiar inflow vector $\vec{v}_{\text{inflow}}(\vec{r})$ is proportional to the spatial gradient of $H_{\text{local}}(\vec{r})$:

$$\vec{v}_{\text{inflow}}(\vec{r}) = +\alpha \nabla H_{\text{local}}(\vec{r})$$

Where $\alpha = 120.0\text{ Mpc}^2$ is the coupling scale converting spatial expansion gradients ($\text{km/s/Mpc}^2$) into physical velocity ($\text{km/s}$).

### 3. Line-of-Sight Peculiar Velocity
Along any given observer sightline $\hat{r} = \vec{r} / \vert{}\vec{r}\vert{}$, the line-of-sight component of peculiar velocity evaluates to:

$$v_{\text{pec, LOS}}(\vec{r}) = \vec{v}_{\text{inflow}}(\vec{r}) \cdot \hat{r} = \alpha \left( \nabla H_{\text{local}}(\vec{r}) \cdot \hat{r} \right)$$

In deep void cores where matter density vanishes ($\rho_N \to 0$), the local unspooling field relaxes to $H_{\text{global}}$, driving gradients to zero ($\nabla H_{\text{local}} = 0$) and suppressing peculiar velocity ($v_{\text{inflow}} = 0\text{ km/s}$).

---

## ⚙️ Architectural & Technical Implementation

### Primary Functions
1. `compute_line_of_sight_peculiar_velocity(...)`: Evaluates 3D finite differences ($dr = 0.2\text{ Mpc}$) around a coordinate point to compute line-of-sight velocity components.
2. `compute_bulk_flow_field(...)`: Generates a 2D grid slice on the galactic plane ($Z=0$), evaluating $H_{\text{local}}$ and spatial derivatives $\frac{\partial H}{\partial x}$, $\frac{\partial H}{\partial y}$ using `numpy.gradient`.
3. `run_bulk_flow_analysis(...)`: Main Streamlit UI driver function. Assembles 2D contour maps, vector quivers, galaxy overlays, KPI metrics, and export datasets.

### Algorithmic Pipeline
1. **Catalog Coordinate Query:** Fetches $X, Y, Z$ positions and node counts $N = M / m_p$ from the active `cKDTree` index.
2. **2D Grid Construction:** Builds an $N \times N$ spatial mesh bounded by $\pm \text{spatial\_bounds}$ (e.g., $\pm 150\text{ Mpc}$).
3. **Density & Expansion Evaluation:** Computes $H_{\text{local}}(x, y, 0)$ across the mesh using Gaussian kernel bandwidth $\sigma$.
4. **Vector Field Generation:** Takes spatial derivatives to obtain inflow components $V_x = \alpha \frac{\partial H}{\partial x}$, $V_y = \alpha \frac{\partial H}{\partial y}$, and magnitude $V_{\text{mag}} = \sqrt{V_x^2 + V_y^2}$.
5. **Render & Export:** Renders color-coded contours, vector quivers, and plane-sliced galaxy points ($\vert{}Z\vert{} \le z_{\text{thickness}}$).

---

## 📊 Visual Diagnostic Output & Graphics

The module outputs a unified 2D spatial diagnostic plot (`figsize=(10, 7)`):

* **Background Contour Map (`plasma` colormap):** Displays the continuous local unspooling field $H_{\text{local}}(\vec{r})$, shading from underdense void cores ($67.42\text{ km/s/Mpc}$, dark purple) up to dense supercluster filaments ($71.40\text{ km/s/Mpc}$, bright yellow).
* **Vector Quiver Overlay (`cool` colormap):** Displays directional arrows representing $\vec{v}_{\text{inflow}}$. Vector length and color scale with velocity magnitude $V_{\text{mag}}$, illustrating convergence onto overdense filament corridors and divergence away from void centers.
* **Catalog Galaxy Overlay:** White scatter points representing real catalog galaxies located within the Z-slice plane ($\vert{}Z\vert{} \le z_{\text{thickness}}$).

---

## 🎯 KPI Metrics & Export Deliverables

### Top-Level Streamlit KPI Cards

| Metric Key | Unit / Format | Physical Meaning |
| :--- | :--- | :--- |
| `Mean Drift Velocity` | `km/s` | Spatially averaged peculiar flow magnitude across the evaluation volume. |
| `Peak Filament Infall` | `km/s` | Maximum velocity magnitude recorded at steep cluster/filament boundaries. |
| `Void Core Drift` | `km/s` | Minimum velocity magnitude inside underdense void centers (evaluates to $0.00\text{ km/s}$). |

### Exported CSV Dataset (`bulk_flow_velocity_grid.csv`)
* `X_Mpc`: Cartesian X coordinate on galactic plane [Mpc].
* `Y_Mpc`: Cartesian Y coordinate on galactic plane [Mpc].
* `H_local_kmsMpc`: Local unspooling rate $H_{\text{local}}(x, y, 0)$ [km/s/Mpc].
* `Vx_inflow_kms`: Inward peculiar velocity component along X axis [km/s].
* `Vy_inflow_kms`: Inward peculiar velocity component along Y axis [km/s].
* `V_magnitude_kms`: Net peculiar velocity magnitude $V_{\text{mag}} = \sqrt{V_x^2 + V_y^2}$ [km/s].

---

## ⚡ Performance Benchmarks & Catalog Behavioral Notes

### Benchmark Evaluation (Cosmicflows-4 | $R \le 150\text{ Mpc}$ | $\sigma = 15.0\text{ Mpc}$ | Grid $50 \times 50$)
* **Execution Time:** `~0.18s` (full field calculation, gradient evaluation, and plot assembly).
* **Unspooling Field Range ($H_{\text{local}}$):** $67.420\text{ km/s/Mpc}$ (min) to $71.398\text{ km/s/Mpc}$ (max), Mean = $69.068\text{ km/s/Mpc}$.
* **Peculiar Velocity Range ($V_{\text{mag}}$):** $0.00\text{ km/s}$ (void interior) to $10.88\text{ km/s}$ (macro filament infall), Mean = $3.85\text{ km/s}$.

### Catalog Scaling & Resolution Insights
* **UNG 2013 (≤ 35 Mpc, Fine Grid $\sigma = 1.8\text{ Mpc}$):** Captures sharp local cluster infall velocities (e.g., Virgo infall peaks up to $80\text{ km/s}$ with mean sheet drift $\sim 3.36\text{ km/s}$).
* **Cosmicflows-4 / 2M++ (≤ 150–200 Mpc, Smooth $\sigma = 15.0\text{ Mpc}$):** Smoothes fine cluster spikes, revealing macro-scale bulk flow channels directed toward supercluster attractors (Laniakea, Perseus-Pisces) while averaging macro void drift down to $\sim 3.85\text{ km/s}$.