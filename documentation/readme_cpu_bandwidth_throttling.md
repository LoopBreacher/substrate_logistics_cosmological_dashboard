# CPU Bandwidth Throttling & Gravitational Time Dilation Engine (`cpu_bandwidth_throttling.py`)

## 🌟 Executive Summary
The **CPU Bandwidth Throttling & Gravitational Time Dilation Engine** replaces standard General Relativity's continuous spacetime curvature with hardware metric processing latency[cite: 38, 39]. Under Logistical Relativity, time is the experiential observation of localized metric state updating[cite: 39]. Governed by the Pythagorean Bandwidth Allocation ($c^2 = v^2 + u^2$), external spatial velocity ($v_\Omega$) and gravitational address deletion ($g_\Omega$) consume the grid's finite processing capacity[cite: 38, 39]. This throttles the internal update rate ($u_{\text{clock}}$), dephasing atomic clocks without physical time bending[cite: 38, 39]. At the Event Horizon ($R_{\text{EH}}$), $g_\Omega^2 = 1.0$, dropping residual update capacity to zero ($u = 0.0000$) and triggering a $100\%$ CPU Kernel Lock[cite: 38, 39].

---

## 🔬 Physical Foundation & First Principles

### 1. Primary Equations

* **Pythagorean Bandwidth Allocation ($c^2 = v^2 + u^2$):**
  $$c^2 = v^2 + u^2 \implies u = c \sqrt{1 - \frac{v^2}{c^2}}$$

* **Normalized Dynamic Hardware Loads ($g_\Omega^2, v_\Omega^2$):**
  $$g_\Omega^2 = \frac{v_{\text{escape}}^2}{c^2} = \frac{2 N \mathcal{K}_\Omega}{c^2 R}$$
  $$v_\Omega^2 = \frac{v^2}{c^2}$$

* **Observable Clock Residual Capacity ($u_{\text{clock}}$):**
  $$u_{\text{clock}} = \sqrt{1.0 - (g_\Omega^2 + v_\Omega^2)}$$

* **Absolute Frequency Shift ($\Delta f / f_0$) vs. Unthrottled Deep Space ($u = 1.0$):**
  $$\left(\frac{\Delta f}{f_0}\right)_{\text{space}} = u_{\text{clock}} - 1.0$$

* **Relative Frequency Shift ($\Delta f / f_0$) vs. Earth Ground ($u_{\text{Earth}}$):**
  $$\left(\frac{\Delta f}{f_0}\right)_{\text{Earth}} = \frac{u_{\text{target}} - u_{\text{Earth}}}{u_{\text{Earth}}}$$

* **Daily Clock Drift [$\mu\text{s/day}$]:**
  $$\text{Drift}_{\text{daily}} = \left(\frac{\Delta f}{f_0}\right) \times 8.64 \times 10^{10}\ \mu\text{s/day}$$

### 2. Reference Frame Disambiguation (Deep Space vs. Earth Ground)

* **Deep Space Baseline ($r \to \infty, v=0$):** Both Earth ground and GPS satellite clocks experience absolute metric throttling relative to unthrottled deep space ($\Delta f_{\text{space}} < 0$). For GPS orbit ($20,180\text{ km}$ altitude), $\Delta f_{\text{space}} / f_0 = -2.5040 \times 10^{-10}$ (losing $-21.63\,\mu\text{s/day}$).
* **Earth Ground Baseline ($r = R_\oplus$):** Because Earth ground clocks are deeper in Earth's gravitational potential ($\Delta f_{\text{space}} = -6.9701 \times 10^{-10}$), the GPS clock experiences less throttling than the ground clock. Relative to an Earth ground observer, the GPS clock runs fast by $+4.4661 \times 10^{-10}$ ($+38.59\,\mu\text{s/day}$).

---

## ⚙️ Architectural & Technical Implementation

* **Primary Functions:** `compute_bandwidth_throttling_point(...)`, `compute_bandwidth_throttling_profile(...)`, & `run_cpu_bandwidth_throttling_analysis(...)`
* **Input Parameters:**
  * `preset_key` (`str`): Target benchmark preset (`"Earth Surface"`, `"GPS Satellite Orbit"`, `"ISS Orbit"`, `"Solar Photosphere"`, `"Sirius B"`, `"Crab Pulsar"`, `"Sagittarius A*"`, or `"Custom"`)[cite: 38].
  * `observer_mode` (`str`): `"Static (v = 0)"` or `"Orbital (v = v_orbit)"`[cite: 38].
  * `alpha_tension` (`float`, default `0.0`): Static internal cell structural tension slider[cite: 38].
* **Algorithmic Pipeline:**
  1. **Dynamic Load Evaluation:** Solves $g_\Omega^2 = \frac{2 N \mathcal{K}_\Omega}{c^2 R}$ and $v_\Omega^2 = \frac{v^2}{c^2}$[cite: 38].
  2. **Clock Rate Extraction:** Evaluates $u_{\text{clock}} = \sqrt{1.0 - (g_\Omega^2 + v_\Omega^2)}$ and frequency shifts $\Delta f_{\text{space}}$ and $\Delta f_{\text{Earth}}$[cite: 38].
  3. **Radial Profile Generation:** Computes a 1D radial profile across $N_{\text{pts}} = 250$ spatial points[cite: 38].
  4. **Ground Baseline Calibration:** Compares target clock drift against Earth ground clock baseline to derive net relative drift ($+38.59\,\mu\text{s/day}$ for GPS)[cite: 38].

---

## 📊 Visual Diagnostic Output & Graphics

* **Panel 1 (CPU Bandwidth Allocation Stacked Bar Chart):** Visualizes system resource partitioning: Gravitational Load $g_\Omega^2$ (red), Kinematic Load $v_\Omega^2$ (blue), Internal Tension $\alpha_\Omega^2$ (yellow), and Residual CPU Budget $u_\Omega$ (green)[cite: 38]. Displays an inline annotation box detailing exact budget percentages or $100\%$ CPU Kernel Lock status[cite: 38].
* **Panel 2 (Radial Clock Frequency Dephasing Profile):** Plots daily clock drift in $\mu\text{s/day}$ as a function of radial distance $r / R_{\text{surface}}$ relative to deep space unthrottled baseline ($0\,\mu\text{s/day}$) and highlights the evaluation surface[cite: 38].

---

## 🎯 KPI Metrics & Export Deliverables

| Metric Key | Format / Unit | Physical Meaning |
| :--- | :--- | :--- |
| `System Status` | Status Badge | Operational capacity percentage or `🔒 KERNEL LOCK (0.00% CPU)`[cite: 38]. |
| `Integer Node Load (N)` | Scientific (`nodes`) | Total compiled proton node count $N = M / m_p$[cite: 38]. |
| `Gravitational Load (g_Ω²)` | Normalized Ratio | Gravitational address deletion load $2GM / c^2 R$[cite: 38]. |
| `Kinematic Load (v_Ω²)` | Normalized Ratio | Kinematic data routing load $v^2 / c^2$[cite: 38]. |
| `Shift vs Deep Space (Δf/f₀)` | Scientific Ratio | Absolute fractional frequency shift relative to unthrottled deep space ($-2.5040 \times 10^{-10}$ for GPS)[cite: 38]. |
| `Shift vs Earth Ground (Δf/f₀)` | Scientific Ratio | Relative fractional frequency shift observed from Earth ground ($+4.4661 \times 10^{-10}$ for GPS)[cite: 38]. |
| `Relative Drift vs Earth` | `μs/day` | Daily clock drift relative to an Earth surface ground clock ($+38.59\,\mu\text{s/day}$ for GPS)[cite: 38]. |

* **CSV Export Schema (`cpu_bandwidth_throttling_profile.csv`):**
  * `Radius [m]`: Distance from center in meters[cite: 38].
  * `Altitude [km]`: Altitude above surface in kilometers[cite: 38].
  * `r / R_surface`: Radial distance ratio $r / R_{\text{surface}}$[cite: 38].
  * `g_Omega^2`: Gravitational load $g_\Omega^2$[cite: 38].
  * `v_Omega^2`: Kinematic load $v_\Omega^2$[cite: 38].
  * `External Load Omega`: Combined dynamic load $g_\Omega^2 + v_\Omega^2$[cite: 38].
  * `Observable u_clock`: Local clock update rate $u_{\text{clock}}$[cite: 38].
  * `Delta_f / f0 (vs Space)`: Fractional frequency shift $\Delta f / f_0$ relative to deep space.
  * `Drift vs Space [μs/day]`: Daily clock drift relative to deep space in $\mu\text{s/day}$.
  * `Delta_f / f0 (vs Earth)`: Fractional frequency shift $\Delta f / f_0$ relative to Earth ground observer.
  * `Drift vs Earth [μs/day]`: Daily clock drift relative to Earth surface ground clock in $\mu\text{s/day}$.

---

## ⚡ Performance Benchmarks & Hardware Behavioral Notes

* **Execution Benchmarks (Analytical 1D Solver):**
  * **Earth Ground Observer:** `< 0.2 ms`[cite: 38]
  * **GPS Constellation Satellite Orbit ($20,180\text{ km}$):** `< 0.05 ms`[cite: 38]
  * **ISS Orbit ($400\text{ km}$):** `< 0.05 ms`[cite: 38]
  * **Solar Photosphere ($1.988 \times 10^{30}\text{ kg}$):** `< 0.05 ms`[cite: 38]
  * **Sirius B White Dwarf ($1.018 M_\odot$):** `< 0.05 ms`[cite: 38]
  * **Crab Pulsar Neutron Star ($1.4 M_\odot$):** `< 0.05 ms`[cite: 38]
  * **Sagittarius A\* Horizon ($4.154 \times 10^6 M_\odot$):** `< 0.05 ms`[cite: 38]

* **Validation against Real-World Atomic Clock Measurements:**
  * **GPS Satellite Constellation Calibration:** GPS clocks run at $+38.59\,\mu\text{s/day}$ faster than Earth ground clocks ($\Delta f / f_0 = +4.4661 \times 10^{-10}$) because weak gravitational load at $20,180\text{ km}$ altitude frees up more bandwidth than speed throttling consumes[cite: 38].
  * **ISS Low Earth Orbit Calibration:** ISS clocks run at $-24.62\,\mu\text{s/day}$ slower than Earth ground clocks ($\Delta f / f_0 = -2.8495 \times 10^{-10}$) because high LEO orbital speed ($7.67\text{ km/s}$) consumes more bandwidth than the slight elevation gain frees up[cite: 38].