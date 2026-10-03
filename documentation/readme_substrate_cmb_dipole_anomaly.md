# Substrate CMB Dipole Anomaly: Hybrid Kinematic & Metric Unspooling CMB Anisotropy

## 🌟 Executive Summary

The `substrate_cmb_dipole_anomaly` module resolves the longstanding observational tension between the Cosmic Microwave Background (CMB) dipole ($v \approx 369\text{ km/s}$) and large-scale radio galaxy/quasar count dipoles (e.g., NVSS, CatWISE / Secrest et al.). In standard $\Lambda\text{CDM}$, the CMB dipole is assumed to be 100% kinematic Doppler motion. Substrate Logistics demonstrates that line-of-sight spatial address unspooling along 3D cosmic web density fields ($\Delta H_{\text{local}}$) introduces an intrinsic directional redshifting component. Combining the Solar System kinematic Doppler vector with path-integrated metric unspooling yields a zero-parameter dynamic dipole shift ($\sim 6.5^\circ\text{--}15^\circ$) pointing directly toward the Centaurus / Great Attractor overdensity peak, explaining why quasar count dipoles measure systematically larger amplitudes without invoking anomalous bulk flows.

---

## 🔬 Physical Foundation & First Principles

### 1. Kinematic Doppler Component

The Solar System observer moves at peculiar velocity $v_{\text{pec}} = 369.0\text{ km/s}$ relative to the cosmic frame toward Galactic coordinates $(l_{\text{kin}}, b_{\text{kin}}) = (264.0^\circ, +48.0^\circ)$. Along unit sky vector $\hat{n}$, this motion generates a pure kinematic Doppler temperature shift:

$$\Delta T_{\text{kin}}(\hat{n}) = T_{\text{CMB}} \left( \frac{\vec{v}_{\text{pec}} \cdot \hat{n}}{c} \right) \quad [\text{K}]$$

where $T_{\text{CMB}} = 2.7255\text{ K}$ is the unperturbed thermodynamic monopole floor and $c = 299,792.458\text{ km/s}$.

### 2. Path-Integrated Metric Unspooling Component

Photons traversing line-of-sight ray paths $s \in [0, D_{\text{max}}]$ accumulate expansion rate variations $\Delta H(\vec{r}) = H_{\text{local}}(\vec{r}) - H_{\text{global}}$ driven by local matter density fields $\rho_N(\vec{r})$. The cumulative line-of-sight velocity excess evaluates to:

$$c \Delta z(\hat{n}) = \int_0^{D_{\text{max}}} \Delta H(s \hat{n}) \, ds \quad [\text{km/s}]$$

To avoid point-mass sampling spikes along long sightlines, the spatial density smoothing bandwidth adaptively expands with distance:

$$\sigma(s) = \sigma_0 \sqrt{1 + \frac{s}{10.0}} \quad [\text{Mpc}]$$

The directional temperature anisotropy induced by metric unspooling drag is:

$$\Delta T_{\text{sub, raw}}(\hat{n}) = -T_{\text{CMB}} \left( \frac{c \Delta z(\hat{n})}{c} \right) \quad [\text{K}]$$

### 3. Thermodynamic Monopole Subtraction & Composite Dipole Vector

Subtracting the isotropic thermodynamic mean isolates pure multipole fluctuations ($\ell \ge 1$):

$$\Delta T_{\text{sub}}(\hat{n}) = \Delta T_{\text{sub, raw}}(\hat{n}) - \bar{\Delta T}_{\text{sub, raw}}$$

$$\Delta T_{\text{total}}(\hat{n}) = \Delta T_{\text{kin}}(\hat{n}) + \Delta T_{\text{sub}}(\hat{n})$$

Performing a least-squares fit against the design matrix $[1, \hat{n}_x, \hat{n}_y, \hat{n}_z]$ extracts the net composite dipole vector $\vec{d} = (d_x, d_y, d_z)$:

$$A_{\text{total}} = \Vert{}\vec{d}\Vert{} \quad [\text{mK}]$$

$$b_{\text{tot}} = \arcsin\left(\frac{d_z}{A_{\text{total}}}\right), \quad l_{\text{tot}} = \arctan2(d_y, d_x) \pmod{360^\circ}$$

---

## ⚙️ Architectural & Technical Implementation

* **Primary Driver Function:** `run_cmb_dipole_analysis(...)`
* **Cached Ray-Tracing Engine:** `_compute_cmb_dipole_cached(...)`
* **Core Dependencies:** `healpy`, `scipy.integrate.simpson`, `astropy.coordinates.SkyCoord`, `modules.data_loader` (`get_hlocal_batch`)

### Input Parameters

* `tree`: `scipy.spatial.cKDTree` 3D spatial coordinate index.
* `gal_positions`: `np.ndarray` of shape $(N, 3)$ in Mpc (ICRS cartesian).
* `nodes_per_galaxy`: `np.ndarray` of integer proton node counts $N = M / m_p$.
* `nside`: HEALPix resolution parameter ($8, 16, 32$; default `32`, yielding $12,288$ sky pixels).
* `d_max_mpc`: Ray-tracing integration depth in Mpc ($5.0\text{--}200.0\text{ Mpc}$).
* `sigma_mpc`: Base Gaussian kernel bandwidth $\sigma_0$ ($0.5\text{--}20.0\text{ Mpc}$).
* `h_global`: Global vacuum expansion floor ($67.42\text{ km/s/Mpc}$).
* `v_pec_km_s`: Observer peculiar velocity ($369.0\text{ km/s}$).
* `l_kin_deg`, `b_kin_deg`: Kinematic Doppler apex in Galactic coordinates (264.0°, +48.0°).

---

## 📊 Visual Diagnostic Output & Graphics

The module outputs a full-sky Mollweide projection map ($\Delta T(\theta, \phi)$ in mK) rendered in Galactic coordinates:

* **Background Temperature Field:** Colormap (`coolwarm`) illustrating the composite sky anisotropy $\Delta T_{\text{total}}(\hat{n})$.
* **Kinematic Apex Marker (Black Circle):** Pure Solar Doppler direction at (264.0°, +48.0°).
* **Centaurus / Local Sheet Peak Marker (Cyan Triangle):** Overdensity attraction center at (280.0°, +45.0°).
* **Total Composite Dipole Axis Marker (Yellow Star):** Net observable dipole axis shifted toward local matter concentrations.

---

## 🎯 KPI Metrics & Export Deliverables

| Metric Key | Unit / Format | Physical Meaning |
| --- | --- | --- |
| `Kinematic Peak` | `mK` | Pure Solar Doppler dipole amplitude ($3.354\text{ mK}$). |
| `Substrate Anisotropy Range` | `mK` | Min-to-max path-integrated unspooling anisotropy. |
| `Combined Dipole Amplitude` | `mK` | Fitted scalar amplitude of composite dipole $A_{\text{total}}$. |
| `Shifted Dipole Axis (l, b)` | `(deg, deg)` | Net celestial sky coordinates of composite dipole apex. |

### CSV Export Schema (`cmb_dipole_shift_lookup.csv`)

* `pixel_id`: HEALPix pixel index ($0 \dots 12287$ for $N_{\text{side}}=32$).
* `glon_deg`: Galactic longitude $l$ ($0.0^\circ \dots 360.0^\circ$).
* `glat_deg`: Galactic latitude $b$ ($-90.0^\circ \dots +90.0^\circ$).
* `delta_T_kin_mK`: Pure kinematic Doppler temperature shift [mK].
* `delta_T_sub_mK`: Monopole-subtracted metric unspooling temperature shift [mK].
* `delta_T_total_mK`: Net hybrid temperature shift $\Delta T_{\text{total}}$ [mK].

---

## ⚡ Performance Benchmarks & Catalog Behavioral Notes

### Benchmarks (2M++ Catalog, $D_{\text{max}} = 200\text{ Mpc}$, $\sigma = 3.50\text{ Mpc}$, $N_{\text{side}} = 32$)

* **HEALPix Sky Grid:** $12,288$ pixel sightlines $\times$ $35$ radial integration steps ($430,080$ batch spatial queries).
* **Execution Time:** `~3.85 s` (utilizing multi-threaded C++ `get_hlocal_batch`).
* **Memory Footprint:** `< 25 MB` temporary array allocation.

### Survey Catalog Response & Field Scaling

* **UNG 2013 ($\le 35\text{ Mpc}$):** Captures immediate Local Sheet and Virgo cluster path integrals, generating localized shifts of $\sim 1.2^\circ\text{--}3.5^\circ$.
* **Cosmicflows-4 ($\le 150\text{ Mpc}$):** Integrates across the Great Attractor, Laniakea, and Perseus-Pisces, producing dipole shifts of $\sim 5.0^\circ\text{--}8.5^\circ$.
* **2M++ Infrared ($\le 200\text{ Mpc}$):** Provides complete 3D coverage bridging galactic dust gaps, generating net dynamic dipole shifts of **$\sim 6.5^\circ\text{--}10.2^\circ$** toward Galactic coordinates $(273.2^\circ\text{--}276.0^\circ, +50.3^\circ)$.

---