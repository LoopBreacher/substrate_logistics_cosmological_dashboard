# 🌌 Substrate Logistics: Cosmological Engine & Observatory

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://substrate-logistics.streamlit.app/)
[![GitHub](https://img.shields.io/badge/GitHub-Repository-181717?style=flat&logo=github)](https://github.com/LoopBreacher/substrate_logistics_cosmological_dashboard)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0003--8413--6027-A6CE39?style=flat&logo=orcid&logoColor=white)](https://orcid.org/0009-0003-8413-6027)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An interactive, first-principles cosmological simulation engine and observatory. **Substrate Logistics** replaces empirical curve-fitting parameters (Dark Energy $\Omega_\Lambda$, Cold Dark Matter halos, and FLRW single-scalar constraints) with a zero-parameter metric unspooling framework driven directly by 3D compiled matter density fields $\rho_N(\vec{r})$ across real astronomical catalogs (`UNG 2013`, `Cosmicflows-4`, `2M++`) via a **4-tier Multi-Band & Color-Dependent Stellar Mass Pipeline** (near-IR $K_s$, WISE W1, Bell+03 $B-V$, with legacy constant $M/L$ comparison toggle) and the zero-parameter **Two-Condition Compilation Gate**.

---

## 📖 Key Highlight: Resolving the Hubble Tension

Before exploring the individual modules, read our **[Comparative Architecture Guide](documentation/Comparative%20Architecture%20Guide.md)**. It provides a detailed theoretical comparison between $\Lambda\text{CDM}$ and Substrate Logistics, demonstrating why local $H_0 \approx 73\text{ km/s/Mpc}$ measurements and CMB $H_0 \approx 67.4\text{ km/s/Mpc}$ ($H_{\text{global}} = 67.4275\text{ km/s/Mpc}$) observations reflect two distinct, physically valid regimes rather than an observational crisis.

---

## ⚙️ Physical Invariants & Anchor Constants

| Parameter | Symbol | Value | Units / Description |
| :--- | :--- | :--- | :--- |
| **Global Vacuum Floor** | $H_{\text{global}}$ | `67.4275` | $\text{km/s/Mpc}$ (Derived closed-form, zero $G$) |
| **Substrate Anchor Constant** | $\mathcal{K}_\Omega$ | `1.1159160239e-37` | $\text{m}^3 / (\text{node}\cdot\text{s}^2)$ |
| **Universal Attenuation Tensor** | $\mathcal{T}_\Omega$ | `1.3014688631e+19` | dimensionless Tensor |
| **Topological Drag** | $\mathcal{C}_\delta$ | `4.72136e-4` | Exact algebraic $(\sqrt{5}-2)/500$ |
| **Mass per Integer Node** | $m_p$ | `1.6726219e-27` | $\text{kg}$ (CODATA 2022 proton mass) |
| **QCD Bag Floor ($B_{\text{QCD}}$)** | $`P_{\delta,\text{micro}}`$ | `47.624` | $\text{MeV/fm}^3$ ($7.6302 \times 10^{33}\text{ Pa}$ micro floor) |
| **Macro Drag Pressure** | $`P_{\delta,\text{macro}}`$ | `2.8339e-13` | $\text{Pa}$ ($0.2834\text{ pPa}$ cosmic floor) |
| **Structure Growth** | $S_8$ | `0.7542 / 0.8327` | Void / Global 120-cell partition |
| **Disk Truncation Cliff** | $`\frac{1}{2}\ln(1/\mathcal{C}_\delta)`$ | `3.8291` | $R_d$ (2D exponential disk cutoff) |
| **BAO Natural Ruler** | $r_s$ | `149.25` | $\text{Mpc}$ (First-principles acoustic scale) |
| **Bohr Phase-Lock Scale** | $`E_{\text{Ryd}} / T_{\text{lock}}`$ | `13.6057 eV / 105,258 K` | Zero-parameter EM unlock ($\alpha^{-1}$, $m_p/m_e$) |

---

## 🚀 Quickstart & Installation

```bash
# Clone the repository
git clone https://github.com/LoopBreacher/substrate_logistics_cosmological_dashboard.git
cd substrate_logistics_cosmological_dashboard

# Create and activate a virtual environment
python -m venv venv
source venv/bin/activate
# On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Launch the Streamlit dashboard
streamlit run app.py
```

---

## 🔬 Observatory Modules & Domain Breakdown

### Domain 1: Distances & $H_0$ Tension
Focuses on line-of-sight ray tracing, directional anisotropy, and zero-parameter supernova calibration.

| Module | Visual Diagnostic Preview | Physical Principle & Highlights | Detailed Docs |
| :--- | :--- | :--- | :--- |
| **Sky Anisotropy Maps** | <img src="documentation/plot_substrate_directional_correction_field.png" width="300" alt="Sky Anisotropy Maps" /> | Full-sky $H_0(\hat{n})$ unspooling maps and Euclid/LSST distance modulus anomaly forecasts $\Delta \mu(\hat{n})$. | [Read Documentation](documentation/readme_substrate_directional_correction_field.md) |
| **Line-of-Sight Kinematics** | <img src="documentation/plot_line_of_sight_kinematics.png" width="300" alt="Line-of-Sight Kinematics" /> | Decouples pure cosmological expansion ($cz_{\text{exp}}$) from peculiar velocity ($v_{\text{pec}}$) along TRGB/Cepheid target paths. | [Read Documentation](documentation/readme_line_of_sight_kinematics.md) |
| **Redshift Path Integral** | <img src="documentation/plot_redshift_path_integral.png" width="300" alt="Redshift Path Integral" /> | Integrates redshift accumulation through dense filament corridors versus empty void sightlines. | [Read Documentation](documentation/readme_redshift_path_integral.md) |
| **Supernova Fit ($\Omega_\Lambda = 0$)** | <img src="documentation/plot_substrate_hubble_diagram_fitter.png" width="300" alt="Supernova Fit" /> | Eliminates Dark Energy ($\Omega_\Lambda = 0$), reproducing Pantheon+ high-<i>z</i> dimming via photon memory drag attenuation. | [Read Documentation](documentation/readme_substrate_hubble_diagram_fitter.md) |
| **KBC Supervoid Profile** | <img src="documentation/plot_substrate_kbc_void_profile.png" width="300" alt="KBC Supervoid Profile" /> | Models radial expansion relaxation from local filament cores ($73\text{ km/s/Mpc}$) down to the global floor ($67.4275\text{ km/s/Mpc}$). | [Read Documentation](documentation/readme_substrate_kbc_void_profile.md) |

---

### Domain 2: Cosmic Web & Velocity Dynamics
Focuses on spatial velocity gradients, void wall dynamics, 2D RSD correlations, rotation curves, 3D web rendering, and 3D Dark Energy field theory.

| Module | Visual Diagnostic Preview | Physical Principle & Highlights | Detailed Docs |
| :--- | :--- | :--- | :--- |
| **Bulk Flow Field** | <img src="documentation/plot_substrate_bulk_flow_field.png" width="300" alt="Bulk Flow Field" /> | Maps 2D spatial vector fields driven by expansion gradients <b><i>v</i></b><sub>inflow</sub> = α · ∇<i>H</i><sub>local</sub>(<b>r</b>). | [Read Documentation](documentation/readme_substrate_bulk_flow_field.md) |
| **Void Boundary Dynamics** | <img src="documentation/plot_substrate_void_boundary_dynamics.png" width="300" alt="Void Boundary Dynamics" /> | Evaluates vacuum floor locking in void cores and velocity shear spikes ($dv_{\text{rec}}/dr$) across wall corridors. | [Read Documentation](documentation/readme_substrate_void_boundary_dynamics.md) |
| **Redshift Space Distortions** | <img src="documentation/plot_substrate_redshift_space_distortions.png" width="300" alt="Redshift Space Distortions" /> | Zero-parameter simultaneous extraction of Kaiser coherent squashing and cluster-core Fingers-of-God (FoG) elongation. | [Read Documentation](documentation/readme_substrate_redshift_space_distortions.md) |
| **Galactic Rotation Curves** | <img src="documentation/plot_substrate_galactic_rotation_curve.png" width="300" alt="Galactic Rotation Curves" /> | Eradicates Dark Matter halos via Freeman (1970) thin disk & spherical profiles, the 3.8291 <i>R</i><sub>d</sub> truncation cliff, and flat floor $`v_{\text{flat}} = \sqrt[4]{N_{\text{total}}\mathcal{K}_\Omega a_\Omega}`$. | [Read Documentation](documentation/readme_substrate_rotation_curves.md) |
| **3D Cosmic Web Viewer** | <img src="documentation/plot_substrate_3d_density_viewer.png" width="300" alt="3D Cosmic Web Viewer" /> | Interactive Plotly 3D Cartesian rendering of active survey galaxies and volumetric filament web isosurfaces. | [Read Documentation](documentation/readme_substrate_3d_density_viewer.md) |
| **3D Dark Energy Field** | <img src="documentation/plot_substrate_3d_dark_energy.png" width="300" alt="3D Dark Energy Field" /> | Derives the 3D spatial Dark Energy equation of state <i>w</i><sub>0</sub>(<b>r</b>) and pointer recycling drift <i>w</i><sub>a</sub>(<b>r</b>) from compiled matter density, enforcing the cosmic sum rule <i>w</i><sub>0</sub> + <i>w</i><sub>a</sub> = -1.0. | [Read Documentation](documentation/readme_substrate_3d_dark_energy.md) |

---

### Domain 3: Astrophysical & Horizon Tests
Focuses on bare-metal geometric rulers, hybrid CMB dipole shifts, and optical light deflection.

| Module | Visual Diagnostic Preview | Physical Principle & Highlights | Detailed Docs |
| :--- | :--- | :--- | :--- |
| **BAO Ruler Warping** | <img src="documentation/plot_substrate_bao_ruler_warping.png" width="300" alt="BAO Ruler Warping" /> | Alcock-Paczynski anisotropic distortion of the $149.25\text{ Mpc}$ acoustic sphere across local filament structures. | [Read Documentation](documentation/readme_substrate_bao_ruler_warping.md) |
| **CMB Dipole Anomaly** | <img src="documentation/plot_substrate_cmb_dipole_anomaly.png" width="300" alt="CMB Dipole Anomaly" /> | Hybrid Solar motion + line-of-sight redshifting vector sum (J2000 Galactic-to-ICRS aligned), shifting the observable dipole axis toward Centaurus. | [Read Documentation](documentation/readme_substrate_cmb_dipole_anomaly.md) |
| **Gravitational Lensing** | <img src="documentation/plot_substrate_galactic_lensing.png" width="300" alt="Gravitational Lensing" /> | Zero-dark-matter light deflection separating $1/b$ baryonic address deletion from the linear unspooling tax $4 a_\Omega b / c^2$ and $(1-b) = 0.8203$ cluster bias. | [Read Documentation](documentation/readme_substrate_galactic_lensing.md) |
| **$S_8$ Weak Lensing Tension** | <img src="documentation/plot_substrate_s8_tension.png" width="300" alt="S8 Weak Lensing Tension" /> | Resolves the weak lensing vs. Planck tension via the **Two-Condition Compilation Gate** ($`\Theta(\rho_N - 4)\Theta(\Lambda_{\text{lock}} - 1)`$): 3D catalog $S_8$ field ($`S_{8,\text{void}} = 0.7542`$, $`S_{8,\text{global}} = 0.8327`$, $`S_{8,\text{core}} = 0.8603`$, $(1-b)=0.8203$) plus 16-regime phase diagram & Cool-Core `r_lock = 5.42 kpc` ICM explorer. | [Read Documentation](documentation/readme_substrate_s8_tension.md) |

---

### Domain 4: Master Regulatory Protocol & Bandwidth Allocation
Focuses on spatial memory management, coordinate deletion rates, hardware clock throttling, and universal pressure boundaries.

| Module | Visual Diagnostic Preview | Physical Principle & Highlights | Detailed Docs |
| :--- | :--- | :--- | :--- |
| **Spatial Garbage Collection** | <img src="documentation/plot_spatial_garbage_collection_engine.png" width="300" alt="Spatial Garbage Collection" /> | Reclassifies gravity as $E=0$ spatial address deletion ($\dot{N}_{\text{pixels}}$), evaluating escape inflow, $R_{\text{quantum}}$, magnetospheric shells, and $1.0c$ Kernel Lock flushes. | [Read Documentation](documentation/readme_spatial_garbage_collection_engine.md) |
| **CPU Bandwidth Throttling** | <img src="documentation/plot_cpu_bandwidth_throttling.png" width="300" alt="CPU Bandwidth Throttling" /> | Explains gravitational time dilation as metric processing latency under the Pythagorean bandwidth invariant ($c^2 = v^2 + u^2$). | [Read Documentation](documentation/readme_cpu_bandwidth_throttling.md) |
| **Topological Drag Pressure Floor ($P_\delta$)** | <img src="documentation/plot_substrate_topological_drag_floor.png" width="300" alt="Topological Drag Pressure Floor" /> | 4-regime scale invariance spanning 38 orders of magnitude: MIT Bag / pion ($47.624\text{ MeV/fm}^3$, $138.31\text{ MeV}$), magnetospheres ($6.784R$), heliosphere ($89.44/120.74\text{ AU}$), $3.8291 R_d$ disk cliffs, and $R_{\text{GC}}$. | [Read Documentation](documentation/readme_substrate_topological_drag_floor.md) |

---

## 👤 Author & Citation

**Marco Lindenbeck**  
ORCID: [0009-0003-8413-6027](https://orcid.org/0009-0003-8413-6027)  
Repository: [LoopBreacher/substrate_logistics_cosmological_dashboard](https://github.com/LoopBreacher/substrate_logistics_cosmological_dashboard)

If you use this software or framework in your research, please cite:

```bibtex
@software{lindenbeck2026substrate,
  author    = {Lindenbeck, Marco},
  title     = {Substrate Logistics: Cosmological Engine \& Observatory},
  year      = {2026},
  publisher = {GitHub},
  url       = {https://github.com/LoopBreacher/substrate_logistics_cosmological_dashboard}
}
```

---

## Theoretical Foundation & Citations

This cosmological simulation engine and observatory is grounded in the **Substrate Logistics** theoretical manuscripts by **Marco Lindenbeck** ([ORCID: 0009-0003-8413-6027](https://orcid.org/0009-0003-8413-6027)):

1. **Lindenbeck, M. (2026).** *Substrate Logistics: The Master Regulatory Protocol.* Zenodo. [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.22127234-007ACC?style=flat)](https://doi.org/10.5281/zenodo.22127234)
2. **Lindenbeck, M. (2026).** *Substrate Logistics: The Garbage Collection Threshold.* Zenodo. [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.22774810-007ACC?style=flat)](https://doi.org/10.5281/zenodo.22774810)
3. **Lindenbeck, M. (2026).** *Substrate Logistics: The* $`S_8`$ *Weak Lensing Tension.* Zenodo. [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.23476628-007ACC?style=flat)](https://doi.org/10.5281/zenodo.23476628)
4. **Lindenbeck, M. (2026).** *Substrate Logistics: The Topological Drag Pressure Floor* ($`P_\delta`$). Zenodo. [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24042263-007ACC?style=flat)](https://doi.org/10.5281/zenodo.24042263)
5. **Lindenbeck, M. (2026).** *Substrate Logistics: A Synthesis of the Universal Attenuation Tensor ($\mathcal{T}_\Omega$) from First Principles.* Zenodo. [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.22738464-007ACC?style=flat)](https://doi.org/10.5281/zenodo.22738464)

### BibTeX
```bibtex
@paper{lindenbeck_2026_22127234,
  author       = {Lindenbeck, Marco},
  title        = {Substrate Logistics: The Master Regulatory Protocol},
  month        = aug,
  year         = 2026,
  publisher    = {Lindenbeck, Marco},
  doi          = {10.5281/zenodo.22127234},
  url          = {https://doi.org/10.5281/zenodo.22127234},
}

@paper{lindenbeck_2026_22774810,
  author       = {Lindenbeck, Marco},
  title        = {Substrate Logistics: The Garbage Collection Threshold},
  month        = sep,
  year         = 2026,
  publisher    = {Lindenbeck, Marco},
  doi          = {10.5281/zenodo.22774810},
  url          = {https://doi.org/10.5281/zenodo.22774810},
}

@paper{lindenbeck_2026_23476628,
  author       = {Lindenbeck, Marco},
  title        = {Substrate Logistics: The $S_8$ Weak Lensing Tension},
  month        = oct,
  year         = 2026,
  publisher    = {Lindenbeck, Marco},
  doi          = {10.5281/zenodo.23476628},
  url          = {https://doi.org/10.5281/zenodo.23476628},
}

@paper{lindenbeck_2026_24042263,
  author       = {Lindenbeck, Marco},
  title        = {Substrate Logistics: The Topological Drag Pressure Floor ($P_\delta$) — Scale-Invariant Boundary Mechanics Across 38 Orders of Magnitude: Deriving the QCD Bag Constant, the 138.3 MeV Pion Mass, the 120.7 AU Heliopause, and the 3.83 R_d Galactic Disk Truncation Law},
  month        = nov,
  year         = 2026,
  publisher    = {Lindenbeck, Marco},
  doi          = {10.5281/zenodo.24042263},
  url          = {https://doi.org/10.5281/zenodo.24042263},
}

@paper{lindenbeck_2026_22738464,
  author       = {Lindenbeck, Marco},
  title        = {Substrate Logistics: A Synthesis of the Universal Attenuation Tensor ($\mathcal{T}_\Omega$) from First Principles},
  month        = sep,
  year         = 2026,
  publisher    = {Lindenbeck, Marco},
  doi          = {10.5281/zenodo.22738464},
  url          = {https://doi.org/10.5281/zenodo.22738464},
}
```

---

## 📡 Multi-Band Stellar Mass Pipeline & Local Volume Audit

The shared catalog engine ([`modules/data_loader.py`](modules/data_loader.py)) removes the legacy single-band constant `M/L_B = 0.44` dial via a vectorized 4-tier photometric mass estimator (`compute_multiband_stellar_nodes`) persisted in `*_v2.npz` caches (with graceful offline fallback to `v1` caches and a live sidebar toggle between **Multi-Band / Color-Dependent** and **Legacy Constant M/L**):

1. **Tier 1 — Near-IR $K_s$-band (`M/L_K = 0.60`, $M_{K,\odot} = 3.27$):** McGaugh & Schombert (2014) near-IR mass-to-light calibration (full coverage in `UNG 2013` and `2M++` `Ksmag`).
2. **Tier 2 — Near-IR WISE W1 (`M/L_W1 = 0.50`, $M_{W1,\odot} = 3.24$ Vega):** Lelli et al. (2016) calibration; for `Cosmicflows-4`, cross-matched by `PGC` ID against Kourkchi et al. (2020, `J/ApJ/902/145`) with Jarrett et al. (2011) `W1_Vega = W1_AB - 2.699` conversion.
3. **Tier 3 — Optical $B - V$ Color (Bell et al. 2003, Table 7):** $\log_{10}(M/L_B) = -0.942 + 1.737(B-V) - 0.15$ (Chabrier/Kroupa IMF offset), clamped to `[0.10, 4.00]`.
4. **Tier 4 — Single-Band Legacy Fallback:** Falls back to single-band `default_m_l` when multi-band photometry is unavailable.

| Catalog | Mode | N | Photometry Provenance | $`\rho_N(\le 10\text{ Mpc})`$ [nodes/m³] | $`H_{\text{local}}(\le 10\text{ Mpc})`$ [km/s/Mpc] | Module 17 `<S8>` |
| :--- | :--- | ---: | :--- | ---: | ---: | ---: |
| UNG 2013 | `legacy` | 846 | Legacy Constant `M/L` (`v2`) | `0.00251` | `68.92` | `0.7796` |
| UNG 2013 | `multiband` | 846 | `100% K-band` (`v2`) | `0.01279` | `70.80` | `0.7890` |
| Cosmicflows-4 | `legacy` | 17,511 | Legacy Constant `M/L` (`v2`) | `0.00107` | `68.40` | `0.8317` |
| Cosmicflows-4 | `multiband` | 17,511 | `27% WISE W1 \| 73% Fallback` (`v2`) | `0.00228` | `68.85` | `0.8313` |

*Inside the `10 Mpc` Local Volume, switching `UNG 2013` from legacy `B`-band to 100% near-IR `K`-band lifts `H_local` from `68.92` to `70.80 km/s/Mpc` (`+1.00 km/s/Mpc` vs. CCHP TRGB `69.8 ± 1.9`, `-2.24 km/s/Mpc` vs. SH0ES `73.04 ± 1.04`).*

---

## 🧪 Verification Benchmark & CLI

Run the headless 18-module verification harness directly from the repository root using [`tools/benchmark_dashboard.py`](tools/benchmark_dashboard.py):

```bash
python tools/benchmark_dashboard.py
```

This benchmark validates:
1. Exact closed-form first-principles constants derived in [`modules/data_loader.py`](modules/data_loader.py) without CODATA $G$.
2. Orthogonal J2000 Galactic-to-ICRS sky coordinate alignment (`R_GAL_TO_ICRS`) verified against `astropy.coordinates.SkyCoord` (`< 1e-14` error).
3. Headless execution across all **18 computational modules** under `mass_mode="multiband"` (with real local catalog caches or deterministic fallback), verifying non-empty output tables, zero `NaN`/`Inf` numeric values, and clean Matplotlib renders alongside physics regression checks.
4. Zero-parameter **Two-Condition Compilation Gate** regression checks verifying `E_Ryd = 13.6057 eV`, `T_lock = 105,258.34 K`, Cool-Core Perseus `r_lock = 5.42 kpc`, and exact `8 PASSIVE / 8 COMPILED` verdicts across all 16 astrophysical regimes.
5. **4-Tier Photometric Mass Audit** unit tests and side-by-side `legacy` vs. `multiband` Local Volume `rho_N(<=10 Mpc)`, `H_local`, and `<S8>` comparisons.

---

## 🤝 Development Methodology & AI Collaboration Disclosure

* **Theoretical Physics & Equations:** All first-principles cosmological equations, 120-cell hardware partitions, algebraic invariants, and falsifiable astrophysical benchmarks in this repository originate from the **Substrate Logistics** research program by **Marco Lindenbeck**.
* **Software Engineering & AI Assistance:** Generative AI coding tools were utilized strictly as software engineering assistants for vectorized NumPy/`cKDTree` optimization, Matplotlib dark-theme visualization layout, Streamlit UI scaffolding, headless benchmark test wiring, and Markdown documentation formatting under direct human verification.

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.