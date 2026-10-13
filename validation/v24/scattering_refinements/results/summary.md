### Radius-integration and Legendre convergence (convergence.py)

| Case (model, um, r_e) | nodes V24 / reference / V23 | V24 vs 4x finer: ext, g, phase 5-150, 150-180 | V23 grid vs 4x finer: ext, g, phase 150-180 | 3.5 vs 5 r_e: ext | L (2 Nq) | max dchi at 2 Nq | reconstruction | noise |
|---|---|---|---|---|---|---|---|---|
| liquid 240 0.469 5 | 9377 / 37505 / 3349 | 1.4e-05, 1.2e-05, 2.2e-03, 6.9e-03 | 1.4e-03, 1.1e-03, 1.2e-01 | 5.9e-07 | 487 (487) | 3e-13 | 3e-09 | 2e-10 |
| liquid 240 0.646 21 | 28593 / 114369 / 2431 | 3.0e-06, 2.7e-06, 1.4e-03, 1.3e-03 | 1.9e-05, 1.7e-05, 1.4e-01 | 6.2e-07 | 1509 (1449) | 3e-12 | 4e-09 | 5e-09 |
| liquid 240 0.646 39 | 53073 / 212289 / 2431 | 6.5e-07, 5.4e-07, 1.2e-03, 6.2e-04 | 2.7e-04, 6.1e-06, 5.3e-02 | 6.5e-07 | 2751 (2751) | 6e-12 | 9e-09 | 1e-08 |
| liquid 240 0.646 5 | 6817 / 27265 / 2431 | 1.7e-05, 1.4e-05, 2.4e-03, 2.9e-03 | 2.4e-03, 1.8e-03, 1.9e-01 | 5.8e-07 | 359 (359) | 2e-13 | 9e-10 | 9e-11 |
| liquid 240 0.858 11 | 11281 / 45121 / 1830 | 1.2e-05, 9.7e-06, 1.7e-03, 2.3e-03 | 7.6e-04, 6.3e-04, 1.4e-01 | 5.5e-07 | 581 (581) | 2e-13 | 3e-09 | 3e-10 |
| liquid 240 11.03 1 | 113 / 449 / 200 | 0.0e+00, 2.2e-08, 7.6e-08, 4.5e-08 | 8.9e-03, 2.7e-04, 8.1e-04 | 2.1e-06 | 10 (10) | 2e-16 | 4e-11 | 1e-14 |
| liquid 240 11.03 21 | 2353 / 9409 / 200 | 0.0e+00, 0.0e+00, 5.9e-12, 2.7e-12 | 5.7e-07, 0.0e+00, 1.0e-07 | 5.7e-07 | 92 (92) | 2e-14 | 6e-08 | 1e-12 |
| liquid 240 2.13 21 | 8657 / 34625 / 737 | 1.2e-07, 0.0e+00, 8.7e-05, 2.5e-05 | 4.4e-05, 1.0e-04, 7.9e-02 | 6.0e-07 | 447 (447) | 3e-13 | 1e-09 | 1e-10 |
| liquid 240 3.79 1 | 241 / 961 / 414 | 0.0e+00, 0.0e+00, 4.7e-09, 3.8e-10 | 1.7e-04, 2.9e-05, 2.2e-04 | 1.4e-06 | 17 (17) | 4e-16 | 3e-09 | 6e-14 |
| liquid 240 3.79 21 | 4865 / 19457 / 414 | 0.0e+00, 0.0e+00, 1.1e-07, 5.8e-08 | 2.8e-04, 2.2e-04, 6.2e-02 | 5.8e-07 | 255 (255) | 3e-14 | 1e-09 | 2e-11 |
| ice sph 0.646 41 | 19441 / 77761 / 2431 | 1.8e-06, 1.7e-06, 3.2e-03, 3.0e-03 | 6.9e-06, 6.3e-06, 9.6e-02 | - | 2031 (2031) | 6e-12 | 2e-08 | 2e-08 |
| ice sph 0.646 9 | 19441 / 77761 / 2431 | 3.2e-05, 2.7e-05, 4.2e-03, 2.9e-03 | 1.4e-04, 1.1e-04, 6.6e-02 | - | 2031 (839) | 5e-13 | 7e-10 | 2e-09 |
| ice sph 11.03 93 | 1593 / 6369 / 200 | 7.9e-08, 0.0e+00, 8.8e-07, 6.7e-10 | 6.4e-06, 0.0e+00, 4.8e-08 | - | 134 (134) | 1e-13 | 2e-08 | 2e-11 |
| ice sph 2.13 41 | 5889 / 23553 / 737 | 6.4e-08, 0.0e+00, 9.3e-05, 3.6e-05 | 3.8e-07, 6.0e-08, 5.8e-02 | - | 615 (615) | 5e-13 | 4e-09 | 7e-10 |

### DISORT with NMOM = max(1000, L) (disort_nmom.py; original Intel kernel, NSTR 60)

| Case (um, r_e) | L | max(1000,L) = L bitwise | L+512 padded bitwise | max(L,61) bitwise (L<1000) | fluxes bitwise under truncation | max dR at N = 1000 / 1500 / 2000 / 2500 | s/call at 1000 / production |
|---|---|---|---|---|---|---|---|
| 0.469 21 | 2057 | True | True | n/a | True | 1.0e-01 / 5.0e-04 / 0.0e+00 / - | 0.087 / 0.100 |
| 0.646 21 | 1509 | True | True | n/a | True | 1.4e-03 / 0.0e+00 / - / - | 0.086 / 0.092 |
| 0.646 39 | 2751 | True | True | n/a | True | 1.4e+00 / 4.7e-02 / 8.3e-04 / 3.6e-06 | 0.089 / 0.108 |
| 0.858 11 | 581 | True | True | True | True | - / - / - / - | 0.090 / 0.095 |
| 3.79 21 | 255 | True | True | True | True | - / - / - / - | 0.100 / 0.097 |

### Class scattering properties, V24 - V23 (compare_optics.py)

| Model | max abs rel. extinction | max abs rel. extinction ratio | max abs dw | max abs dg | L range | max abs dchi (l < 1000) | consistency |
|---|---|---|---|---|---|---|---|
| liquid-water_240 | 9.1e-03 | 1.1e-02 | 8.3e-04 | 3.7e-03 | 9-3789 | 3.4e-02 | w in [0,1] True, chi_0 = 1 True, max abs chi 1, abs(chi_1 - g) <= 0e+00 |
| liquid-water_old | 9.0e-03 | 9.8e-03 | 8.4e-04 | 2.7e-03 | 9-3789 | 3.4e-02 | w in [0,1] True, chi_0 = 1 True, max abs chi 1, abs(chi_1 - g) <= 0e+00 |
| water-ice_sph | 2.3e-03 | 2.1e-03 | 3.5e-04 | 7.7e-04 | 26-2793 | 1.9e-01 | w in [0,1] True, chi_0 = 1 True, max abs chi 1, abs(chi_1 - g) <= 0e+00 |

### Compact LUTs, V24 - V23 (compare_compact.py; identical Grid B nodes and settings)

**liquid-water_240**: variables True, dimensions True, coordinates True, V24 finite True, bitwise identical variables 38, isolated R_0v call outliers 2

| Quantity | median | p99 | max | location of max | max rel. (V23 >= 0.01) | max by r_e (1, 5, 11, 21, 39 / ice 5, 41, 93) |
|---|---|---|---|---|---|---|
| extinction_coefficient | 1.2e-01 | 3.8e+00 | 5.0e+00 | {'effective_radius': 11.0, 'channels': 20} | 9.1e-03 | 2e-02, 6e-01, 5e+00, 5e-01, 2e+00 |
| extinction_coefficient_ratio | 3.1e-04 | 6.7e-03 | 9.1e-03 | {'effective_radius': 11.0, 'channels': 20} | 1.1e-02 | 2e-03, 3e-03, 9e-03, 3e-04, 3e-04 |
| single_scatter_albedo | 7.2e-07 | 6.5e-04 | 8.3e-04 | {'effective_radius': 21.0, 'channels': 20} | 1.1e-02 | 4e-04, 3e-04, 3e-04, 8e-04, 2e-04 |
| asymmetry_parameter | 4.6e-05 | 3.4e-03 | 3.7e-03 | {'effective_radius': 11.0, 'channels': 20} | 4.6e-03 | 2e-03, 3e-03, 4e-03, 5e-04, 3e-04 |
| R_0v | 1.2e-04 | 1.6e-02 | 1.1e-01 | {'relative_azimuth': 0.0, 'satellite_zenith': 30.0, 'solar_zenith': 30.0, 'optical_depth': 256.0, 'effective_radius': 39.0, 'solar_channels': 3} | 4.1e+00 | 6e-03, 3e-02, 9e-02, 1e-02, 1e-01 |
| R_0d | 1.6e-04 | 3.9e-03 | 6.0e-03 | {'solar_zenith': 0.0, 'optical_depth': 8.0, 'effective_radius': 5.0, 'solar_channels': 2} | 3.1e-02 | 2e-03, 6e-03, 5e-03, 8e-04, 3e-03 |
| T_0d | 7.8e-05 | 3.7e-03 | 6.2e-03 | {'solar_zenith': 0.0, 'optical_depth': 8.0, 'effective_radius': 11.0, 'solar_channels': 20} | 7.5e-02 | 2e-03, 6e-03, 6e-03, 2e-03, 3e-03 |
| T_00 | 3.5e-09 | 9.9e-04 | 3.0e-03 | {'solar_zenith': 0.0, 'optical_depth': 1.0, 'effective_radius': 11.0, 'solar_channels': 20} | 3.4e-02 | 7e-04, 9e-04, 3e-03, 9e-05, 1e-04 |
| R_dv | 6.0e-05 | 3.7e-03 | 6.0e-03 | {'satellite_zenith': 0.0, 'optical_depth': 8.0, 'effective_radius': 5.0, 'channels': 2} | 3.2e-02 | 2e-03, 6e-03, 5e-03, 8e-04, 3e-03 |
| T_dv | 3.9e-05 | 3.7e-03 | 6.2e-03 | {'satellite_zenith': 0.0, 'optical_depth': 8.0, 'effective_radius': 11.0, 'channels': 20} | 7.5e-02 | 4e-03, 6e-03, 6e-03, 2e-03, 3e-03 |
| R_dd | 5.9e-05 | 3.4e-03 | 4.9e-03 | {'optical_depth': 8.0, 'effective_radius': 5.0, 'channels': 2} | 1.7e-02 | 2e-03, 5e-03, 4e-03, 5e-04, 3e-03 |
| T_dd | 3.8e-05 | 3.9e-03 | 4.9e-03 | {'optical_depth': 8.0, 'effective_radius': 5.0, 'channels': 2} | 7.4e-02 | 4e-03, 5e-03, 4e-03, 1e-03, 3e-03 |
| E_md | 2.1e-05 | 3.7e-03 | 3.9e-03 | {'satellite_zenith': 0.0, 'optical_depth': 8.0, 'effective_radius': 1.0, 'thermal_channels': 32} | 1.4e-02 | 4e-03, 3e-03, 4e-03, 2e-03, 3e-04 |

**liquid-water_old**: variables True, dimensions True, coordinates True, V24 finite True, bitwise identical variables 38, isolated R_0v call outliers 0

| Quantity | median | p99 | max | location of max | max rel. (V23 >= 0.01) | max by r_e (1, 5, 11, 21, 39 / ice 5, 41, 93) |
|---|---|---|---|---|---|---|
| extinction_coefficient | 4.6e-02 | 3.4e+00 | 4.0e+00 | {'effective_radius': 11.0, 'channels': 20} | 9.0e-03 | 1e-02, 3e-01, 4e+00, 2e+00, 2e+00 |
| extinction_coefficient_ratio | 3.9e-04 | 6.0e-03 | 7.1e-03 | {'effective_radius': 11.0, 'channels': 20} | 9.8e-03 | 1e-03, 4e-03, 7e-03, 8e-04, 1e-04 |
| single_scatter_albedo | 1.4e-06 | 7.3e-04 | 8.4e-04 | {'effective_radius': 11.0, 'channels': 20} | 1.1e-02 | 4e-04, 2e-04, 8e-04, 3e-04, 2e-05 |
| asymmetry_parameter | 4.4e-05 | 2.4e-03 | 2.7e-03 | {'effective_radius': 11.0, 'channels': 6} | 3.2e-03 | 6e-04, 2e-03, 3e-03, 5e-04, 9e-05 |
| R_0v | 1.1e-04 | 1.2e-02 | 7.3e-02 | {'relative_azimuth': 0.0, 'satellite_zenith': 30.0, 'solar_zenith': 30.0, 'optical_depth': 256.0, 'effective_radius': 39.0, 'solar_channels': 3} | 3.6e-01 | 1e-02, 2e-02, 3e-02, 4e-03, 7e-02 |
| R_0d | 1.1e-04 | 7.1e-03 | 9.2e-03 | {'solar_zenith': 0.0, 'optical_depth': 256.0, 'effective_radius': 11.0, 'solar_channels': 6} | 2.7e-02 | 9e-04, 7e-03, 9e-03, 2e-03, 3e-03 |
| T_0d | 7.3e-05 | 3.7e-03 | 1.1e-02 | {'solar_zenith': 0.0, 'optical_depth': 8.0, 'effective_radius': 11.0, 'solar_channels': 6} | 2.0e-01 | 9e-04, 4e-03, 1e-02, 2e-03, 3e-03 |
| T_00 | 6.4e-09 | 1.4e-03 | 2.3e-03 | {'solar_zenith': 0.0, 'optical_depth': 1.0, 'effective_radius': 11.0, 'solar_channels': 20} | 2.8e-02 | 4e-04, 2e-03, 2e-03, 3e-04, 5e-05 |
| R_dv | 3.6e-05 | 6.9e-03 | 9.2e-03 | {'satellite_zenith': 0.0, 'optical_depth': 256.0, 'effective_radius': 11.0, 'channels': 6} | 2.8e-02 | 9e-04, 7e-03, 9e-03, 2e-03, 3e-03 |
| T_dv | 4.5e-05 | 3.6e-03 | 1.1e-02 | {'satellite_zenith': 0.0, 'optical_depth': 8.0, 'effective_radius': 11.0, 'channels': 6} | 2.0e-01 | 4e-03, 4e-03, 1e-02, 2e-03, 3e-03 |
| R_dd | 4.8e-05 | 5.9e-03 | 8.5e-03 | {'optical_depth': 256.0, 'effective_radius': 11.0, 'channels': 6} | 1.4e-02 | 7e-04, 6e-03, 9e-03, 2e-03, 3e-03 |
| T_dd | 5.1e-05 | 3.2e-03 | 9.0e-03 | {'optical_depth': 8.0, 'effective_radius': 11.0, 'channels': 6} | 2.0e-01 | 3e-03, 3e-03, 9e-03, 2e-03, 3e-03 |
| E_md | 2.2e-05 | 3.5e-03 | 3.6e-03 | {'satellite_zenith': 30.0, 'optical_depth': 8.0, 'effective_radius': 1.0, 'thermal_channels': 32} | 1.4e-02 | 4e-03, 8e-04, 4e-03, 4e-04, 5e-05 |

**water-ice_agg**: variables True, dimensions True, coordinates True, V24 finite True, bitwise identical variables 52, isolated R_0v call outliers 0

| Quantity | median | p99 | max | location of max | max rel. (V23 >= 0.01) | max by r_e (1, 5, 11, 21, 39 / ice 5, 41, 93) |
|---|---|---|---|---|---|---|
| extinction_coefficient | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'effective_radius': 5.0, 'channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| extinction_coefficient_ratio | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'effective_radius': 5.0, 'channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| single_scatter_albedo | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'effective_radius': 5.0, 'channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| asymmetry_parameter | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'effective_radius': 5.0, 'channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| R_0v | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'relative_azimuth': 0.0, 'satellite_zenith': 0.0, 'solar_zenith': 0.0, 'optical_depth': 0.125, 'effective_radius': 5.0, 'solar_channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| R_0d | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'solar_zenith': 0.0, 'optical_depth': 0.125, 'effective_radius': 5.0, 'solar_channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| T_0d | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'solar_zenith': 0.0, 'optical_depth': 0.125, 'effective_radius': 5.0, 'solar_channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| T_00 | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'solar_zenith': 0.0, 'optical_depth': 0.125, 'effective_radius': 5.0, 'solar_channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| R_dv | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'satellite_zenith': 0.0, 'optical_depth': 0.125, 'effective_radius': 5.0, 'channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| T_dv | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'satellite_zenith': 0.0, 'optical_depth': 0.125, 'effective_radius': 5.0, 'channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| R_dd | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'optical_depth': 0.125, 'effective_radius': 5.0, 'channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| T_dd | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'optical_depth': 0.125, 'effective_radius': 5.0, 'channels': 1} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |
| E_md | 0.0e+00 | 0.0e+00 | 0.0e+00 | {'satellite_zenith': 0.0, 'optical_depth': 0.125, 'effective_radius': 5.0, 'thermal_channels': 20} | 0.0e+00 | 0e+00, 0e+00, 0e+00 |

### Mie-step cost, 24 Grid B liquid radii (timing.py)

| um | V23 | + 3.5 r_e | + 3.5 r_e + refinement | + 3.5 r_e + adaptive Legendre | V24 (all three) | V24 / V23 |
|---|---|---|---|---|---|---|
| 0.469 | 241 s | 125 s | 2009 s | 369 s | 5990 s | 24.9 |
| 0.646 | 131 s | 69 s | 1089 s | 156 s | 2512 s | 19.2 |
| 1.64 | 24 s | 13 s | 213 s | 18 s | 284 s | 11.7 |
| 3.79 | 5 s | 3 s | 46 s | 3 s | 43 s | 8.3 |
| 11.03 | 1 s | 1 s | 11 s | 1 s | 8 s | 6.5 |
