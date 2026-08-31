# Spinach code index: examples

- Source root: `/home/kuprov/.openclaw/workspace/Spinach`
- Source commit: `f053e432a61d7144f3946d73d0a672e3ccfc3fc5`
- Source tree state: `clean`
- Path set: tracked MATLAB files from `git ls-files '*.m'`; untracked MATLAB files are excluded.
- Files indexed: **756** MATLAB files
- Generated: 2026-08-30T02:52:31

| File | Signature | Summary | LOC |
|---|---|---|---:|
| `examples/benchmarks/comm_gpu.m` | `comm_gpu(n)` | GPU communications benchmark. Adapted from example code in Matlab documentation. Data for IK's favourite NVIDIA cards on | 97 |
| `examples/benchmarks/fft_gpu.m` | `fft_gpu()` | GPU arithmetic benchmark -3D Fourier transforms. | 54 |
| `examples/benchmarks/iserstep_bench_hiord.m` | `iserstep_bench_hiord()` | Benchmarks iserstep higher-order methods on a chirped-frequency oscillator with radiation damping, that has a state-depe | 191 |
| `examples/benchmarks/mult_gpu.m` | `mult_gpu(precision)` | CPU and GPU matrix arithmetic benchmark. Set the argument to either 'single' or 'double'. IK's workstation output: GPU,  | 69 |
| `examples/benchmarks/parallelization_1.m` | `parallelization_1()` | Parallelization test: multi-threaded evaluation of observables in Hilbert space time propagation. For further informatio | 67 |
| `examples/benchmarks/parallelization_2.m` | `parallelization_2()` | Parallelization test: multi-threaded evaluation of observables in Hilbert space time propagation for pyrene radical spin | 47 |
| `examples/benchmarks/polyadic_bench.m` | `polyadic_bench()` | A benchmark for the polyadic object. | 100 |
| `examples/dnp_liq/ccdnp/freq_scan_main_text.m` | `freq_scan_main_text()` | Steady state nuclear magnetisation as a function of microwave frequency offset and the magnet field in a DNP experiment  | 101 |
| `examples/dnp_liq/ccdnp/freq_scan_si_sys_a.m` | `freq_scan_si_sys_a()` | Steady state nuclear magnetisation as a function of microwave frequency offset and the magnet field in a DNP experiment  | 105 |
| `examples/dnp_liq/ccdnp/freq_scan_si_sys_b.m` | `freq_scan_si_sys_b()` | Steady state nuclear magnetisation as a function of microwave frequency offset and the magnet field in a DNP experiment  | 105 |
| `examples/dnp_liq/ccdnp/rates_main_text.m` | `rates_main_text()` | Self-and cross-relaxation rates in cross-correlated DNP, considering a system with two electrons connected by exchange c | 99 |
| `examples/dnp_liq/ccdnp/rates_si_sys_a.m` | `rates_si_sys_a()` | Self-and cross-relaxation rates in cross-correlated DNP, considering a system with two electrons connected by exchange c | 99 |
| `examples/dnp_liq/ccdnp/rates_si_sys_b.m` | `rates_si_sys_b()` | Self-and cross-relaxation rates in cross-correlated DNP, considering a system with two electrons connected by exchange c | 100 |
| `examples/dnp_liq/ccdnp/states_vs_tauc_main_text.m` | `states_vs_tauc_main_text()` | Steady state populations of various spin states a function of rotational correlation time in a DNP experiment with two e | 143 |
| `examples/dnp_liq/ccdnp/states_vs_tauc_si_sys_a.m` | `states_vs_tauc_si_sys_a()` | Steady state populations of various spin states a function of rotational correlation time in a DNP experiment with two e | 143 |
| `examples/dnp_liq/ccdnp/states_vs_tauc_si_sys_b.m` | `states_vs_tauc_si_sys_b()` | Steady state populations of various spin states a function of rotational correlation time in a DNP experiment with two e | 143 |
| `examples/dnp_liq/ccdnp/tau_scan_main_text.m` | `tau_scan_main_text()` | Steady state nuclear magnetisation as a function of microwave frequency offset and the rotational correlation time in a  | 107 |
| `examples/dnp_liq/ccdnp/tau_scan_si_sys_a.m` | `tau_scan_si_sys_a()` | Steady state nuclear magnetisation as a function of microwave frequency offset and the rotational correlation time in a  | 104 |
| `examples/dnp_liq/ccdnp/tau_scan_si_sys_b.m` | `tau_scan_si_sys_b()` | Steady state nuclear magnetisation as a function of microwave frequency offset and the rotational correlation time in a  | 104 |
| `examples/dnp_liq/jdnp/energy_levels.m` | `energy_levels()` | Energy level diagram transition from the Zeeman limit to the exchange coupling limit in a two-electron system. | 52 |
| `examples/dnp_liq/jdnp/fig_1_exch_and_field_scan.m` | `fig_1_exch_and_field_scan()` | Matching condition plot for JDNP -proton polarisation at a particular time as a function of the external field and the i | 96 |
| `examples/dnp_liq/jdnp/fig_2_tau_and_field_traject.m` | `fig_2_tau_and_field_traject()` | Time evolution plot for JDNP: proton polarisation as a function of time for specific external fields and rotational corr | 104 |
| `examples/dnp_liq/jdnp/fig_3_time_dep_bot_row.m` | `fig_3_time_dep_bot_row()` | Time evolution plot for JDNP: proton polarisation as a function of time for specific external fields. The inter-electron | 86 |
| `examples/dnp_liq/jdnp/fig_3_time_dep_top_row.m` | `fig_3_time_dep_top_row()` | A demonstration that the JDNP effect vanishes when the second electron is removed from the system. Proton polarisation a | 87 |
| `examples/dnp_liq/jdnp/fig_4_state_amplitudes.m` | `fig_4_state_amplitudes()` | Time evolution of the individual states in the basis set built from the singlet-triplet basis on the two electrons and C | 126 |
| `examples/dnp_liq/jdnp/fig_5_spatial_distribution.m` | `fig_5_spatial_distribution()` | An illustration of the fact that JDNP effect does not vanish on position and orientation averaging in liquid phase. The  | 156 |
| `examples/dnp_liq/jdnp/fig_6_microwave_free.m` | `fig_6_microwave_free()` | A demonstration of Maria Grazia Concilio's microwave-free JDNP effect where a field ramp in combination with unequal rel | 134 |
| `examples/dnp_liq/jdnp/system_specification.m` | `[sys,inter,bas,parameters]=system_specification()` | Parameters of the 2e1n system used for the simulations reported in https://doi.org/10.1039/d1cp04186j | 50 |
| `examples/dnp_liq/odnp_liquid_1.m` | `odnp_liquid_1()` | Overhauser type DNP in liquid phase at room temperature, using a continu- ous on-resonance CW irradiation of the electro | 69 |
| `examples/dnp_liq/odnp_liquid_2.m` | `odnp_liquid_2()` | Overhauser type DNP in liquid phase at room temperature, after a perfect inversion pulse on the electron ESR signal. The | 80 |
| `examples/dnp_liq/odnp_liquid_3.m` | `odnp_liquid_3()` | Steady state nuclear magnetisation as a function of microwave frequency offset and the magnet field in a DNP experiment  | 95 |
| `examples/dnp_liq/sdnp/scalar_dnp.m` | `scalar_dnp()` | Field dependence of the couping factor between 13C of CHCl3 and the electron spin of a nitroxide radical. Further partic | 116 |
| `examples/dnp_mas/cross_effect_mas_dynam.m` | `cross_effect_mas_dynam()` | A MAS DNP simulation performed as described in Fred Mentink- Vigier's paper (Spinach rotation conventions are different) | 94 |
| `examples/dnp_mas/cross_effect_mas_enlev.m` | `cross_effect_mas_enlev()` | A MAS DNP simulation performed as described in Fred Mentink- Vigier's paper (Spinach rotation conventions are different) | 72 |
| `examples/dnp_mas/cross_effect_mas_powder.m` | `cross_effect_mas_powder()` | A MAS DNP simulation performed as described in Fred Mentink- Vigier's paper (Spinach rotation conventions are different) | 69 |
| `examples/dnp_mas/cross_effect_mas_steady.m` | `cross_effect_mas_steady()` | A MAS DNP simulation performed as described in Fred Mentink- Vigier's paper (Spinach rotation conventions are different) | 112 |
| `examples/dnp_mas/solid_effect_mas_dynam.m` | `solid_effect_mas_dynam()` | A MAS DNP simulation performed as described in Fred Mentink- Vigier's paper (Spinach rotation conventions are different) | 90 |
| `examples/dnp_mas/solid_effect_mas_enlev.m` | `solid_effect_mas_enlev()` | A MAS DNP simulation performed as described in Fred Mentink- Vigier's paper (Spinach rotation conventions are different) | 62 |
| `examples/dnp_mas/solid_effect_mas_powder.m` | `solid_effect_mas_powder()` | A MAS DNP simulation performed as described in Fred Mentink- Vigier's paper (Spinach rotation conventions are different) | 65 |
| `examples/dnp_mas/solid_effect_mas_steady.m` | `solid_effect_mas_steady()` | A MAS DNP simulation performed as described in Fred Mentink- Vigier's paper (Spinach rotation conventions are different) | 108 |
| `examples/dnp_sol/beam_dnp/beam_contact_curve.m` | `beam_contact_curve()` | The transformation of -E_z into I_z during the contact time of the BEAM DNP experiment. Further information in: Calculat | 60 |
| `examples/dnp_sol/beam_dnp/beam_field_profile.m` | `beam_field_profile()` | Field profile of a BEAM DNP experiment. <I_z> after a fixed contact time is calculated as a function of electron pulse a | 81 |
| `examples/dnp_sol/beam_dnp/beam_parameter_scan.m` | `beam_parameter_scan()` | 2D parameter scan of a BEAM DNP experiment. <I_z> after a set contact time is calculated as a function of electron pulse | 99 |
| `examples/dnp_sol/cross_effect_field_scan_1.m` | `cross_effect_field_scan_1()` | Magnetic field sweep cross effect DNP experiment -steady-state proton magnetisation under microwave iradiation as a func | 82 |
| `examples/dnp_sol/cross_effect_freq_scan_1.m` | `cross_effect_freq_scan_1()` | A simple TOTAPOL based Cross Effect DNP system. Set to repro- duce Figure 2c from Intensity differences are due to a dif | 72 |
| `examples/dnp_sol/cross_effect_freq_scan_2.m` | `cross_effect_freq_scan_2()` | A simple TOTAPOL based Cross Effect DNP system. Set to repro- duce Figure 2a from Intensity differences are due to a dif | 72 |
| `examples/dnp_sol/cross_effect_freq_scan_3.m` | `cross_effect_freq_scan_3()` | A simple TOTAPOL based Cross Effect DNP system. Set to repro- duce Figure 2b from Intensity differences are due to a dif | 72 |
| `examples/dnp_sol/crosspol_powder_static_1.m` | `crosspol_powder_static_1()` | E-15N cross-polarization experiment in the doubly rotating frame. Static powder simulation. Calculation time: seconds | 53 |
| `examples/dnp_sol/novel_dnp/novel_contact_curve.m` | `novel_contact_curve()` | The transformation of -E_z into I_z during the contact time of the NOVEL Solid Effect DNP experiment. Further informatio | 62 |
| `examples/dnp_sol/novel_dnp/novel_field_profile.m` | `novel_field_profile()` | Field profile of a NOVEL DNP experiment. <I_z> on 1H after a 0.25 us contact time is calculated as a function of electro | 83 |
| `examples/dnp_sol/novel_dnp/novel_parameter_scan.m` | `novel_parameter_scan()` | 2D parameter scan of a NOVEL DNP experiment. <I_z> on 1H after a 0.25 us contact time is calculated as a function of ele | 104 |
| `examples/dnp_sol/solid_effect_field_scan_1.m` | `solid_effect_field_scan_1()` | Magnetic field sweep DNP experiment involving a gadolinium ion, steady-state polarisation of a 15N nucleus is computed a | 72 |
| `examples/dnp_sol/solid_effect_freq_scan_1.m` | `solid_effect_freq_scan_1()` | A scan through the microwave frequency range in a steady state DNP experiment for a single 15N labelled urea mole- cule  | 89 |
| `examples/dnp_sol/solid_effect_freq_scan_2.m` | `solid_effect_freq_scan_2()` | A scan through the microwave frequency range in a steady state DNP experiment for a single 15N labelled urea mole- cule  | 89 |
| `examples/dnp_sol/solid_effect_timedep_1.m` | `solid_effect_timedep_1()` | A simulation of solid effect DNP for a tilted linear chain of three protons positioned at distances 7, 10 and 14 Angstro | 91 |
| `examples/dnp_sol/solid_effect_timedep_2.m` | `solid_effect_timedep_2()` | A simulation of solid effect DNP for a tilted linear chain of protons positioned at distances of 4+n^2 Angstrom with n=2 | 101 |
| `examples/dnp_sol/steady_state/novel_x_rep_time_ensemble_b1.m` | `novel_x_rep_time_ensemble_b1()` | Simulation of NOVEL DNP repetition time scan in the steady state with distributions in microwave B1 field. Calculation t | 130 |
| `examples/dnp_sol/steady_state/novel_x_rep_time_ensemble_b1_r.m` | `novel_x_rep_time_ensemble_b1_r()` | Simulation of NOVEL DNP repetition time scan in the steady state with distributions in electron-proton distance and micr | 138 |
| `examples/dnp_sol/steady_state/novel_x_rep_time_ensemble_r.m` | `novel_x_rep_time_ensemble_r()` | Simulation of NOVEL DNP repetition time scan in the steady state with distributions in electron-proton distance. Calcula | 123 |
| `examples/dnp_sol/steady_state/novel_x_rep_time_single.m` | `novel_x_rep_time_single()` | Simulation of NOVEL DNP repetition time scan in the steady state. Calculation time: minutes. | 114 |
| `examples/dnp_sol/steady_state/top_q_con_time_ensemble_b1.m` | `top_q_con_time_ensemble_b1()` | Simulation of TOP DNP contact time dependence in the steady state with electron Rabi frequency ensemble. Calculation tim | 156 |
| `examples/dnp_sol/steady_state/top_q_con_time_ensemble_b1_r.m` | `top_q_con_time_ensemble_b1_r()` | Simulation of TOP DNP contact time dependence in the steady state with electron-proton distance and elec- tron Rabi freq | 160 |
| `examples/dnp_sol/steady_state/top_q_con_time_ensemble_r.m` | `top_q_con_time_ensemble_r()` | Simulation of TOP DNP contact time dependence in the steady state with electron-proton distance ensembles. Calculation t | 126 |
| `examples/dnp_sol/steady_state/top_q_con_time_single.m` | `top_q_con_time_single()` | Simulation of TOP DNP contact time dependence in the steady state. Calculation time: hours. | 111 |
| `examples/dnp_sol/steady_state/top_q_nutation_ensemble_b1_r.m` | `top_q_nutation_ensemble_b1_r()` | Simulation of nutation frequency dependence of TOP DNP field profiles in the steady state with electron-proton distance  | 131 |
| `examples/dnp_sol/steady_state/top_q_rep_time_ensemble_b1.m` | `top_q_rep_time_ensemble_b1()` | Simulation of TOP DNP repetition time scan in the steady state with distributions in microwave B1 field. Calculation tim | 114 |
| `examples/dnp_sol/steady_state/top_q_rep_time_ensemble_b1_r.m` | `top_q_rep_time_ensemble_b1_r()` | Simulation of TOP DNP repetition time scan in the steady state with distributions in electron-proton distance and microw | 121 |
| `examples/dnp_sol/steady_state/top_q_rep_time_ensemble_r.m` | `top_q_rep_time_ensemble_r()` | Simulation of TOP DNP repetition time scan in the steady state with distributions in electron-proton distance. Calculati | 109 |
| `examples/dnp_sol/steady_state/top_q_rep_time_single.m` | `top_q_rep_time_single()` | Simulation of TOP DNP repetition time scan in the steady state. Calculation time: seconds. | 101 |
| `examples/dnp_sol/steady_state/tppm_q_con_time_ensemble_b1.m` | `tppm_q_con_time_ensemble_b1()` | Simulation of TPPM DNP contact time dependence in the steady state with electron Rabi frequency ensemble. Calculation ti | 117 |
| `examples/dnp_sol/steady_state/tppm_q_con_time_ensemble_b1_r.m` | `tppm_q_con_time_ensemble_b1_r()` | Simulation of TPPM DNP contact time dependence in the steady state with electron-proton distance and elec- tron Rabi fre | 124 |
| `examples/dnp_sol/steady_state/tppm_q_con_time_ensemble_r.m` | `tppm_q_con_time_ensemble_r()` | Simulation of TPPM DNP contact time dependence in the steady state with electron-proton distance ensembles. Calculation  | 111 |
| `examples/dnp_sol/steady_state/tppm_q_con_time_single.m` | `tppm_q_con_time_single()` | Simulation of TPPM DNP contact time dependence in the steady state. Calculation time: minutes. | 94 |
| `examples/dnp_sol/steady_state/tppm_q_rep_time_ensemble_b1.m` | `tppm_q_rep_time_ensemble_b1()` | Simulation of TPPM DNP repetition time scan in the steady state with distributions in microwave B1 field. Calculation ti | 113 |
| `examples/dnp_sol/steady_state/tppm_q_rep_time_ensemble_b1_r.m` | `tppm_q_rep_time_ensemble_b1_r()` | Simulation of TPPM DNP repetition time scan in the steady state with distributions in electron-proton distance and micro | 120 |
| `examples/dnp_sol/steady_state/tppm_q_rep_time_ensemble_r.m` | `tppm_q_rep_time_ensemble_r()` | Simulation of TPPM DNP repetition time scan in the steady state with distributions in electron-proton distance. Calculat | 108 |
| `examples/dnp_sol/steady_state/tppm_q_rep_time_single.m` | `tppm_q_rep_time_single()` | Simulation of TPPM DNP repetition time scan in the steady state. Calculation time: seconds. | 100 |
| `examples/dnp_sol/steady_state/xix_q_con_time_ensemble_b1.m` | `xix_q_con_time_ensemble_b1()` | Simulation of XiX DNP contact time dependence in the steady state with electron Rabi frequency ensemble. Calculation tim | 116 |
| `examples/dnp_sol/steady_state/xix_q_con_time_ensemble_b1_r.m` | `xix_q_con_time_ensemble_b1_r()` | Simulation of XiX DNP contact time dependence in the steady state with electron-proton distance and elec- tron Rabi freq | 123 |
| `examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r.m` | `xix_q_con_time_ensemble_r()` | Simulation of XiX DNP contact time dependence in the steady state with electron-proton distance ensembles. Calculation t | 111 |
| `examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T1e.m` | `xix_q_con_time_ensemble_r_T1e()` | Simulation of T1e dependence of XiX DNP contact curves in the steady state with electron-proton distance ensemble. Calcu | 131 |
| `examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T1n.m` | `xix_q_con_time_ensemble_r_T1n()` | Simulation of T1n dependence of XiX DNP contact curves in the steady state with electron-proton distance ensemble. Calcu | 129 |
| `examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T2e.m` | `xix_q_con_time_ensemble_r_T2e()` | Simulation of T2e dependence of XiX DNP contact curves in the steady state with electron-proton distance ensemble. Calcu | 131 |
| `examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T2n.m` | `xix_q_con_time_ensemble_r_T2n()` | Simulation of T2n dependence of XiX DNP contact curves in the steady state with electron-proton distance ensemble. Calcu | 131 |
| `examples/dnp_sol/steady_state/xix_q_con_time_single.m` | `xix_q_con_time_single()` | Simulation of XiX DNP contact time dependence in the steady state. Calculation time: minutes. | 94 |
| `examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_b1.m` | `xix_q_field_profile_ensemble_b1()` | Simulation of XiX DNP field profile in the steady state with electron Rabi frequency ensemble averaging. Calculation tim | 100 |
| `examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_b1_r.m` | `xix_q_field_profile_ensemble_b1_r()` | Simulation of XiX DNP field profile in the steady state with averaging over electron-proton distance and electron Rabi f | 111 |
| `examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r.m` | `xix_q_field_profile_ensemble_r()` | Simulation of XiX DNP field profile in the steady state with electron-proton distance ensemble averaging. Calculation ti | 96 |
| `examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T1e.m` | `xix_q_field_profile_ensemble_r_T1e()` | Simulation of T1e dependence of XiX DNP field profiles in the steady state with electron-proton distance ensemble. Calcu | 116 |
| `examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T1n.m` | `xix_q_field_profile_ensemble_r_T1n()` | Simulation of T1n dependence of XiX DNP field profiles in the steady state with electron- proton distance ensemble. Calc | 114 |
| `examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T2e.m` | `xix_q_field_profile_ensemble_r_T2e()` | Simulation of T2e dependence of XiX DNP field profiles in the steady state with electron-proton distance ensemble. Calcu | 116 |
| `examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T2n.m` | `xix_q_field_profile_ensemble_r_T2n()` | Simulation of T2n dependence of XiX DNP field profiles in the steady state with electron- proton distance ensemble. Calc | 116 |
| `examples/dnp_sol/steady_state/xix_q_field_profile_single.m` | `xix_q_field_profile_single()` | Simulation of XiX DNP field profile in the steady state, a single spin system without ensemble averaging. Calculation ti | 84 |
| `examples/dnp_sol/steady_state/xix_q_nutation_ensemble_b1_r.m` | `xix_q_nutation_ensemble_b1_r()` | Simulation of nutation frequency dependence of XiX DNP field profiles in the steady state with electron-proton distance  | 134 |
| `examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_b1.m` | `xix_q_rep_time_ensemble_b1()` | Simulation of XiX DNP repetition time scan in the steady state with distributions in microwave B1 field. Calculation tim | 113 |
| `examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_b1_r.m` | `xix_q_rep_time_ensemble_b1_r()` | Simulation of XiX DNP repetition time scan in the steady state with distributions in electron-proton distance and microw | 120 |
| `examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r.m` | `xix_q_rep_time_ensemble_r()` | Simulation of XiX DNP repetition time scan in the steady state with distributions in electron-proton distance. Calculati | 108 |
| `examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r_T1e.m` | `xix_q_rep_time_ensemble_r_T1e()` | Simulation of T1e dependent XiX DNP optimisation of repetition time in the steady state with electron- proton distance e | 128 |
| `examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r_T1n.m` | `xix_q_rep_time_ensemble_r_T1n()` | Simulation of T1n dependent XiX DNP optimisation of repetition time in the steady state with electron- proton distance e | 126 |
| `examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r_T2e.m` | `xix_q_rep_time_ensemble_r_T2e()` | Simulation of T2e dependence of XiX DNP repetition time profile in the steady state with electron-proton distan- ce ense | 128 |
| `examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r_T2n.m` | `xix_q_rep_time_ensemble_r_T2n()` | Simulation of T2n dependent XiX DNP optimisation of repetition time in the steady state with electron- proton distance e | 128 |
| `examples/dnp_sol/steady_state/xix_q_rep_time_single.m` | `xix_q_rep_time_single()` | Simulation of XiX DNP repetition time scan in the steady state. Calculation time: seconds. | 100 |
| `examples/dnp_sol/steady_state/xix_w_field_profile_ensemble_b1.m` | `xix_w_field_profile_ensemble_b1()` | Simulation of XiX DNP field profile in the steady state with electron Rabi frequency ensemble averaging. Calculation tim | 98 |
| `examples/dnp_sol/steady_state/xix_w_field_profile_ensemble_b1_r.m` | `xix_w_field_profile_ensemble_b1_r()` | Simulation of XiX DNP field profile in the steady state with averaging over electron-proton distance and electron Rabi f | 109 |
| `examples/dnp_sol/steady_state/xix_w_field_profile_ensemble_r.m` | `xix_w_field_profile_ensemble_r()` | Simulation of XiX DNP field profile in the steady state with electron-proton distance ensemble averaging. Calculation ti | 94 |
| `examples/dnp_sol/steady_state/xix_w_field_profile_single.m` | `xix_w_field_profile_single()` | Simulation of XiX DNP field profile in the steady state, a single spin system without ensemble averaging. Calculation ti | 82 |
| `examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_b1.m` | `xix_w_pulse_dur_ensemble_b1()` | 2D parameter scan of XiX DNP in the steady state with electron Rabi frequency ensemble. Calculation time: hours. | 115 |
| `examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_b1_r.m` | `xix_w_pulse_dur_ensemble_b1_r()` | 2D parameter scan of XiX DNP in the steady state with electron-proton distance and electron Rabi frequency ensembles. Ca | 126 |
| `examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_r.m` | `xix_w_pulse_dur_ensemble_r()` | 2D parameter scan of XiX DNP in the steady state with electron-proton distance ensemble. Calculation time: hours. | 113 |
| `examples/dnp_sol/steady_state/xix_w_pulse_dur_single.m` | `xix_w_pulse_dur_single()` | 2D parameter scan of XiX DNP in the steady state. Calculation time: hours. | 102 |
| `examples/dnp_sol/top_dnp/top_contact_curve.m` | `top_contact_curve()` | The transformation of -E_z into I_z during the contact time of the time-optimised pulsed DNP experiment. Further informa | 63 |
| `examples/dnp_sol/top_dnp/top_field_profile.m` | `top_field_profile()` | Field profile of a TOP DNP experiment. <I_z> after a fixed contact time is calculated as a function of electron pulse am | 82 |
| `examples/dnp_sol/top_dnp/top_parameter_scan.m` | `top_parameter_scan()` | 2D parameter scan of a TOP DNP experiment. <I_z> after a set contact time is calculated as a function of electron pulse  | 100 |
| `examples/dnp_sol/tppm_dnp/tppm_contact_curve.m` | `tppm_contact_curve()` | The transformation of -E_z into I_z during the contact time of the TPPM DNP experiment. Further information in: Calculat | 61 |
| `examples/dnp_sol/tppm_dnp/tppm_field_profile.m` | `tppm_field_profile()` | Field profile of a TPPM DNP experiment. <I_z> after a fixed contact time is calculated as a function of electron pulse a | 82 |
| `examples/dnp_sol/tppm_dnp/tppm_parameter_scan.m` | `tppm_parameter_scan()` | 2D parameter scan of a TPPM DNP experiment. <I_z> after a set contact time is calculated as a function of electron pulse | 100 |
| `examples/dnp_sol/xix_dnp/xix_contact_curve.m` | `xix_contact_curve()` | The transformation of -E_z into I_z during the contact time of the X-inverse-X DNP experiment. Further information in: C | 61 |
| `examples/dnp_sol/xix_dnp/xix_field_profile.m` | `xix_field_profile()` | Field profile of a XiX DNP experiment. <I_z> after a fixed contact time is calculated as a function of electron pulse am | 82 |
| `examples/dnp_sol/xix_dnp/xix_parameter_scan.m` | `xix_parameter_scan()` | 2D parameter scan of a XiX DNP experiment. <I_z> after a set contact time is calculated as a function of electron pulse  | 100 |
| `examples/esr_liq_pulsed/data_import/gaussian_import_example.m` | `gaussian_import_example()` | Methyl radical simulation, Gaussian import. The uncommon signal intensity pattern comes from g-HFC cross-correlation. | 54 |
| `examples/esr_liq_pulsed/data_import/orca_import_example.m` | `orca_import_example()` | Methyl radical simulation, ORCA import. The uncommon signal intensity pattern comes from g-HFC cross-correlation. | 65 |
| `examples/esr_liq_pulsed/endor_benzoquinone.m` | `endor_benzoquinone()` | CW ENDOR on 2-methoxy-1,4-benzoquinone radical in liquid state. Set to reproduce Figure 2 in http://dx.doi.org/10.1002/m | 57 |
| `examples/esr_liq_pulsed/endor_methyl.m` | `endor_methyl()` | Mims ENDOR spectrum of a methyl radical in liquid state. Magnetic parameters taken from a DFT calculation. Calculation t | 53 |
| `examples/esr_liq_pulsed/endor_nitroxide.m` | `endor_nitroxide()` | Mims ENDOR on a 15N-labelled nitroxide radical in liquid state. Magnetic parameters taken from a DFT calculation. Calcul | 53 |
| `examples/esr_liq_pulsed/endor_phenyl.m` | `endor_phenyl()` | Mims ENDOR on a phenyl radical in liquid state. The g-factor and the isotropic proton hyperfine couplings are specified  | 69 |
| `examples/esr_liq_pulsed/pulse_acquire_benzoquinone.m` | `pulse_acquire_benzoquinone()` | Pulse-acquire FFT ESR on 2-methoxy-1,4-benzoquinone radical in liquid state. Set to reproduce Figure 1 in Simple common  | 75 |
| `examples/esr_liq_pulsed/pulse_acquire_biaryl.m` | `pulse_acquire_biaryl()` | A time-domain pulse-acquire version of the EasySpin biaryl test file, with acknowledgements to Stefan Stoll. The Spinach | 85 |
| `examples/esr_liq_pulsed/pulse_acquire_chrysene.m` | `pulse_acquire_chrysene()` | W-band pulse-acquire FFT ESR spectrum of a chrysene cation radical in a non-viscous liquid. Simple common line width is  | 68 |
| `examples/esr_liq_pulsed/pulse_acquire_methyl.m` | `pulse_acquire_methyl()` | X-band pulse-acquire FFT ESR spectrum of methyl radical. Simple common line width is used as a relaxation model. Set to  | 64 |
| `examples/esr_liq_pulsed/pulse_acquire_phenyl.m` | `pulse_acquire_phenyl()` | W-band pulse-acquire FFT ESR spectrum of phenyl radical. Simple fixed line width is used as a relaxation model. Calculat | 61 |
| `examples/esr_liq_pulsed/rapidscan_nitroxide.m` | `rapidscan_nitroxide()` | Rapid scan ESR spectrum of a nitroxide radical. Calculation time: seconds | 55 |
| `examples/esr_liq_pulsed/relaxation_bisnitroxide.m` | `relaxation_bisnitroxide()` | X-band pulse-acquire FFT ESR spectrum of a bisnitroxide radical, using explicit time domain simulation with Redfield rel | 79 |
| `examples/esr_liq_pulsed/relaxation_fremysalt.m` | `relaxation_fremysalt()` | Pulse-acquire FFT ESR version of the EasySpin Fremy salt test file, with acknowledgements to Stefan Stoll. The Spinach s | 72 |
| `examples/esr_liq_pulsed/relaxation_nitroxide.m` | `relaxation_nitroxide()` | W-band pulse-acquire FFT ESR spectrum of a nitroxide radical, using explicit time domain simulation with Redfield relaxa | 60 |
| `examples/esr_liq_pulsed/relaxation_parafluoronitrobenzene.m` | `relaxation_parafluoronitrobenzene()` | A pulse-acquire FFT version of the EasySpin parafluoronitrobenzene test file, with acknowledgements to Stefan Stoll. The | 85 |
| `examples/esr_liq_pulsed/relaxation_parafluorotoluene.m` | `relaxation_parafluorotoluene()` | X-band pulse-acquire FFT ESR spectrum of parafluorotoluene radical, simulated using explicit time-domain propagation inc | 65 |
| `examples/esr_sol_pulsed/endor_davies_nox_crystal.m` | `endor_davies_nox_crystal()` | Davies ENDOR simulation for a nitroxide radical at a single orientation. Soft pulses are simulated using Fokker-Planck f | 116 |
| `examples/esr_sol_pulsed/endor_davies_nox_powder.m` | `endor_davies_nox_powder()` | Davies ENDOR simulation for a nitroxide radical. Soft pulses are simulated using Fokker-Planck formalism. This is a pain | 79 |
| `examples/esr_sol_pulsed/endor_mims_bdpa.m` | `endor_mims_bdpa()` | Mims ENDOR pulse sequence on BDPA with ideal electron pulses, reproducing Figure 10 from Calculation time: hours, much f | 67 |
| `examples/esr_sol_pulsed/endor_mims_echo_bdpa.m` | `endor_mims_echo_bdpa()` | Stimulated echo stage of the Mims ENDOR pulse sequence on BDPA. The nuclear pulse is not applied, this is echo dia- gnos | 57 |
| `examples/esr_sol_pulsed/endor_mims_nox_powder.m` | `endor_mims_nox_powder()` | Mims ENDOR simulation for a nitroxide radical powder. Ideal hard pulses are assumed. Calculation time: seconds. | 61 |
| `examples/esr_sol_pulsed/eseem_methyl_crystal.m` | `eseem_methyl_crystal()` | Two-pulse X-band ESEEM spectrum of a methyl radical at a specific orien- tation relative to the lab frame. Magnetic para | 59 |
| `examples/esr_sol_pulsed/eseem_nitroxide_crystal.m` | `eseem_nitroxide_crystal()` | Two-pulse X-band ESEEM spectrum of a nitroxide radical at a specific orientation relative to the lab frame. Magnetic par | 68 |
| `examples/esr_sol_pulsed/eseem_nitroxide_powder.m` | `eseem_nitroxide_powder()` | Powder-averaged two-pulse ESEEM on a 14N nitroxide radical. Time-domain simulation in Liouville space with powder averag | 69 |
| `examples/esr_sol_pulsed/eseem_phenyl_crystal.m` | `eseem_phenyl_crystal()` | Two-pulse X-band ESEEM spectrum of a phenyl radical at a specific orientation relative to the lab frame. Magnetic parame | 59 |
| `examples/esr_sol_pulsed/hard_3_pulse_deer_cu.m` | `hard_3_pulse_deer_cu()` | Three-pulse DEER on a Cu(II)-NO two electron system at X-band. The numerical calculation is done by brute-force time pro | 72 |
| `examples/esr_sol_pulsed/hard_3_pulse_deer_exchange.m` | `hard_3_pulse_deer_exchange()` | Three-pulse DEER on a Cu(II)-Cu(II) system in a linked porphyrin complex with a strong exchange coupling between the ele | 88 |
| `examples/esr_sol_pulsed/hard_3_pulse_deer_gd_1.m` | `hard_3_pulse_deer_gd_1()` | Gadolinium(III) DEER experiment at W-band using ideal pulses. Set to reproduce Figure 2b from the paper by Otting and co | 94 |
| `examples/esr_sol_pulsed/hard_3_pulse_deer_gd_2.m` | `hard_3_pulse_deer_gd_2()` | Gadolinium(III) DEER experiment. The calculation is done by brute- force time propagation and powder averaging. Outermos | 91 |
| `examples/esr_sol_pulsed/hard_3_pulse_deer_no.m` | `hard_3_pulse_deer_no()` | Nitroxide spin label DEER experiment at X-band. Two nitroxide radicals are positioned at a distance of 25 Angstroms. The | 56 |
| `examples/esr_sol_pulsed/hard_3_pulse_echo_cu.m` | `hard_3_pulse_echo_cu()` | Three-pulse DEER echo on a Cu(II)-NO two electron system at X-band. The calculation is done by brute-force time propagat | 57 |
| `examples/esr_sol_pulsed/hard_3_pulse_echo_gd.m` | `hard_3_pulse_echo_gd()` | Gadolinium(III) DEER echo experiment. The calculation is done by brute-force time propagation and powder averaging. Oute | 68 |
| `examples/esr_sol_pulsed/hard_3_pulse_echo_no.m` | `hard_3_pulse_echo_no()` | DEER spin echo for a pair of nitroxide radicals at X-band. Two nit- roxide radicals are positioned at a distance of 25 A | 58 |
| `examples/esr_sol_pulsed/holeburn_gd_dota_powder.m` | `holeburn_gd_dota_powder()` | A hole burning simulation for a gadolinium ion. The soft pulse is simulated using Fokker-Planck formalism. Zero-field sp | 93 |
| `examples/esr_sol_pulsed/holeburn_nitroxide_powder.m` | `holeburn_nitroxide_powder()` | A hole burning simulation for a nitroxide radical. The soft pulse is simulated using Fokker-Planck formalism; it is foll | 91 |
| `examples/esr_sol_pulsed/hpa_gd_dota_powder.m` | `hpa_gd_dota_powder()` | Powder averaged W-band pulsed ESR spectrum of Gd(III) DOTA complex. Ideal pulse with a large numerical powder grid is us | 61 |
| `examples/esr_sol_pulsed/hpa_nitroxide_powder.m` | `hpa_nitroxide_powder()` | Powder averaged pulse-acquire W-band Fourier ESR spectrum of nitroxide radical. An ideal pulse is assumed. Calculation t | 70 |
| `examples/esr_sol_pulsed/hpa_triplet.m` | `hpa_triplet()` | Hypothetical powder averaged X-band pulse-acquire ESR spectrum of photogenerated pentacene triplet state. Calculation ti | 73 |
| `examples/esr_sol_pulsed/hyscore_nitroxide_powder.m` | `hyscore_nitroxide_powder()` | Powder-averaged HYSCORE on a 14N nitroxide radical. Time-domain simulation in Liouville space. Set to reproduce Figure 2 | 65 |
| `examples/esr_sol_pulsed/oop_eseem_photosystem_1.m` | `oop_eseem_photosystem_1()` | Powder-averaged two-pulse out-of-phase ESEEM on the [P700+,A1-] spin-correlated electron pair in Photosystem I. Time-dom | 75 |
| `examples/esr_sol_pulsed/ridme_cu_nitroxide.m` | `ridme_cu_nitroxide()` | RIDME on a Cu(II)-NO two electron system at Q-band. The numerical calculation is done by brute-force time propaga- tion  | 90 |
| `examples/esr_sol_pulsed/sifter_nitroxide_powder.m` | `sifter_nitroxide_powder()` | An example of the SIFTER sequence. Calculation time: minutes. | 72 |
| `examples/esr_sol_pulsed/soft_3_pulse_deer_2e.m` | `soft_3_pulse_deer_2e()` | Three-pulse DEER simulation for a two-electron system. Soft pulses are simulated using the Fokker-Planck formalism. Calc | 79 |
| `examples/esr_sol_pulsed/soft_3_pulse_deer_3e.m` | `soft_3_pulse_deer_3e()` | Three-pulse DEER simulation for a three-electron system. Soft pulses are simulated using the Fokker-Planck formalism. Ca | 82 |
| `examples/esr_sol_pulsed/soft_4_pulse_deer_2e.m` | `soft_4_pulse_deer_2e()` | Four-pulse DEER simulation for a two-electron system. Soft pulses are simulated using the Fokker-Planck formalism. Calcu | 80 |
| `examples/esr_sol_pulsed/soft_4_pulse_deer_3e.m` | `soft_4_pulse_deer_3e()` | Four-pulse DEER simulation for a three-electron system. Soft pulses are simulated using the Fokker-Planck formalism. Cal | 83 |
| `examples/esr_sol_pulsed/spa_gd_dota_powder.m` | `spa_gd_dota_powder()` | A soft pulse simulation for a gadolinium ion. The soft pulse is simulated using Fokker-Planck formalism. Zero-field spli | 83 |
| `examples/esr_sol_pulsed/spa_nitroxide_powder.m` | `spa_nitroxide_powder()` | A soft pulse simulation for a nitroxide radical powder. The soft pulse is simulated using the Fokker-Planck formalism; i | 73 |
| `examples/esr_sol_swept/fieldsweep_gd_dota.m` | `fieldsweep_gd_dota()` | Powder averaged W-band field-swept ESR spectrum of Gd(III) DOTA complex. Exact diagonalisation is used. Calculation time | 52 |
| `examples/esr_sol_swept/fieldsweep_nitroxide.m` | `fieldsweep_nitroxide()` | Field swept EPR spectrum of nitroxide, computed by finding resonance fields and transition moments. Calculation time: se | 56 |
| `examples/esr_sol_swept/fieldsweep_porphyrin.m` | `fieldsweep_porphyrin()` | Field swept EPR spectrum of copper porphyrin complex, computed by finding resonance fields and transition moments. Calcu | 69 |
| `examples/esr_sol_swept/fieldsweep_triplet.m` | `fieldsweep_triplet()` | Powder averaged X-band field-swept ESR spectrum of photo- generated pentacene triplet state. Calculation time: seconds. | 64 |
| `examples/esr_sol_swept/temperature_gd_dota.m` | `temperature_gd_dota()` | Powder averaged W-band field-swept ESR spectrum of Gd(III) DOTA complex. Exact diagonalisation is used and a tempera- tu | 67 |
| `examples/extremes/difluoroheptane.m` | `difluoroheptane()` | 19F NMR spectrum of anti-3,4-difluoroheptane (16 spins) by explicit time-domain evolution in Liouville space. WARNING: n | 122 |
| `examples/extremes/fluoroisooctane.m` | `fluoroisooctane()` | A deliberately adversarial example from Art Bochevarov at Schodinger Inc. In this case, IK-2 approximation in Liou- vill | 96 |
| `examples/extremes/high_symmetry_1.m` | `high_symmetry_1()` | 1H NMR spectrum of a large and highly symmetric spin system with two tert-butyl groups supplied by Eberhard Matern. Done | 89 |
| `examples/extremes/high_symmetry_2.m` | `high_symmetry_2()` | 31P NMR spectrum of a large and highly symmetric spin system with two tert-butyl groups supplied by Eberhard Matern. Don | 89 |
| `examples/extremes/perfluoropyrene.m` | `perfluoropyrene()` | X-band pulsed ESR spectrum of perfluoropyrene cation radical, computed using brute force operator algebra in the full 4, | 64 |
| `examples/extremes/ph_enc_3d_highres.m` | `ph_enc_3d_highres()` | Slice selection in 3D followed by phase-encoded imaging of the resulting slice. This simulation fills up a sys- tem with | 124 |
| `examples/extremes/phosphorus_cluster.m` | `phosphorus_cluster()` | Phosphorus system simulation for Gerhard Hagele. Done by brute force Liouville space time propagation. WARNING: needs 32 | 112 |
| `examples/fitting/fluoroalkanes/anti_difluorobutane.m` | `anti_difluorobutane()` | Fitting of 1H NMR spectrum of anti-2,3-difluorobutane with respect to J-couplings. See our paper for further details: Ca | 158 |
| `examples/fitting/fluoroalkanes/anti_difluoroheptane.m` | `anti_difluoroheptane()` | Fitting of 1H NMR spectrum of anti-3,5-difluoroheptane with respect to J-couplings. See our paper for further details: M | 208 |
| `examples/fitting/fluoroalkanes/anti_difluoropentane.m` | `anti_difluoropentane()` | Fitting of 1H NMR spectrum of syn-2,4-difluoropentane with respect to J-couplings. See our paper for further details: Ca | 193 |
| `examples/fitting/fluoroalkanes/difluoropropane.m` | `difluoropropane()` | Fitting of 1H and 19F NMR spectrums of 1,3-difluoropropane with respect to J-couplings. See our paper for further detail | 175 |
| `examples/fitting/fluoroalkanes/fluorobutane.m` | `fluorobutane()` | Fitting of 1H NMR spectrum of 2-fluoropentane with respect to J-couplings. See our paper for further details: Calculatio | 196 |
| `examples/fitting/fluoroalkanes/syn_difluorobutane.m` | `syn_difluorobutane()` | Fitting of 1H NMR spectrum of syn-2,3-difluorobutane with respect to J-couplings. See our paper for further details: Cal | 155 |
| `examples/fitting/fluoroalkanes/syn_difluoroheptane.m` | `syn_difluoroheptane()` | Fitting of 1H NMR spectrum of syn-3,5-difluoroheptane with respect to J-couplings. See our paper for further details: Me | 208 |
| `examples/fitting/fluoroalkanes/syn_difluoropentane.m` | `syn_difluoropentane()` | Fitting of 1H NMR spectrum of syn-2,4-difluoropentane with respect to J-couplings. See our paper for further details: Ca | 193 |
| `examples/fitting/fumarate_global.m` | `fumarate_global()` | Simultaneous fitting of 1H and 13C NMR spectra of a slightly asymmetric fumarate diester. Calculation time: minutes | 134 |
| `examples/fitting/maleate_global.m` | `maleate_global()` | Simultaneous fitting of 1H and 13C NMR spectra of a slightly asymmetric maleate diester. Calculation time: hours | 134 |
| `examples/fitting/nmr_kinetics/glucose_exsy_a.m` | `glucose_exsy_a()` | Fitting of 2,2,3,3-tetrafluoroglucose NOESY with respect to the reaction rates in a chemical exchange problem and the ro | 201 |
| `examples/fitting/nmr_kinetics/glucose_exsy_b.m` | `glucose_exsy_b()` | Fitting of 3,3-difluoroglucose NOESY with respect to the reaction rates in a chemical exchange and the rotational correl | 172 |
| `examples/fitting/nmr_rdc/saupe_example.m` | `saupe_example()` | Extracting Saupe order matrix from NH RDC data. Experimental measurements kindly provided by Andras Boeszoermenyi, Thiba | 50 |
| `examples/fundamentals/clebsch_gordan/cg_test_large.m` | `cg_test_large()` | Compares the output of Spinach Clebsch-Gordan function with the arbitrary precision results returned by Mathematica. | 35 |
| `examples/fundamentals/clebsch_gordan/cg_test_small.m` | `cg_test_small()` | Compares the output of Spinach Clebsch-Gordan function with the arbitrary precision results returned by Mathematica. | 35 |
| `examples/fundamentals/convention_tests/blinv_test.m` | `blinv_test()` | Internal consistency of Blicharski invariants and their polarisation relationships with inner products of irre- ducible  | 43 |
| `examples/fundamentals/convention_tests/euler_sup_test.m` | `euler_sup_test()` | Euler angle superposition tests. | 100 |
| `examples/fundamentals/convention_tests/invariants.m` | `invariants()` | A test of Equation 3 in http://dx.doi.org/10.1002/chem.200902300 | 60 |
| `examples/fundamentals/convention_tests/nqi_test.m` | `nqi_test()` | Test of the reverse decomposition of spin-1 Hamiltonians. | 39 |
| `examples/fundamentals/convention_tests/rotations_1.m` | `rotations_1()` | Tests the internal consistency of kernel rotation functions. | 129 |
| `examples/fundamentals/convention_tests/rotations_2.m` | `rotations_2()` | A rotations test comparing the Hamiltonians for a manually rotated (at the interaction specification level) spin system  | 80 |
| `examples/fundamentals/convention_tests/spsk_test.m` | `spsk_test()` | Test of the span-skew interaction convention. | 34 |
| `examples/fundamentals/convention_tests/stevens_test.m` | `stevens_test()` | Tests of Spinach Stevens operator function against explicit expressions from the literature. | 123 |
| `examples/fundamentals/convention_tests/tensors.m` | `tensors()` | Test the conversion from Stevens operator coefficients to irreducible spherical tensor coefficients. | 45 |
| `examples/fundamentals/correlation_function_1.m` | `correlation_function_1()` | Computes rotational correlation functions using a Monte-Carlo method and compares them to the analytical results returne | 96 |
| `examples/fundamentals/correlation_function_2.m` | `correlation_function_2()` | Computes rotational correlation functions using a Monte-Carlo method and compares them to the analytical results returne | 96 |
| `examples/fundamentals/correlation_function_3.m` | `correlation_function_3()` | Computes rotational correlation functions using a Monte-Carlo method and compares them to the analytical results returne | 97 |
| `examples/fundamentals/correlation_function_4.m` | `correlation_function_4()` | Computes rotational correlation functions using a Monte-Carlo method and compares them to the analytical results returne | 96 |
| `examples/fundamentals/correlation_function_5.m` | `correlation_function_5()` | Computes the following rotational correlation function G(k,m,p,q)=<R(k,m)*R(p,q)> where R is the 3D Cartesian rotation m | 50 |
| `examples/fundamentals/derivative_tests/auxmat_test.m` | `auxmat_test()` | Testing IK's favourite equation numerically -auxiliary matrix expression against a high-accuracy finite diffe- rence app | 61 |
| `examples/fundamentals/derivative_tests/difdiff_bs_rect.m` | `difdiff_bs_rect()` | Directional derivative test for Cartesian GRAPE with Bloch-Siegert corrections. | 151 |
| `examples/fundamentals/derivative_tests/dirdiff_1.m` | `dirdiff_1()` | Test of matrix exponential differentiation routines. Analytical derivatives are compared to central finite differences. | 58 |
| `examples/fundamentals/derivative_tests/dirdiff_2.m` | `dirdiff_2()` | Test of matrix exponential differentiation of second order Magnus product quadrature (trapdiff.m) with the result com- p | 56 |
| `examples/fundamentals/derivative_tests/dirdiff_3_rect.m` | `dirdiff_3_rect()` | Directional derivative test for the Cartesian GRAPE module, rectangles integrator. | 82 |
| `examples/fundamentals/derivative_tests/dirdiff_3_trap.m` | `dirdiff_3_trap()` | Directional derivative test for the Cartesian GRAPE module, trapezium integrator. | 82 |
| `examples/fundamentals/derivative_tests/dirdiff_4_rect.m` | `dirdiff_4_rect()` | Directional derivative test for the phase-modulated GRAPE module, rectangles integrator. | 83 |
| `examples/fundamentals/derivative_tests/dirdiff_4_trap.m` | `dirdiff_4_trap()` | Directional derivative test for the phase-modulated GRAPE module, trapezium integrator. | 83 |
| `examples/fundamentals/derivative_tests/dirdiff_5_rect.m` | `dirdiff_5_rect()` | GRAPE Hessian internal consistency test: Newton against Goodwin algorithm. | 71 |
| `examples/fundamentals/derivative_tests/dirdiff_6_rect.m` | `dirdiff_6_rect()` | GRAPE Hessian test against finite-differenced gradients. | 85 |
| `examples/fundamentals/derivative_tests/dirdiff_7_rect.m` | `dirdiff_7_rect()` | Directional derivative test for the phase-modulated GRAPE module, with an ensemble of chained linear and non-linear inst | 88 |
| `examples/fundamentals/derivative_tests/dirdiff_8_rect.m` | `dirdiff_8_rect()` | GRAPE phase Hessian test against finite-differenced gradients, rectangles integrator. | 88 |
| `examples/fundamentals/derivative_tests/dirdiff_test_system.m` | `[spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalism)` | Spin system generator for directional derivative tests. Syntax: [spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(forma | 80 |
| `examples/fundamentals/exchange_coupling/yamaguchi.m` | `yamaguchi()` | Yamaguchi equation estimate of exchange coupling from a broken-symmetry DFT calculation on a bistrityl bira- dical with  | 18 |
| `examples/fundamentals/high_spin_system_1.m` | `high_spin_system_1()` | Pulse-acquire NMR spectrum in a system with a hypothetical scalar coupling to a 235U nucleus. The spectral lines should  | 52 |
| `examples/fundamentals/nuclear_structure/woods_saxon.m` | `woods_saxon(mass_number,level_number)` | A loose implementation of single-nucleon Hamiltonian eigenfunction calculation in the three-dimensional Woods-Saxon pote | 58 |
| `examples/fundamentals/nutation_dist_test.m` | `nutation_dist_test()` | Recovery of an RF field distribution from a nutation curve measured with the same coil used for excitation and detection | 98 |
| `examples/fundamentals/operator_tests/commutation_1.m` | `commutation_1()` | Commutators of simple operators and superoperators. The test calculation is performed three times in the three formalism | 42 |
| `examples/fundamentals/operator_tests/commutation_2.m` | `commutation_2()` | Commutators of simple operators and superoperators. The test calculation is performed three times in the three formalism | 42 |
| `examples/fundamentals/operator_tests/commutation_3.m` | `commutation_3()` | Commutators of simple operators and superoperators. The test calculation is performed three times in the three formalism | 51 |
| `examples/fundamentals/operator_tests/commutation_4.m` | `commutation_4()` | Commutators of simple operators and superoperators. The test calculation is performed three times in the three formalism | 44 |
| `examples/fundamentals/operator_tests/commutation_5.m` | `commutation_5()` | Tests of basic product and commutation relations between bosonic operators. | 62 |
| `examples/fundamentals/operator_tests/commutation_6.m` | `commutation_6()` | Commutation and product relations for Pauli and central transition operators. | 91 |
| `examples/fundamentals/operator_tests/commutation_7.m` | `commutation_7()` | Expansion relations for operator basis transforms. | 126 |
| `examples/fundamentals/operator_tests/commutation_8.m` | `commutation_8()` | Delicate action and commutation tests for Hilbert-Liouville conversion and Stevens operators. | 59 |
| `examples/fundamentals/operator_tests/expansions.m` | `expansions()` | Product action tests for irreducible spherical tensor operators and orthogonalised bosonic monomials. | 83 |
| `examples/fundamentals/paramag_test.m` | `paramag_test()` | Paramagnetic chemical shift module tests. | 32 |
| `examples/fundamentals/perturb_theory.m` | `perturb_theory()` | Rayleigh-Schrodinger and Van Vleck perturbation theory modules test. Eigenvector representations differ in the two theor | 81 |
| `examples/fundamentals/pfg_test_1.m` | `pfg_test_1()` | A test of the explicit gradient pulse function that uses the auxiliary matrix formalism to compute sample volume integra | 77 |
| `examples/fundamentals/pfg_test_2.m` | `pfg_test_2()` | Demonstrate the use of the auxiliary matrix algorithm in generating a gradient sandwich multiple-quantum filter. For fur | 86 |
| `examples/fundamentals/quadratures/aht_benchmark.m` | `aht_benchmark()` | Accuracy benchmark for the Hamiltonian period propagator caclulation using Lie group integrators. Bacause the mo- dulati | 135 |
| `examples/fundamentals/quadratures/expmint2_test.m` | `expmint2_test()` | Verification of expmint2 against numerical integration. | 48 |
| `examples/fundamentals/quadratures/grid_diagrams.m` | `grid_diagrams()` | Spherical grid diagrams for IK's book. | 55 |
| `examples/fundamentals/quadratures/grid_leb_vs_etc.m` | `grid_leb_vs_etc()` | Heuristic vs Lebedev spherical quadrature bake-off, illus- trating the fact that, well... heuristic grids suck. | 68 |
| `examples/fundamentals/quadratures/grid_quality.m` | `grid_quality()` | Performance analysis for the spherical and SO(3) integration grids supplied with Spinach kernel. | 90 |
| `examples/fundamentals/quadratures/mas_benchmark.m` | `mas_benchmark()` | Integrating the Lioville -von Neumann equation through one period of the MAS rotor using the piecewise-constant Hamilton | 129 |
| `examples/fundamentals/quadratures/product_quadratures_1.m` | `product_quadratures_1()` | Accuracy test for Lie-group product quadratures as a function of discretisation step in the E1000B Veshtort-Griffin puls | 135 |
| `examples/fundamentals/quadratures/product_quadratures_2.m` | `product_quadratures_2()` | A test of Lie-group product quadratures on a chirped frequency oscillator with radiation damping that has a state-depend | 140 |
| `examples/fundamentals/roof_effect.m` | `roof_effect()` | Roof effect in a strongly J-coupled two-spin system. | 67 |
| `examples/fundamentals/spin_lock.m` | `spin_lock()` | A spin-locking experiment on a two-spin system. | 66 |
| `examples/fundamentals/state_spaces_1.m` | `state_spaces_1()` | Correlation order dynamics in a pulse-acquire experiment on strychnine. Set to reproduce Figure 4 from our state space r | 60 |
| `examples/fundamentals/state_spaces_2.m` | `state_spaces_2()` | Contributions from different orders of spin correlation to the system trajectory in the pulse-acquire 1H NMR simulation  | 116 |
| `examples/fundamentals/state_spaces_3.m` | `state_spaces_3()` | Transverse magnetisation dynamics in a pulse-acquire experiment on a fatty acid. This example looks at how the magnetisa | 64 |
| `examples/fundamentals/state_spaces_4.m` | `state_spaces_4()` | Trajectory analysis for a MAS simulation of isotopically labelled glycine powder, starting from L+ on protons. Calculati | 56 |
| `examples/fundamentals/state_tests/normalization_1.m` | `normalization_1()` | Internal consistency test for the state vectors and matrices in each of the three formalisms supported by Spinach. | 47 |
| `examples/fundamentals/state_tests/normalization_2.m` | `normalization_2()` | Internal consistency test for the state vectors and matrices. Checks that the inner products are consistent between the  | 52 |
| `examples/fundamentals/state_tests/state_consistency_1.m` | `state_consistency_1()` | Test of internal consistency for state and operator generation across the three formalisms supported by Spinach. Two-and | 76 |
| `examples/fundamentals/state_tests/state_consistency_2.m` | `state_consistency_2()` | Test of consistency in the projection between spherical tensor basis set and Zeeman basis set. | 72 |
| `examples/fundamentals/state_tests/state_consistency_3.m` | `state_consistency_3()` | Deuterium pair singlet, triplet, and quintet state internal consistency test. | 81 |
| `examples/fundamentals/state_tests/thermal_equilibrium_1.m` | `thermal_equilibrium_1()` | Observables at thermal equilibrium using the three formalisms supported by Spinach kernel, tested against known answers. | 57 |
| `examples/fundamentals/state_tests/thermal_equilibrium_2.m` | `thermal_equilibrium_2()` | Thermal equilibrium states, using all the different formalisms supported by Spinach kernel. | 47 |
| `examples/fundamentals/state_tests/thermal_equilibrium_3.m` | `thermal_equilibrium_3()` | Test of the thermal equilibrium functionality against the textbook expressions for the Boltzmann populations. | 73 |
| `examples/fundamentals/state_tests/thermal_equilibrium_4.m` | `thermal_equilibrium_4()` | Test of the invariance of the thermal equilibrium state under the thermalised relaxation superoperator. | 72 |
| `examples/fundamentals/state_tests/thermal_equilibrium_5.m` | `thermal_equilibrium_5()` | Cross-formalism test of state recovery towards the thermodynamic equilibrium. | 86 |
| `examples/fundamentals/strong_coupling.m` | `strong_coupling()` | A garden variety strongly coupled two-spin system. | 53 |
| `examples/fundamentals/symmetry_1.m` | `symmetry_1()` | Liouvillian symmetrization for a radical pair with four equivalent nuclei under the S4 permutation group. | 52 |
| `examples/fundamentals/symmetry_2.m` | `symmetry_2()` | Pulse-acquire NMR spectrum of a highly symmetric spin system provided by Andres Castillo. Uses the fully sym- metric irr | 68 |
| `examples/fundamentals/symmetry_3.m` | `symmetry_3()` | 1H NMR spectrum of valine. Uses the fully symmetric irreducible representation of the S3(x)S3 group. | 57 |
| `examples/fundamentals/symmetry_4.m` | `symmetry_4()` | Hamiltonian symmetrization for a radical pair with four equivalent nuclei under the S4 permutation group. | 50 |
| `examples/fundamentals/tensor_structures/amensolve_test_1.m` | `amensolve_test_1()` | Detailed unit test for ttclass/amensolve against dense references. The test builds structured positive-definite tensor-t | 226 |
| `examples/fundamentals/tensor_structures/amensum_test_1.m` | `amensum_test_1()` | Detailed unit test for ttclass/amensum against dense references. The test uses buffered rank-one tensor trains, compares | 150 |
| `examples/fundamentals/tensor_structures/kronm_test_1.m` | `kronm_test_1()` | Tests for the kron-times-matrix infrastructure. | 47 |
| `examples/fundamentals/tensor_structures/polyadic_test_1.m` | `polyadic_test_1()` | Unit tests for the polyadic object. | 120 |
| `examples/fundamentals/tensor_structures/polyadic_test_2.m` | `polyadic_test_2()` | Unit tests for advanced polyadic functionality. | 110 |
| `examples/fundamentals/tensor_structures/ttrain_test_1.m` | `ttrain_test_1()` | A simple test of ttclass object arithmetic. | 31 |
| `examples/giant_spin/dy_lft_single_1.m` | `dy_lft_single_1()` | Reproduction of MOLCAS results with the Ligand Field Theory model for a single Dy(III) ion. Calculation time: seconds | 97 |
| `examples/giant_spin/dy_lft_single_2.m` | `dy_lft_single_2()` | A demonstration that most lanthanide complexes are in the ZFS limit for the purposes of relaxation theory. One of the fi | 91 |
| `examples/giant_spin/lanthanide_powder.m` | `lanthanide_powder()` | Powder spectrum of Gd(III) with ZFS up to 3rd spherical rank using the giant spin Hamiltonian formalism in a sweepable 4 | 50 |
| `examples/giant_spin/lanthanide_redfield.m` | `lanthanide_redfield()` | Relaxation rate of Gd(III) as a function of zero-field splitting, computed using Redfield theory. The correla- tion time | 63 |
| `examples/giant_spin/nuclear_relaxation_1.m` | `nuclear_relaxation_1()` | Nuclear relaxation rates using the adiabatic elimination method for a rapidly relaxing Dy(III) ion with a user-specified | 148 |
| `examples/giant_spin/quartet_levels.m` | `quartet_levels()` | Energy levels magnetic field scan for a spin-3/2 particle with a zero-field splitting. Calculation time: seconds | 42 |
| `examples/giant_spin/quartet_magn.m` | `quartet_magn()` | Sample magnetisation during a finite-speed magnetic field sweep for a spin-3/2 particle with a zero-field splitting. Cal | 52 |
| `examples/giant_spin/triple_dy_eqmag_field.m` | `triple_dy_eqmag_field()` | Simulation of the field dependence of the magnetisation of a triple- Dy triangular complex -see Figure S27 and S28 in th | 148 |
| `examples/giant_spin/triple_dy_eqmag_temp.m` | `triple_dy_eqmag_temp()` | Simulation of the temperature dependence of the magnetisation of a triple-Dy triangular complex -see Figure S27 and S28  | 152 |
| `examples/giant_spin/triple_dy_levels.m` | `triple_dy_levels()` | Eight lowest energy levels as a function of the applied magnetic fi- eld in a triple Dy triangular complex -see Figure 1 | 124 |
| `examples/giant_spin/triple_dy_magn.m` | `triple_dy_magn()` | Simulation of a finite-speed magnetic field sweep experiment for a single crystal of a triple-Dy triangular complex in a | 133 |
| `examples/giant_spin/triple_tb_eqmag_field.m` | `triple_tb_eqmag_field()` | Simulation of the field dependence of the magnetisation of a triple- Tb triangular complex -see Figure S27 and S28 in th | 222 |
| `examples/giant_spin/triple_tb_eqmag_temp.m` | `triple_tb_eqmag_temp()` | Simulation of the temperature dependence of the magnetisation of a triple-Tb triangular complex -see Figure S27 and S28  | 226 |
| `examples/imaging/bright_fat_effect_cpmg.m` | `bright_fat_effect_cpmg()` | Bright fat effect under CPMG echo train -magnetisation losses are greater in MRI experiments on J-coupled systems becaus | 90 |
| `examples/imaging/bright_fat_effect_udd.m` | `bright_fat_effect_udd()` | Bright fat effect under UDD echo train -magnetisation losses are greater in MRI experiments on J-coupled systems because | 89 |
| `examples/imaging/diffusion_weighted_2d.m` | `diffusion_weighted_2d()` | 2D diffusion weighted image with an arbitrary geometric pattern serving as diffusion coefficient distribution. Simulatio | 80 |
| `examples/imaging/diffusion_weighted_epi_2d.m` | `diffusion_weighted_epi_2d()` | 2D echo planar imaging example in the presence of istropic diffusion. Stejskal-Tanner SE echo planar diffusion-weighted  | 95 |
| `examples/imaging/diffusion_weighted_epi_3d.m` | `diffusion_weighted_epi_3d()` | Three-dimensional echo planar imaging in the presence of realistic diffusion. Simulation time: hours, faster with a Tesl | 133 |
| `examples/imaging/dpfgse_signal_select.m` | `dpfgse_signal_select()` | DPFGSE signal selection example for a solution of GABA in water. Gradients and soft pulses are done explicitly. Simulati | 99 |
| `examples/imaging/dpfgse_signal_suppress.m` | `dpfgse_signal_suppress()` | DPFGSE water suppression example for a solution of GABA in water. Gradients and soft pulses are done explicitly. Simulat | 99 |
| `examples/imaging/echo_planar_2d.m` | `echo_planar_2d()` | Echo planar imaging example in 2D for a brain phantom. Simulation time: seconds, faster with a Tesla V100 GPU. | 90 |
| `examples/imaging/echo_planar_3d.m` | `echo_planar_3d()` | Slice selection in 3D followed by three-dimensional echo planar imaging sequence. Simulation time: hours, faster with a  | 125 |
| `examples/imaging/fast_echo_2d.m` | `fast_echo_2d()` | Fast (in the experiment duration sense) spin echo 2D brain imaging example. Simulation time: hours, faster with a Tesla  | 83 |
| `examples/imaging/gradient_echo_1d.m` | `gradient_echo_1d()` | A gradient echo experiment in the presence of diffusion and flow. Calculation time: seconds. Ahmed Allami Ilya Kuprov | 64 |
| `examples/imaging/gradient_image_1d.m` | `gradient_image_1d()` | 1D imaging experiment with a hard pulse in the presence of diffusion and flow. Calculation time: seconds. Ahmed Allami I | 78 |
| `examples/imaging/phase_encoding_2d.m` | `phase_encoding_2d()` | Simple phase-encoded 2D imaging example. Calculation time: seconds. Ahmed Allami Ilya Kuprov | 75 |
| `examples/imaging/phase_encoding_3d.m` | `phase_encoding_3d()` | Slice selection in 3D followed by phase-encoded imaging of the resulting slice. Simulation time: minutes, faster with a  | 123 |
| `examples/imaging/press_1d_example.m` | `press_1d_example()` | 1D PRESS example. Three independent spin systems are localised in three areas of a 1D sample. The areas are slectively e | 115 |
| `examples/imaging/press_2d_example.m` | `press_2d_example()` | 2D PRESS example. Three independent spin systems are localised in three spots of a 2D sample. The spots are slectively e | 113 |
| `examples/imaging/press_3d_example.m` | `press_3d_example()` | PRESS excitation profile in three dimensions with tilted gradient system. Change the frequency under Pulse Parame- ters  | 79 |
| `examples/imaging/slice_select_1d_shaped.m` | `slice_select_1d_shaped()` | Slice selection example using a one-dimensional sample and a shaped slice selection pulse in the presence of diffusion a | 88 |
| `examples/imaging/slice_select_1d_square.m` | `slice_select_1d_square()` | Slice selection example using a one-dimensional sample and a square slice selection pulse in the presence of diffusion a | 88 |
| `examples/imaging/spin_echo_1d.m` | `spin_echo_1d()` | A spin echo experiment under a gradient in the presence of diffusion and flow. Calculation time: seconds. Ahmed Allami I | 65 |
| `examples/imaging/spiral_2d.m` | `spiral_2d()` | Spiral K-space imaging example in 2D. Calculation time: minutes. Ahmed Allami Ilya Kuprov | 85 |
| `examples/karplus_curves/leu_chi_fit.m` | `leu_chi_fit()` | Karplus coefficients extraction from a DFT dihedral angle scan over one of the chi angles in leucine using Gaussian09. | 18 |
| `examples/karplus_curves/tyr_chi_fit.m` | `tyr_chi_fit()` | Karplus coefficients extraction from a DFT dihedral angle scan over one of the chi angles in tyrosine using Gaussian09. | 18 |
| `examples/kinetics/aziridine_exsy_1.m` | `aziridine_exsy_1()` | NOESY/EXSY experiment on phenylaziridine, including scalar relaxation of the second kind induced by the 14N nucleus, in  | 194 |
| `examples/kinetics/aziridine_exsy_2.m` | `aziridine_exsy_2()` | NOESY/EXSY experiment on phenylaziridine, including scalar relaxation of the second kind induced by the 14N nucleus, in  | 205 |
| `examples/kinetics/exchange_asymmetric.m` | `exchange_asymmetric()` | Two-spin asymmetric chemical exchange pattern. Calculation time: seconds. | 51 |
| `examples/kinetics/exchange_symmetric.m` | `exchange_symmetric()` | Two-spin symmetric chemical exchange pattern. Calculation time: seconds. | 51 |
| `examples/kinetics/flux_asymmetric.m` | `flux_asymmetric()` | Two-spin asymmetric magnetization flux problem. Calculation time: seconds. | 52 |
| `examples/kinetics/flux_symmetric.m` | `flux_symmetric()` | Two-spin symmetric magnetization flux problem. Calculation time: seconds. | 52 |
| `examples/kinetics/frydman_pump_a.m` | `frydman_pump_a()` | Lucio Frydman's water exchange based spin-lock pump, Figure 2 from https://doi.org/10.1016/j.jmr.2021.107083 Calculation | 132 |
| `examples/kinetics/frydman_pump_b.m` | `frydman_pump_b()` | Lucio Frydman's water exchange based spin-lock pump, Figure 9 from https://doi.org/10.1016/j.jmr.2021.107083 Calculation | 173 |
| `examples/kinetics/glucose_exsy_a.m` | `glucose_exsy_a()` | 2D EXSY of transmembrane exchange of 2,2,3,3-tetrafluoroglucose. See the fitting example set for the script that yielded | 168 |
| `examples/kinetics/glucose_exsy_b.m` | `glucose_exsy_b()` | 2D EXSY of transmembrane exchange of 3,3-difluoroglucose. See the fitting example set for the script that yielded the pa | 137 |
| `examples/kinetics/mas_exchange_1.m` | `mas_exchange_1()` | Two-site position exchange for a deuterium nucleus. The sites differ in the chemical shift and the orientation of the qu | 69 |
| `examples/kinetics/nonlinear/diels_alder_spec.m` | `diels_alder_spec()` | Repeated pulse-acquire experiment during the Diels-Alder cyclo- addition of acetylene to butadiene, demonstrating the no | 265 |
| `examples/kinetics/nonlinear/diels_alder_zmag.m` | `diels_alder_zmag()` | Time-domain Z magnetisation dynamics in the Diels-Alder cycloaddition of acetylene to butadiene, demonstrating the non-l | 170 |
| `examples/kinetics/relayed_hyperpol.m` | `relayed_hyperpol()` | Relayed NOE from hyperpolarized water to ALA-GLY dipeptide, generating Figure S7 from Christopher Pötzl | 101 |
| `examples/liquid_crystals/rdc_fourspin.m` | `rdc_fourspin()` | CLIP-HSQC spectrum of a four-spin system in a liquid crystal with a user-specified order matrix. Calculation times: seco | 70 |
| `examples/liquid_crystals/rdc_twospin.m` | `rdc_twospin()` | CLIP-HSQC spectrum of a C-H system in a liquid crystal with a user-specified order matrix. Calculation time: seconds. | 64 |
| `examples/microfluidics/plain_diff.m` | `plain_diff()` | Simple diffusion simulation without spin dynamics. Longitudinal magnetisation is tracked as a function of time. | 115 |
| `examples/microfluidics/plain_flow.m` | `plain_flow()` | Simple flow simulation with no dynamics in the spin subspace: longitudinal magnetisation is tracked as a function of tim | 116 |
| `examples/microfluidics/plain_nmr.m` | `plain_nmr()` | NMR spectrum of the reaction mixture in the absence of chemical kinetics and spatial dynamics. | 52 |
| `examples/microfluidics/plain_reaction.m` | `plain_reaction()` | Non-linear reaction kinetics in a situation when there is no hydrodynamics, diffusion, or spin dynamics. This is in- ten | 52 |
| `examples/microfluidics/reacting_flow.m` | `reacting_flow()` | Flow in the absence of spin dynamics, but presence of two unidirectional second-order chemical reactions. Simulation tim | 140 |
| `examples/microfluidics/reacting_flow_nmr.m` | `reacting_flow_nmr()` | Complete microfluidic simulation: diffusion, flow, two second- order chemical reactions, and NMR detection in a narrow s | 252 |
| `examples/microfluidics/reacting_nmr.m` | `reacting_nmr()` | Non-linear reaction kinetics in combination with spin evolution (repeated pulse-acquire NMR) and relaxation (Redfield th | 213 |
| `examples/microfluidics/show_mesh.m` | `show_mesh()` | Import, Voronoi tessellation, and plotting of the hydrodynamic mesh and velocity field from COMSOL. | 34 |
| `examples/nmr_diffusion/diffusion_test_1.m` | `diffusion_test_1()` | A standard diffusion equation solver with no spin dynamics present. Calculation time: seconds | 61 |
| `examples/nmr_diffusion/diffusion_test_2a.m` | `diffusion_test_2a()` | A standard diffusion equation solver with no spin dynamics present. Isotropic diffusion. Calculation time: minutes. | 62 |
| `examples/nmr_diffusion/diffusion_test_2b.m` | `diffusion_test_2b()` | A standard diffusion equation solver with no spin dynamics present. Non-uniform diffusion coefficient distribution. Calc | 64 |
| `examples/nmr_diffusion/diffusion_test_2c.m` | `diffusion_test_2c()` | A standard diffusion equation solver with no spin dynamics present. Anisotropic diffusion with a pe- riodic boundary con | 63 |
| `examples/nmr_diffusion/flow_test_1a.m` | `flow_test_1a()` | A standard diffusion and flow equation solver with no spin dynamics present and periodic boundary condition. Diffusion c | 62 |
| `examples/nmr_diffusion/flow_test_1b.m` | `flow_test_1b()` | A standard diffusion and flow equation solver with no spin dynamics present and periodic boundary condition. Diffusion c | 62 |
| `examples/nmr_diffusion/flow_test_2.m` | `flow_test_2()` | A combination of diffusion and flow in two dimensions with a periodic boundary condition. Calculation time: minutes. | 62 |
| `examples/nmr_diffusion/flow_test_3.m` | `flow_test_3()` | Circular flow in three-dimensional space in the absence of spin dynamics. Calculation time: minutes, faster on GPU. | 81 |
| `examples/nmr_diffusion/slosh_test_1.m` | `slosh_test_1()` | Probability density sloshing around a harmonic oscillator. Calculation time: seconds. | 36 |
| `examples/nmr_diffusion/slosh_test_2.m` | `slosh_test_2()` | Probability density sloshing around a harmonic oscillator in the presence of a gravitational pull twards the left. Calcu | 37 |
| `examples/nmr_liquids/clip_hsqc_camphor.m` | `clip_hsqc_camphor()` | CLIP-HSQC spectrum of camphor with natural content of 13C isotope. Coordinates, shielding anisotropies and J-couplings c | 87 |
| `examples/nmr_liquids/clip_hsqc_sucrose.m` | `clip_hsqc_sucrose()` | CLIP-HSQC spectrum of sucrose with natural content of 13C isotope. Coordinates, shielding anisotropies and J-couplings c | 91 |
| `examples/nmr_liquids/coloc_test.m` | `coloc_test()` | A simple COLOC pulse sequence example for a two-spin 1H-13C system with a long-range J-coupling. Calculation time: secon | 54 |
| `examples/nmr_liquids/cosy45_rotenone.m` | `cosy45_rotenone()` | Magnitude mode COSY-45 spectrum of rotenone. Calculation time: minutes | 83 |
| `examples/nmr_liquids/cosy90_derome.m` | `cosy90_derome()` | Figure 8.26 from Andrew Derome's "Modern NMR Techniques for Chemistry Research". Calculation time: seconds | 59 |
| `examples/nmr_liquids/cosy90_rotenone.m` | `cosy90_rotenone()` | COSY-90 spectrum of rotenone. Calculation time: minutes | 84 |
| `examples/nmr_liquids/cosy90_strychnine.m` | `cosy90_strychnine()` | COSY spectrum of strychnine. Calculation time: minutes | 56 |
| `examples/nmr_liquids/cosy90_sucrose.m` | `cosy90_sucrose()` | COSY spectrum of sucrose (magnetic parameters computed with DFT). Calculation time: minutes | 56 |
| `examples/nmr_liquids/crazed_test.m` | `crazed_test()` | Long range intermolecular coherences predicted by Warren and co-workers Calculation time: seconds | 60 |
| `examples/nmr_liquids/ct_cosy_2spins.m` | `ct_cosy_2spins()` | CT COSY spectrum for 2 spins. Calculation time: minutes | 53 |
| `examples/nmr_liquids/ct_cosy_rotenone.m` | `ct_cosy_rotenone()` | CT COSY spectrum of rotenone using the assignment reported in Calculation time: minutes | 82 |
| `examples/nmr_liquids/ct_cosy_three_spin.m` | `ct_cosy_three_spin()` | CT-COSY of three spin system. Calculation time: seconds | 55 |
| `examples/nmr_liquids/ct_hsqc_2spins.m` | `ct_hsqc_2spins()` | CT HSQC spectrum of 2 spin system Calculation time: seconds | 73 |
| `examples/nmr_liquids/ct_hsqc_strychnine.m` | `ct_hsqc_strychnine()` | CT HSQC spectrum of strychnine with natural content of 13C isotope. Calculation time: hours. | 76 |
| `examples/nmr_liquids/ct_hsqc_sucrose.m` | `ct_hsqc_sucrose()` | CT HSQC spectrum of sucrose with natural content of 13C isotope (magnetic parameters computed with DFT). Calculation tim | 87 |
| `examples/nmr_liquids/dept_strychnine.m` | `dept_strychnine()` | DEPT135 experiment on strychnine. Calculation time: minutes | 67 |
| `examples/nmr_liquids/deptq_strychnine.m` | `deptq_strychnine()` | DEPTQ135 experiment on strychnine. Calculation time: minutes | 67 |
| `examples/nmr_liquids/dqf_cosy_strychnine.m` | `dqf_cosy_strychnine()` | DQF-COSY spectrum of strychnine. Calculation time: minutes | 61 |
| `examples/nmr_liquids/dqf_cosy_sucrose.m` | `dqf_cosy_sucrose()` | DQF-COSY spectrum of sucrose (magnetic parameters computed with DFT). Calculation time: minutes | 62 |
| `examples/nmr_liquids/ecosy_strychnine.m` | `ecosy_strychnine()` | E.COSY spectrum of strychnine. Calculation time: minutes | 60 |
| `examples/nmr_liquids/gcosy90_strychnine.m` | `gcosy90_strychnine()` | Gradient-selected COSY spectrum of strychnine. Calculation time: minutes | 67 |
| `examples/nmr_liquids/hetcor_strychnine.m` | `hetcor_strychnine()` | HETCOR spectrum of strychnine with natural content of 13C isotope. Calculation time: minutes | 70 |
| `examples/nmr_liquids/hetcor_sucrose.m` | `hetcor_sucrose()` | HETCOR spectrum of sucrose with natural content of 13C isotope (magnetic parameters computed with DFT). Calculation time | 81 |
| `examples/nmr_liquids/hmbc_camphor.m` | `hmbc_camphor()` | HMBC spectrum of camphor with natural content of 13C isotope. Coordinates, shielding anisotropies and J-couplings comput | 74 |
| `examples/nmr_liquids/hmbc_cyprinol.m` | `hmbc_cyprinol()` | HMBC of cyprinol with natrual abundance of 13C isotope. Calculation time: seconds | 71 |
| `examples/nmr_liquids/hmbc_sucrose.m` | `hmbc_sucrose()` | HMBC spectrum of sucrose with natural content of 13C isotope (magnetic parameters computed with DFT). Calculation time:  | 80 |
| `examples/nmr_liquids/hmqc_cyprinol.m` | `hmqc_cyprinol()` | HMQC spectrum of cyprinol with natural abundance of 13C isotope. Calculation time: seconds | 73 |
| `examples/nmr_liquids/hmqc_strychnine.m` | `hmqc_strychnine()` | HMQC spectrum of strychnine with natural content of 13C isotope. Calculation time: minutes | 71 |
| `examples/nmr_liquids/hmqc_sucrose.m` | `hmqc_sucrose()` | HMQC spectrum of sucrose with natural content of 13C isotope (magnetic parameters computed with DFT). Calculation time:  | 81 |
| `examples/nmr_liquids/hoesy_camphor.m` | `hoesy_camphor()` | 13C{1H} HOESY spectrum of camphor with natural content of 13C isotope. Coordinates, shielding anisotropies and J-couplin | 89 |
| `examples/nmr_liquids/hoesy_ftyr_a.m` | `hoesy_ftyr_a()` | (1H) -> (19F) HOESY spectrum of fluorotyrosine, with the magneti- sation transfer direction picked so as to minimise the | 75 |
| `examples/nmr_liquids/hoesy_ftyr_b.m` | `hoesy_ftyr_b()` | (19F) -> (1H) HOESY spectrum of fluorotyrosine. This is not the right way to run this sequence in proteins because aro-  | 74 |
| `examples/nmr_liquids/hoesy_strychnine.m` | `hoesy_strychnine()` | 13C{1H} HOESY spectrum of strychnine with natural content of 13C isotope. Calculation time: minutes | 88 |
| `examples/nmr_liquids/hsqc_cyprinol.m` | `hsqc_cyprinol()` | HSQC spectrum of cyprinol with natural abundance of 13C isotope. Calculation time: seconds | 80 |
| `examples/nmr_liquids/hsqc_strychnine.m` | `hsqc_strychnine()` | HSQC spectrum of strychnine with natural content of 13C isotope. Calculation time: minutes | 78 |
| `examples/nmr_liquids/hsqc_sucrose.m` | `hsqc_sucrose()` | HSQC spectrum of sucrose with natural content of 13C isotope (magnetic parameters computed with DFT). Calculation time:  | 89 |
| `examples/nmr_liquids/inad_cyprinol.m` | `inad_cyprinol()` | INADEQUATE spectrum of cyprinol. The sequence selects double- quantum coherence from coupled 13C pairs and converts it b | 80 |
| `examples/nmr_liquids/inad_sucrose.m` | `inad_sucrose()` | INADEQUATE spectrum of sucrose. The sequence selects double- quantum coherence from coupled 13C pairs and converts it ba | 89 |
| `examples/nmr_liquids/inad_three_spin.m` | `inad_three_spin()` | INADEQUATE spectrum of a three-spin system with J-coupling between two spins only. The sequence selects double-quantum c | 50 |
| `examples/nmr_liquids/inad_three_spin_2d.m` | `inad_three_spin_2d()` | Example of a 2D-INADEQUATE spectrum of a generic three-spin system. Calculation time: seconds. Theresa Hune Christian Gr | 81 |
| `examples/nmr_liquids/inept_strychnine.m` | `inept_strychnine()` | INEPT experiment on strychnine. Calculation time: minutes | 66 |
| `examples/nmr_liquids/inv_rec_strychnine.m` | `inv_rec_strychnine()` | 1H inversion-recovery experiment on strychnine at 250 MHz. Calculation time: minutes | 62 |
| `examples/nmr_liquids/mqs_propanol.m` | `mqs_propanol()` | Multiple-quantum NMR experiment for a propanol spin system. Calculation time: seconds | 101 |
| `examples/nmr_liquids/mqs_six_spin.m` | `mqs_six_spin()` | Multiple-quantum (MQ) NMR experiment for a coupled system of six spins. Calculation time: minutes | 90 |
| `examples/nmr_liquids/noe_four_spin.m` | `noe_four_spin()` | Inversion-recovery NOE effect spectrum on a simple four-spin system, with the rightmost proton signal inverted and a pul | 84 |
| `examples/nmr_liquids/noe_strychnine.m` | `noe_strychnine()` | Inversion-recovery NOE effect spectrum on strychnine, with the rightmost proton signal inverted and a pulse-acquire expe | 84 |
| `examples/nmr_liquids/noe_two_spin_het.m` | `noe_two_spin_het()` | Nuclear overhauser effect in a heteronuclear two-spin system in the short correlation time case. Calculation time: secon | 56 |
| `examples/nmr_liquids/noe_two_spin_hom.m` | `noe_two_spin_hom()` | Nuclear overhauser effect in a homonuclear two-spin system in the long correlation time case. Calculation time: seconds | 56 |
| `examples/nmr_liquids/noe_zq_beats.m` | `noe_zq_beats()` | Zero-quantum beats in the Overhauser effect in a strongly coupled two-spin system. Calculation time: seconds | 60 |
| `examples/nmr_liquids/noesy_methanol.m` | `noesy_methanol()` | NOESY spectrum of 13C methanol. J-couplings from Pecul and Helgaker, CSA tensors from DFT. Note the presence of cross-pe | 88 |
| `examples/nmr_liquids/noesy_strychnine.m` | `noesy_strychnine()` | NOESY spectrum of strychnine. Calculation time: minutes | 71 |
| `examples/nmr_liquids/noesy_sucrose.m` | `noesy_sucrose()` | NOESY spectrum of sucrose (magnetic parameters computed with DFT). Calculation time: minutes | 72 |
| `examples/nmr_liquids/pa_difluoroheptane_anti.m` | `pa_difluoroheptane_anti()` | Pulse-acquire 1H NMR spectrum of anti-3,5-difluoroheptane with a manual basis set specification as a merger of Lie algeb | 127 |
| `examples/nmr_liquids/pa_difluoroheptane_syn.m` | `pa_difluoroheptane_syn()` | Pulse-acquire 1H NMR spectrum of syn-3,5-difluoroheptane with a manual basis set specification as a merger of Lie algebr | 127 |
| `examples/nmr_liquids/pa_menthol.m` | `pa_menthol()` | Menthol NMR spectrum from Damien Jeannerat, including the effect of bad Z1 and Z2 magnet shims. Calculation time: minute | 56 |
| `examples/nmr_liquids/pa_naphtopyranone.m` | `pa_naphtopyranone()` | NMR spectrum of 3-phenylmethylene-1H,3H-naphtho-[1,8-c,d]-pyran-1-one, magnetic parameters from: Calculation time: secon | 73 |
| `examples/nmr_liquids/pa_rotenone.m` | `pa_rotenone()` | 1H NMR spectrum of rotenone using T1/T2 relaxation model, magnetic parameters from: Calculation time: seconds | 91 |
| `examples/nmr_liquids/pa_strychnine.m` | `pa_strychnine()` | 1H NMR spectrum of strychnine including an accurate model of line widths via Redfield superoperator. Calculation time: s | 61 |
| `examples/nmr_liquids/pa_styrene.m` | `pa_styrene()` | 1H NMR spectrum of styrene, the calculation is performed in Hilbert space to demonstrate parallel propagation described  | 64 |
| `examples/nmr_liquids/pa_sucrose.m` | `pa_sucrose()` | 1H NMR spectrum of sucrose (magnetic parameters read in from a DFT calculation), including Redfield relaxation superoper | 61 |
| `examples/nmr_liquids/pansy_dual_ch.m` | `pansy_dual_ch()` | PANSY-COSY spectra of camphor with natural content of 13C isotope. Coordinates, shieldings, and J-couplings computed wit | 79 |
| `examples/nmr_liquids/pansy_triple_ch.m` | `pansy_triple_ch()` | Triple-channel PANSY experiment on glycine with natural content of 13C isotope. Coordinates, shieldings, and J- coupling | 111 |
| `examples/nmr_liquids/roesy_strychnine.m` | `roesy_strychnine()` | ROESY spectrum of strychnine. Calculation time: minutes | 69 |
| `examples/nmr_liquids/roesy_sucrose.m` | `roesy_sucrose()` | ROESY spectrum of sucrose (magnetic parameters computed with DFT). Calculation time: minutes | 70 |
| `examples/nmr_liquids/sat_rec_strychnine.m` | `sat_rec_strychnine()` | 1H saturation-recovery experiment on strychnine at 250 MHz. Calculation time: minutes | 62 |
| `examples/nmr_liquids/tocsy_sucrose.m` | `tocsy_sucrose()` | TOCSY spectrum of sucrose (magnetic parameters computed with DFT). Calculation time: seconds | 65 |
| `examples/nmr_metabol/molecule_a.m` | `molecule_a()` | 1H NMR spectrum of a molecule from the GISSMO database. Calculation time: seconds | 47 |
| `examples/nmr_metabol/molecule_b.m` | `molecule_b()` | 1H NMR spectrum of a molecule from the GISSMO database. Calculation time: seconds | 47 |
| `examples/nmr_metabol/molecule_c.m` | `molecule_c()` | 1H NMR spectrum of a molecule from the GISSMO database. Calculation time: seconds | 47 |
| `examples/nmr_nucleic/expt_data/rna_noesy_expt.m` | `rna_noesy_expt()` | Experimental NOESY spectrum of the Harvard RNA. Shunsuke Imai Scott Robson Gerhard Wagner Zenawi Welderufael Ilya Kuprov | 31 |
| `examples/nmr_nucleic/rna_hsqc_theo.m` | `rna_hsqc_theo()` | 1H-13C HSQC spectrum of the example RNA molecule provided by the Wagner group. Calculation time: minutes Shunsuke Imai S | 79 |
| `examples/nmr_nucleic/rna_noesy_theo.m` | `rna_noesy_theo()` | 1H-1H NOESY spectrum of the example RNA molecule provided by the Gerhard Wagner group at Harvard University. Calculation | 84 |
| `examples/nmr_overtone/cpmas_glycine_accum.m` | `cpmas_glycine_accum()` | Cross-polarization experiment between protons and 14N overtone transition in glycine under MAS. Glycine quadrupolar tens | 93 |
| `examples/nmr_overtone/cpmas_glycine_match_1.m` | `cpmas_glycine_match_1()` | Cross-polarization experiment between protons and 14N overtone transition in glycine under MAS. Glycine quadrupolar tens | 92 |
| `examples/nmr_overtone/cpmas_glycine_match_2.m` | `cpmas_glycine_match_2()` | Cross-polarization experiment between protons and 14N overtone transition in glycine under MAS. Glycine quadrupolar tens | 110 |
| `examples/nmr_overtone/cpmas_glycine_simple.m` | `cpmas_glycine_simple()` | Cross-polarization experiment between protons and 14N overtone transition in glycine under MAS. Glycine quadrupolar tens | 74 |
| `examples/nmr_overtone/cpmas_valine_accum.m` | `cpmas_valine_accum()` | Cross-polarization experiment between protons and 14N overtone transition in N-acetylvaline under MAS, computed using Fo | 93 |
| `examples/nmr_overtone/cpmas_valine_match_1.m` | `cpmas_valine_match_1()` | Cross-polarization experiment between protons and 14N overtone transition in N-acetylvaline under MAS, computed using Fo | 93 |
| `examples/nmr_overtone/cpmas_valine_match_2.m` | `cpmas_valine_match_2()` | Cross-polarization experiment between protons and 14N overtone transition in N-acetylvaline under MAS, computed using Fo | 110 |
| `examples/nmr_overtone/cpmas_valine_simple.m` | `cpmas_valine_simple()` | Cross-polarization experiment between protons and 14N overtone transition in N-acetylvaline under MAS, computed using Fo | 75 |
| `examples/nmr_overtone/dante_glycine.m` | `dante_glycine()` | 14N overtone DANTE spectrum of glycine, computed using Fokker-Planck formalism. Glycine quadrupolar tensor data comes fr | 73 |
| `examples/nmr_overtone/dor_glycine.m` | `dor_glycine()` | Panoramic double rotation overtone 14N spectrum of glycine, simulated as described in our paper (Figure 1B): A short pul | 78 |
| `examples/nmr_overtone/mas_boron_1.m` | `mas_boron_1()` | Overtone Z-detection 10B magic angle spinning NMR spectrum. The sample is spinning in the JEOL direction. Parameters fro | 54 |
| `examples/nmr_overtone/mas_boron_2.m` | `mas_boron_2()` | Overtone 10B magic angle spinning NMR spectrum. The sample is spinning in the JEOL direction. Parameters from Nghia Duon | 66 |
| `examples/nmr_overtone/mas_boron_3.m` | `mas_boron_3()` | Overtone 10B magic angle spinning NMR spectrum. The sample is spinning in the JEOL direction. Parameters from Nghia Duon | 54 |
| `examples/nmr_overtone/mas_glycine_1.m` | `mas_glycine_1()` | Overtone detection 14N magic angle spinning NMR spectrum of glycine, computed using Fokker-Planck formalism. Glycine qua | 76 |
| `examples/nmr_overtone/mas_valine_1.m` | `mas_valine_1()` | Overtone detection 14N magic angle spinning NMR spectrum of N-acetylvaline, computed using Fokker-Planck formalism. Va-  | 72 |
| `examples/nmr_overtone/mas_valine_2.m` | `mas_valine_2()` | Overtone Z-detection 14N magic angle spinning NMR spectrum of N-acetylvaline, computed using Fokker-Planck formalism. Va | 59 |
| `examples/nmr_overtone/powder_glycine_1.m` | `powder_glycine_1()` | Overtone detection 14N powder NMR spectrum of glycine, computed using Fokker-Planck formalism. Glycine quadrupolar tenso | 65 |
| `examples/nmr_overtone/powder_glycine_2.m` | `powder_glycine_2()` | Overtone detection 14N powder NMR spectrum of glycine, computed using Fokker-Planck formalism. Glycine quadrupolar tenso | 60 |
| `examples/nmr_paramag/calbindin/tm_1igv_distr_fit.m` | `tm_1igv_distr_fit()` | Inverse problem for the unpaired electron density distribution. Experimental data from Gottfried Ott- ing (Australian Na | 48 |
| `examples/nmr_paramag/calbindin/tm_1igv_lcurve.m` | `tm_1igv_lcurve()` | Inverse problem for the unpaired electron density distribution. Experimental data from Gottfried Ott- ing (Australian Na | 52 |
| `examples/nmr_paramag/calbindin/tm_1igv_mult_fit.m` | `tm_1igv_mult_fit()` | Electron location and susceptibility tensor recovery from experimental PCS data using point electron model. Experi- ment | 29 |
| `examples/nmr_paramag/calbindin/tm_1igv_point_fit.m` | `tm_1igv_point_fit()` | Electron location and susceptibility tensor recovery from experimental PCS data using point electron model. Experi- ment | 29 |
| `examples/nmr_paramag/carb_anh/s166c_kuprov.m` | `s166c_kuprov()` | Distributed fit for the S166C mutant dataset for human carbonic anhydrase II. The system and the method are described in | 52 |
| `examples/nmr_paramag/carb_anh/s166c_lcurve.m` | `s166c_lcurve()` | L-curves for the S166C mutant dataset for human carbonic anhydrase II. The system and the method are described in: A ste | 57 |
| `examples/nmr_paramag/carb_anh/s166c_mult.m` | `s166c_mult()` | Multipolar fit for the S166C mutant dataset for human carbonic anhydrase II. The system and the method are described in: | 38 |
| `examples/nmr_paramag/carb_anh/s166c_point.m` | `s166c_point()` | Point fit for the S166C mutant dataset for human carbonic anhydrase II. The system and the method are described in: A st | 38 |
| `examples/nmr_paramag/carb_anh/s217c_kuprov.m` | `s217c_kuprov()` | Distributed fit for the S217C mutant dataset for human carbonic anhydrase II. The system and the method are described in | 52 |
| `examples/nmr_paramag/carb_anh/s217c_lcurve.m` | `s217c_lcurve()` | L-curves for the S217C mutant dataset for human carbonic anhydrase II. The system and the method are described in: A ste | 57 |
| `examples/nmr_paramag/carb_anh/s217c_mult.m` | `s217c_mult()` | Multipolar fit for the S217C dataset mutant dataset for human carbonic anhydrase II. The system and the method are descr | 38 |
| `examples/nmr_paramag/carb_anh/s217c_point.m` | `s217c_point()` | Point fit for the S217C dataset mutant dataset for human carbonic anhydrase II. The system and the method are described  | 38 |
| `examples/nmr_paramag/carb_anh/s220c_kuprov.m` | `s220c_kuprov()` | Distributed fit for the S220C mutant dataset for human carbonic anhydrase II. The system and the method are described in | 52 |
| `examples/nmr_paramag/carb_anh/s220c_lcurve.m` | `s220c_lcurve()` | L-curves for the S220C mutant dataset for human carbonic anhydrase II. The system and the method are described in: A ste | 57 |
| `examples/nmr_paramag/carb_anh/s220c_mult.m` | `s220c_mult()` | Multipolar fit for the S220C mutant dataset for human carbonic anhydrase II. The system and the method are described in: | 38 |
| `examples/nmr_paramag/carb_anh/s220c_point.m` | `s220c_point()` | Point fit for the S220C mutant dataset for human carbonic anhydrase II. The system and the method are described in: A st | 38 |
| `examples/nmr_paramag/carb_anh/s50c_kuprov.m` | `s50c_kuprov()` | Distributed fit for the S50C mutant dataset for human carbonic anhydrase II. The system and the method are described in: | 52 |
| `examples/nmr_paramag/carb_anh/s50c_lcurve.m` | `s50c_lcurve()` | L-curves for the S50C mutant dataset for human carbonic anhydrase II. The system and the method are described in: A step | 57 |
| `examples/nmr_paramag/carb_anh/s50c_mult.m` | `s50c_mult()` | Multipolar fit for the S50C mutant dataset for human carbonic anhydrase II. The system and the method are described in:  | 38 |
| `examples/nmr_paramag/carb_anh/s50c_point.m` | `s50c_point()` | Point fit for the S50C mutant dataset for human carbonic anhydrase II. The system and the method are described in: A ste | 38 |
| `examples/nmr_paramag/combi_fit_1.m` | `combi_fit_1()` | Extracting the susceptibility tensor from DFT hyperfine tensors and experimental paramagnetic shifts. A combinatorial pr | 50 |
| `examples/nmr_paramag/combi_fit_2.m` | `combi_fit_2()` | Extracting the susceptibility tensor from DFT hyperfine tensors and experimental paramagnetic shifts. A combinatorial pr | 50 |
| `examples/nmr_paramag/dft_density.m` | `dft_density()` | Simulation of the pseudocontact shift field of the Europium(III) complex of 1,4,7,10-tetrakis(2-pyridylmethyl)-1,4,7,10- | 72 |
| `examples/nmr_paramag/gau_density.m` | `gau_density()` | Simple calculation of the PCS field of a Gaussian distribution of the electron probability density. | 30 |
| `examples/nmr_paramag/point_vs_distr.m` | `point_vs_distr()` | A comparison between a point model fit and a multipole model fit described in in a situation where the probability densi | 117 |
| `examples/nmr_paramag/porphyrin_example_1.m` | `porphyrin_example_1()` | Computing PCS using different models in basic Cu(II) and Co(II) porphyrin complexes. The metal is at the origin. See the | 46 |
| `examples/nmr_paramag/porphyrin_example_2.m` | `porphyrin_example_2()` | Computing PCS using different models in a basic Cu(II) porphyrin complex. The metal is at the origin. See the "getting s | 76 |
| `examples/nmr_paramag/porphyrin_example_3.m` | `porphyrin_example_3()` | Computes PCS using different models in basic Cu(II) porphyrin complex. See the "getting started" manual at The paper des | 74 |
| `examples/nmr_paramag/simple_pcs_1.m` | `simple_pcs_1()` | Pseudocontact shift and Curie relaxation on a proton due to the presence of a point magnetic susceptibility centre. Calc | 50 |
| `examples/nmr_proteins/ct_cosy_gb1.m` | `ct_cosy_gb1()` | Constant-time COSY experiment simulation for the GB1 protein. Simulation time: minutes, faster with a Tesla A100 GPU. | 67 |
| `examples/nmr_proteins/ct_hsqc_gb1.m` | `ct_hsqc_gb1()` | Constant-time HSQC experiment simulation for the GB1 protein. Simulation time: hours, faster with a Tesla A100 GPU. | 74 |
| `examples/nmr_proteins/expt_data/hnco_ubiquitin_expt.m` | `hnco_ubiquitin_expt()` | Experimental HNCO spectrum of human ubiquitin. Donghan Lee (Max Planck Institute) Ilya Kuprov (University of Southampton | 53 |
| `examples/nmr_proteins/expt_data/hsqc_ubiquitin_expt.m` | `hsqc_ubiquitin_expt()` | Experimental HSQC spectrum of human ubiquitin. Donghan Lee (Max Planck Institute) Ilya Kuprov (University of Southampton | 43 |
| `examples/nmr_proteins/expt_data/noesy_ubiquitin_expt.m` | `noesy_ubiquitin_expt()` | Experimental HNCO spectrum of human ubiquitin. Donghan Lee (Max Planck Institute) Ilya Kuprov (University of Southampton | 29 |
| `examples/nmr_proteins/hcanh_gb1.m` | `hcanh_gb1()` | Simulated H(CA)NH spectrum of GB1 protein. It is assumed that only the backbone is 13C,15N-labelled. Calculation time: m | 79 |
| `examples/nmr_proteins/hcanh_simple.m` | `hcanh_simple()` | A minimal example of H(CA)NH pulse sequence simulation. Calculation time: seconds. | 81 |
| `examples/nmr_proteins/hcch_cosy_gb1.m` | `hcch_cosy_gb1()` | 3D HCCH COSY experiment on GB1 protein. Calculation time: hours. | 87 |
| `examples/nmr_proteins/hcch_cosy_simple.m` | `hcch_cosy_simple()` | 3D HCCH COSY experiment on a small protein fragment. Calculation time: seconds. | 81 |
| `examples/nmr_proteins/hnca_gb1.m` | `hnca_gb1()` | Simulated HNCA spectrum of GB1 protein. It is assumed that only the backbone is 13C,15N-labelled. Calculation time: minu | 80 |
| `examples/nmr_proteins/hnca_simple.m` | `hnca_simple()` | A minimal example of HNCA pulse sequence simulation. Calculation time: seconds. | 77 |
| `examples/nmr_proteins/hncaco_gb1.m` | `hncaco_gb1()` | Simulated HNCACO spectrum of GB1 protein. It is assumed that only the backbone is 13C,15N-labelled. Calculation time: mi | 83 |
| `examples/nmr_proteins/hncaco_simple.m` | `hncaco_simple()` | A minimal example of HNCACO pulse sequence simulation. Calculation time: seconds. | 80 |
| `examples/nmr_proteins/hnco_simple.m` | `hnco_simple()` | A minimal example of HNCO pulse sequence simulation. Calculation time: seconds. | 77 |
| `examples/nmr_proteins/hnco_ubiquitin.m` | `hnco_ubiquitin()` | Theoretical HNCO of human ubiquitin. It is assumed that only the backbone is 13C,15N-labelled. Calculation time: minutes | 82 |
| `examples/nmr_proteins/hncoca_simple.m` | `hncoca_simple()` | A minimal HNCOCA pulse sequence simulation. Calculation time: seconds. | 76 |
| `examples/nmr_proteins/hncoca_ubiquitin.m` | `hncoca_ubiquitin()` | Theoretical HN(CO)CA of human ubiquitin. It is assumed that only the backbone is 13C,15N-labelled. Calculation time: min | 81 |
| `examples/nmr_proteins/hsqc_ubiquitin_a.m` | `hsqc_ubiquitin_a()` | 1H-15N HSQC of human ubiquitin, decoupling applied in both dimensions. Calculation time: hours, faster with a Tesla A100 | 72 |
| `examples/nmr_proteins/hsqc_ubiquitin_b.m` | `hsqc_ubiquitin_b()` | 1H-15N HSQC of human ubiquitin, without 1H decoupling in F1 and 15N decoupling in F2. Nitrogen-proton multiplicity is re | 73 |
| `examples/nmr_proteins/methyl_noesy_gb1.m` | `methyl_noesy_gb1()` | 1H-1H NOESY spectrum of GB1 with everything deuterated except methyl groups. Deuteria are kept in the spin system becaus | 84 |
| `examples/nmr_proteins/noesy_ubiquitin.m` | `noesy_ubiquitin()` | 1H-1H NOESY spectrum of ubiquitin with 65 ms mixing time. It is assumed that the protein is not 13C-or 15N-labelled. Cal | 82 |
| `examples/nmr_proteins/noesyhsqc_ubiquitin_deut.m` | `noesyhsqc_ubiquitin_deut()` | 1H-1H-15N NOESY-HSQC spectrum of 15N-labelled ubiquitin at 900 MHz with 90 ms mixing time. It is assumed that the protei | 101 |
| `examples/nmr_proteins/noesyhsqc_ubiquitin_prot.m` | `noesyhsqc_ubiquitin_prot()` | 1H-1H-15N NOESY-HSQC spectrum of 15N-labelled ubiquitin at 900 MHz with 65 ms mixing time. It is assumed that the protei | 94 |
| `examples/nmr_solids/case_studies/akbey_2h_13c_mas/fig2_three_site.m` | `fig2_three_site()` | Three-site position exchange for a deuterium nucleus. The sites differ in the chemical shift and the orientation of the  | 72 |
| `examples/nmr_solids/case_studies/akbey_2h_13c_mas/fig2_two_site.m` | `fig2_two_site()` | Two-site position exchange for a deuterium nucleus. The sites differ in the chemical shift and the orientation of the qu | 70 |
| `examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_adiabatic_vs_optimcon.m` | `cp_adiabatic_vs_optimcon()` | 1H-15N cross-polarisation experiment in the doubly rotating frame using (a) tangent-ramped adiabatic CP; (b) numerically | 147 |
| `examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_square_vs_ramp.m` | `cp_square_vs_ramp()` | 1H-15N cross-polarisation experiment in the doubly rotating frame using (a) fixed amplitude CP; (b) linearly ramped CP;  | 81 |
| `examples/nmr_solids/case_studies/magic_angle_calibration.m` | `magic_angle_calibration()` | Magic angle is usually calibrated using KBr powder. When the angle is not correctly set, the spinning sideband pat- tern | 84 |
| `examples/nmr_solids/case_studies/mas_powder_dd_nqi.m` | `mas_powder_dd_nqi()` | Powder magic angle spinning spectrum of a pair of dipole-coupled quadrupolar nuclei; this is apparently something that o | 67 |
| `examples/nmr_solids/case_studies/mathies_14n_13c/mas_powder_gly_13c.m` | `mas_powder_gly_13c()` | 13C MAS spectrum of glycine powder (assuming decoupling of 1H), computed using the Fokker-Planck MAS formalism and a sph | 103 |
| `examples/nmr_solids/case_studies/mathies_14n_13c/mas_powder_gly_14n.m` | `mas_powder_gly_14n()` | 14N MAS spectrum of glycine powder (assuming decoupling of 1H and 13C), computed using the Fokker-Planck MAS formalism a | 91 |
| `examples/nmr_solids/case_studies/mathies_carbonate/cp_mas_powder_mhc_fplanck_exchange.m` | `cp_mas_powder_mhc_fplanck_exchange()` | Cross-polarisation contact curve under magic angle spinning in the presence of chemical exchange for H1, H4 and C19 in t | 115 |
| `examples/nmr_solids/case_studies/mathies_carbonate/mas_powder_mhc_fplanck.m` | `mas_powder_mhc_fplanck()` | All protons in the unit cell of monohydrocalcite, magic angle spinning NMR simulation. Further details in: Calculation t | 88 |
| `examples/nmr_solids/case_studies/mathies_carbonate/mas_powder_mhc_fplanck_exchange.m` | `mas_powder_mhc_fplanck_exchange()` | Water protons in the unit cell of monohydrocalcite, inc- luding position exchange and MAS. Further details in: Calculati | 87 |
| `examples/nmr_solids/case_studies/mathies_carbonate/sle_nmr_dd_csa_mhc.m` | `sle_nmr_dd_csa_mhc()` | Water protons in the unit cell of monohydrocalcite, inc- luding slow isotropic rotational diffusion and MAS. Fur- ther d | 101 |
| `examples/nmr_solids/cp_acquire_mas_gly.m` | `cp_acquire_mas_gly()` | 1H-13C cross-polarisation followed by acquisition under magic angle spinning in alpha-glycine powder. Reduced Liouville  | 80 |
| `examples/nmr_solids/cp_contact_mas_nh_floquet.m` | `cp_contact_mas_nh_floquet()` | Cross-polarisation experiment in the doubly rotating frame. A single nitrogen-15 and a single proton. Spinning powder si | 60 |
| `examples/nmr_solids/cp_contact_mas_nh_fplanck.m` | `cp_contact_mas_nh_fplanck()` | Cross-polarisation experiment in the doubly rotating frame. A single nitrogen-15 and a single proton. Spinning powder si | 60 |
| `examples/nmr_solids/cp_contact_mas_nh_gridfree.m` | `cp_contact_mas_nh_gridfree()` | Cross-polarisation experiment in the doubly rotating frame. A single nitrogen-15 and a single proton. Spinning powder si | 64 |
| `examples/nmr_solids/cp_contact_mas_nhh.m` | `cp_contact_mas_nhh()` | Cross-polarisation experiment in the doubly rotating frame. A single nitrogen-15 in a bath of 8 protons scattered on a 2 | 77 |
| `examples/nmr_solids/cp_crystal_static_nh.m` | `cp_crystal_static_nh()` | 1H-15N cross-polarisation experiment in the doubly rotating frame. Static single crystal simulation. Calculation time: s | 51 |
| `examples/nmr_solids/cp_crystal_static_nhh.m` | `cp_crystal_static_nhh()` | Cross-polarisation experiment in the doubly rotating frame. A single nitrogen-15 in a bath of 8 protons scattered on a 2 | 67 |
| `examples/nmr_solids/cp_matching_1.m` | `cp_matching_1()` | Hartmann-Hahn matching condition test for a cross-polarisation experiment between a proton and a 15N nucleus under MAS.  | 68 |
| `examples/nmr_solids/cp_matching_2.m` | `cp_matching_2()` | Hartmann-Hahn matching condition test for a cross-polarisation experiment between a proton and a 15N nucleus. The experi | 70 |
| `examples/nmr_solids/cp_matching_3.m` | `cp_matching_3()` | Hartmann-Hahn matching condition test for a cross-polarisation experiment between a proton and a 15N nucleus. A 2D scan  | 81 |
| `examples/nmr_solids/cp_matching_4.m` | `cp_matching_4()` | Hartmann-Hahn matching condition test for a cross-polarisation experiment between a proton and a 15N nucleus in the pres | 77 |
| `examples/nmr_solids/cp_powder_static_nh.m` | `cp_powder_static_nh()` | 1H-15N cross-polarisation experiment in the doubly rotating frame. Static powder simulation. Calculation time: seconds | 51 |
| `examples/nmr_solids/cp_powder_static_nhh.m` | `cp_powder_static_nhh()` | Cross-polarisation experiment in the doubly rotating frame. A single nitrogen-15 in a bath of 8 protons scattered on a 2 | 67 |
| `examples/nmr_solids/cp_respiration.m` | `cp_respiration()` | 1H-13C RESPIRATION-CP experiment in the doubly rotating frame. Magic angle spinning simulation using Fokker-Planck forma | 57 |
| `examples/nmr_solids/dor_powder_nav_fplanck_freq.m` | `dor_powder_nav_fplanck_freq()` | Double angle spinning spectrum of N-acetylvaline 14N nucleus using 1D Fokker-Planck equation and a spherical grid. The c | 63 |
| `examples/nmr_solids/dor_powder_nav_fplanck_time.m` | `dor_powder_nav_fplanck_time()` | Double angle spinning spectrum of N-acetylvaline 14N nucleus using 1D Fokker-Planck equation and a spherical grid. The c | 66 |
| `examples/nmr_solids/fitting/bromide_csa_nqi/kbr_mas_fitting.m` | `kbr_mas_fitting()` | Fitting of a 79Br MAS NMR spectrum of potassium bromide with respect to the quadrupole coupling constant. The spectrum c | 113 |
| `examples/nmr_solids/fitting/vanadium_csa_nqi/vfit_4sim_ik.m` | `vfit_4sim_ik()` | Simultaneous fitting of multiple 51V MAS NMR spectra with respect to the chemical shielding anisotropy and quadrupole co | 168 |
| `examples/nmr_solids/fslghetcor_mas_gly.m` | `fslghetcor_mas_gly()` | FSLG-HETCOR of alpha-glycine powder under MAS. Calculation time: hours on NVidia Tesla A100, much longer on CPU | 91 |
| `examples/nmr_solids/hmqc_mas_dq.m` | `hmqc_mas_dq()` | Powder magic angle spinning CN2D experiment (rotor-synchronized de- tection) on a 14N-1H spin pair using 1D Fokker-Planc | 74 |
| `examples/nmr_solids/hmqc_mas_sq.m` | `hmqc_mas_sq()` | Powder magic angle spinning CN2D experiment (rotor-synchronized de- tection) on a 14N-1H spin pair using 1D Fokker-Planc | 74 |
| `examples/nmr_solids/mas_powder_ala_floquet.m` | `mas_powder_ala_floquet()` | 13C MAS spectrum of alanine powder (assuming decoupling of 1H), computed using the Floquet MAS formalism. Calculation ti | 60 |
| `examples/nmr_solids/mas_powder_ala_fplanck.m` | `mas_powder_ala_fplanck()` | 13C MAS spectrum of alanine powder (assuming decoupling of 1H), computed using the Fokker-Planck MAS formalism and a sph | 58 |
| `examples/nmr_solids/mas_powder_ala_gridfree.m` | `mas_powder_ala_gridfree()` | 13C MAS spectrum of alanine powder (assuming decoupling of 1H), computed using the grid-free Fokker-Planck MAS formalism | 63 |
| `examples/nmr_solids/mas_powder_csa_floquet.m` | `mas_powder_csa_floquet()` | Powder magic angle spinning spectrum of a single anisotropically shielded proton spin using a Floquet theory based forma | 56 |
| `examples/nmr_solids/mas_powder_csa_fplanck.m` | `mas_powder_csa_fplanck()` | Powder magic angle spinning spectrum of a pair of anisotropically shielded protons using a Fokker-Planck theory based fo | 51 |
| `examples/nmr_solids/mas_powder_csa_gridfree.m` | `mas_powder_csa_gridfree()` | Powder magic angle spinning spectrum of a pair of anisotropically shielded proton spins using grid- free Fokker-Planck e | 55 |
| `examples/nmr_solids/mas_powder_dip_floquet.m` | `mas_powder_dip_floquet()` | Spinning powder pulse-acquire experiment on a two-spin system with a dipolar coupling using Floquet theory. Calculation  | 56 |
| `examples/nmr_solids/mas_powder_dip_fplanck.m` | `mas_powder_dip_fplanck()` | Spinning powder pulse-acquire experiment on a two-spin system with a dipolar coupling using Fokker-Planck formalism: Cal | 54 |
| `examples/nmr_solids/mas_powder_dip_gridfree.m` | `mas_powder_dip_gridfree()` | Powder magic angle spinning spectrum of a pair of dipole-coupled proton spins using grid-free Fokker-Planck equation for | 54 |
| `examples/nmr_solids/mas_powder_gly_floquet.m` | `mas_powder_gly_floquet()` | 13C MAS spectrum of glycine powder (assuming decoupling of 1H), computed using the Floquet MAS formalism. Calculation ti | 61 |
| `examples/nmr_solids/mas_powder_gly_fplanck.m` | `mas_powder_gly_fplanck()` | 13C MAS spectrum of glycine powder (assuming decoupling of 1H), computed using the Fokker-Planck MAS formalism. Calculat | 57 |
| `examples/nmr_solids/mas_powder_gly_gridfree.m` | `mas_powder_gly_gridfree()` | 13C MAS spectrum of glycine powder (assuming decoupling of 1H), computed using the grid-free Fokker-Planck MAS formalism | 61 |
| `examples/nmr_solids/mas_powder_nqi_floquet.m` | `mas_powder_nqi_floquet()` | Powder magic angle spinning spectrum of a single quadrupolar deuterium nucleus using Floquet theory. Perturbative correc | 56 |
| `examples/nmr_solids/mas_powder_nqi_fplanck.m` | `mas_powder_nqi_fplanck()` | Powder magic angle spinning spectrum of a single quadrupolar deuterium nucleus using Fokker-Planck theory. Perturbative  | 56 |
| `examples/nmr_solids/mas_powder_nqi_gridfree.m` | `mas_powder_nqi_gridfree()` | Powder magic angle spinning spectrum of a single quadrupolar deuterium nucleus using grid-free Fokker-Planck MAS formali | 56 |
| `examples/nmr_solids/mas_powder_suc_floquet.m` | `mas_powder_suc_floquet()` | 13C MAS spectrum of sucrose powder (assuming decoupling of 1H), computed using the Floquet MAS formalism. Chemical shiel | 61 |
| `examples/nmr_solids/mas_powder_suc_fplanck.m` | `mas_powder_suc_fplanck()` | 13C MAS spectrum of sucrose powder (assuming decoupling of 1H), computed using the Fokker-Planck MAS formalism. Chemical | 58 |
| `examples/nmr_solids/mas_powder_suc_gridfree.m` | `mas_powder_suc_gridfree()` | 13C MAS spectrum of sucrose powder (assuming decoupling of 1H), computed using the grid-free Fokker-Planck MAS formalism | 66 |
| `examples/nmr_solids/mas_powder_trp_floquet.m` | `mas_powder_trp_floquet()` | 13C MAS spectrum of tryptophan powder (assuming decoupling of 1H), computed using the Floquet MAS formalism. Isotropic c | 109 |
| `examples/nmr_solids/mas_powder_trp_fplanck.m` | `mas_powder_trp_fplanck()` | 13C MAS spectrum of tryptophan powder (assuming decoupling of 1H), computed using the Fokker-Planck MAS formalism. Isotr | 106 |
| `examples/nmr_solids/mas_powder_trp_gridfree.m` | `mas_powder_trp_gridfree()` | 13C MAS spectrum of tryptophan powder (assuming decoupling of 1H), computed using the grid-free Fokker-Planck MAS formal | 118 |
| `examples/nmr_solids/mqmas_nqi.m` | `mqmas_nqi()` | Rotor-synchronous MQMAS spectrum of a 87Rb compound, transmitter set to the isotropic chemical shift. Calculation time:  | 59 |
| `examples/nmr_solids/pdsd_simple.m` | `pdsd_simple()` | 13C 2D PDSD spectrum of a simple test spin system. Calculation time: minutes, much faster on GPU. | 74 |
| `examples/nmr_solids/redor_curve.m` | `redor_curve()` | REDOR dephasing curve for a simple 13C-15N spin pair using Fokker-Planck MAS formalism. Calculation time: seconds | 54 |
| `examples/nmr_solids/rframe_nqi_dante.m` | `rframe_nqi_dante()` | DANTE MAS spectrum of a single quadrupolar 14N nucleus using 1D Fokker-Planck equation and a spherical grid. The calcula | 65 |
| `examples/nmr_solids/rframe_nqi_fplanck.m` | `rframe_nqi_fplanck()` | Powder magic angle spinning spectrum (rotor-synchronized detection) of a single quadrupolar 14N nucleus using 1D Fokker- | 57 |
| `examples/nmr_solids/static_powder_ala.m` | `static_powder_ala()` | 13C NMR spectrum of static alanine powder. Protons are assumed to be decoupled. Calculation time: seconds | 58 |
| `examples/nmr_solids/static_powder_csa.m` | `static_powder_csa()` | Two-spin static CSA powder pattern. Calculation time: seconds | 54 |
| `examples/nmr_solids/static_powder_dip.m` | `static_powder_dip()` | Static powder pulse-acquire experiment on a two-spin system with a dipolar coupling. Calculation time: seconds | 55 |
| `examples/nmr_solids/static_powder_gly.m` | `static_powder_gly()` | 13C NMR spectrum of glycine powder. Magnetic parameters loaded from a DFT simulation. Protons assumed to be decoupled. C | 61 |
| `examples/nmr_solids/static_powder_nqi_a.m` | `static_powder_nqi_a()` | Static quadrupolar 14N powder pattern of L-valyl-L-alanine using very large numerical orientation grid, set to reproduce | 54 |
| `examples/nmr_solids/static_powder_nqi_b.m` | `static_powder_nqi_b()` | Static powder 79Br NMR spectrum of potassium bromide. At least 3 quadrupolar tensors are necessary to reproduce the expe | 69 |
| `examples/nmr_solids/static_powder_suc.m` | `static_powder_suc()` | 13C NMR spectrum of static sucrose powder, protons are assumed to be decoupled. Calculation time: hours | 61 |
| `examples/nmr_solids/static_powder_trp.m` | `static_powder_trp()` | 13C NMR spectrum of tryptophan powder. Isotropic chemical shifts come from the experimental data. Coordinates and CSAs a | 73 |
| `examples/nmr_solids/wise_mas_gly.m` | `wise_mas_gly()` | WISE of alpha-glycine powder under MAS. Calculation time: hours, much faster on GPU | 82 |
| `examples/nmr_spen/conv_test.m` | `conv_test()` | Convergence and accuracy test of the spatial dynamics during a pulse field gradientspin (PFG) echo sequence. The accurac | 190 |
| `examples/nmr_spen/dosy_oneshot_1.m` | `dosy_oneshot_1()` | Oneshot DOSY pulse sequence for a system of three coupled spins with different relaxation rates. Timing: minutes on NVid | 88 |
| `examples/nmr_spen/idosyzs_test_1.m` | `idosyzs_test_1()` | Diffusion attenuation during soft pulses in a simplified model sequence of the Zangger-Sterk pure shift iDOSY with fitti | 104 |
| `examples/nmr_spen/psyche_rotenone.m` | `psyche_rotenone()` | PSYCHE pure-shift NMR spectrum of rotenone. Calculation time: hours, faster on a GPU. | 123 |
| `examples/nmr_spen/psyche_strychnine.m` | `psyche_strychnine()` | PSYCHE pure-shift NMR spectrum of strychnine. Calculation time: hours, faster on a GPU. | 101 |
| `examples/nmr_spen/psycosy_acrolein.m` | `psycosy_acrolein()` | PSYCOSY of Acrolein. Calculation time: hours, faster on a GPU. | 92 |
| `examples/nmr_spen/psycosy_dbpa.m` | `psycosy_dbpa()` | PSYCOSY of DBPA (dibromopropionic acid) ring. Calculation time: minutes on NVidia Tesla A100, much longer on CPU | 90 |
| `examples/nmr_spen/psycosy_salsalate.m` | `psycosy_salsalate()` | PSYCOSY of one salsalate ring. Calculation time: minutes on NVidia Tesla A100, much longer on CPU | 92 |
| `examples/nmr_spen/ufcosy_2spin.m` | `ufcosy_2spin()` | Ultrafast COSY for a coupled two-spin system. Calculation time: minutes on NVidia Tesla A100, much longer on CPU Jean-Ni | 75 |
| `examples/nmr_spen/ufcosy_3spin.m` | `ufcosy_3spin()` | Ultrafast COSY for a coupled three-spin system. Calculation time: minutes on NVidia Tesla A100, much longer on CPU Jean- | 77 |
| `examples/nmr_spen/ufdosy_1spin.m` | `ufdosy_1spin()` | Ultrafast DOSY for one spin. Calculation time: seconds on NVidia Tesla A100, much longer on CPU Ludmilla Guduff Jean-Nic | 87 |
| `examples/nmr_spen/ufdosy_2spin.m` | `ufdosy_2spin()` | Ultrafast DOSY for two coupled spins with additional complications like DD and CSA relaxation, and spatial flow. Calcula | 102 |
| `examples/nmr_spen/ufdosycosy_2spin.m` | `ufdosycosy_2spin()` | Ultrafast 3D DOSY-COSY for a two-spin system. Calculation time: hours on NVidia Tesla A100, much longer on CPU Ludmilla  | 100 |
| `examples/nmr_spen/ufmq_2spin.m` | `ufmq_2spin()` | 2Q ultrafast MaxQ NMR spectrum for a coupled two-spin system in the presence of realistic diffusion. Calculation time: m | 106 |
| `examples/nmr_spen/ufmq_4spin.m` | `ufmq_4spin()` | 4Q ultrafast MaxQ NMR spectrum for a coupled four-spin system in the presence of realistic diffusion. Calculation time:  | 114 |
| `examples/nmr_spen/ufmq_6spin.m` | `ufmq_6spin()` | 6Q ultrafast MaxQ NMR spectrum for a coupled six-spin system in the presence of realistic diffusion. Calculation time: h | 125 |
| `examples/nmr_stochastic/snmr_gb1.m` | `snmr_gb1()` | A Primas-style stochastic NMR experiment on GB1 protein. The calculation requires a terabyte of RAM and NVidia A100 GPU. | 165 |
| `examples/nmr_stochastic/snmr_strychnine.m` | `snmr_strychnine()` | A Primas-style stochastic NMR experiment on strychnine. The calculation requires an NVidia Titan V GPU at a minimum. Cal | 107 |
| `examples/nmr_zerofield/earth_field_dfp.m` | `earth_field_dfp()` | Earth's field NMR Simulation for 2,6-difluoropyridine; replicates simulated spectra in Figure 7 of without the weighted  | 102 |
| `examples/nmr_zerofield/field_drop_acetonitrile.m` | `field_drop_acetonitrile()` | Zero-field NMR spectroscopy -acetonitrile. The simulation proceeds by computing the exact thermal equilibrium state and  | 76 |
| `examples/nmr_zerofield/small_field_acetonitrile.m` | `small_field_acetonitrile()` | Small-field NMR spectroscopy -acetonitrile with 13C on the methyl group. Set to reproduce Figure 3 from Calculation time | 59 |
| `examples/nmr_zerofield/small_field_formic_acid.m` | `small_field_formic_acid()` | Zero-field NMR spectroscopy -15N pyridine. Set to reproduce Fig 2 from http://dx.doi.org/10.1103/PhysRevLett.107.107601  | 52 |
| `examples/nmr_zerofield/zero_field_benzene.m` | `zero_field_benzene()` | Zero-field NMR spectroscopy -benzene with one 13C nucleus. Set to reproduce Figure 2 from Calculation time: seconds | 74 |
| `examples/nmr_zerofield/zero_field_methanol.m` | `zero_field_methanol()` | Zero-field NMR spectroscopy -13C methanol. Set to reproduce Figure 1 from http://dx.doi.org/10.1016/j.cplett.2013.06.042 | 55 |
| `examples/nmr_zerofield/zero_field_pyridine.m` | `zero_field_pyridine()` | Zero-field NMR spectroscopy -15N pyridine. Set to reproduce Figure 3 from http://dx.doi.org/10.1021/ja2112405 Calculatio | 66 |
| `examples/nqr/nutation_nqr_iodine.m` | `nutation_nqr_iodine()` | Powder NQR nutation curve for a system with a single 127I nucleus. Calculation time: seconds | 69 |
| `examples/nqr/pure_nqr_iodine.m` | `pure_nqr_iodine()` | Powder NQR spectrum of a system with a single 127I nucleus. Calculation time: seconds | 53 |
| `examples/nqr/pure_nqr_nitrogen.m` | `pure_nqr_nitrogen()` | Powder NQR spectrum of a system with a single 14N nucleus. Calculation time: seconds | 54 |
| `examples/optimal_control/bloch_siegert/bloch_siegert_a.m` | `bloch_siegert_a()` | Bloch-Siegert shift compensation functionality demo. The script optimises a 90-degree pulse (Lz -> Lx) for a sing- le sp | 106 |
| `examples/optimal_control/bloch_siegert/bloch_siegert_b.m` | `bloch_siegert_b()` | Bloch-Siegert shift compensation functionality demo. The script optimises a universal rotation pulse for a range of reso | 108 |
| `examples/optimal_control/bloch_siegert/coote_badcop.m` | `coote_badcop()` | Reproduction of BADCOP-style selective decoupling logic from Coote et al. with Bloch-Siegert corrections enabled in the  | 175 |
| `examples/optimal_control/bloch_siegert/coote_goodcop.m` | `coote_goodcop()` | Reproduction of the GOODCOP pulse design logic from Coote et al. with Bloch-Siegert corrections enabled in the optimiser | 147 |
| `examples/optimal_control/bloch_siegert/ramsey_shifts.m` | `ramsey_shifts()` | Ramsey shifts of other spins under an off-resonant drive. A proton channel drives a three-spin system in which 13C and 1 | 107 |
| `examples/optimal_control/bloch_siegert/yusuke_14n_broadening_demo.m` | `yusuke_14n_broadening_demo()` | Reduced effective-model illustration of the trade-off discussed by Nehra, Agarwal, and Nishiyama for 14N decoupling unde | 205 |
| `examples/optimal_control/bloch_siegert/yusuke_1h_14n_optimal_vs_cw_demo.m` | `yusuke_1h_14n_optimal_vs_cw_demo()` | Reduced heteronuclear control demonstration inspired by the low-power offset-tolerant 14N decoupling papers of Nehra, Ag | 264 |
| `examples/optimal_control/bloch_siegert/yusuke_optimal_vs_cw_demo.m` | `yusuke_optimal_vs_cw_demo()` | Bloch-Siegert-aware phase optimisation compared to a simple constant- phase low-power cycle. This is the control-side co | 240 |
| `examples/optimal_control/case_studies/Tosner_JMR_2009/bb_inversion_pulse.m` | `bb_inversion_pulse()` | Broadband inversion pulse design for liquid-state NMR. Reprodu- ces, using Spinach, the second example from: A single pr | 101 |
| `examples/optimal_control/case_studies/Tosner_JMR_2009/bb_refocusing_pulse.m` | `bb_refocusing_pulse()` | Spinach implementation of the broadband refocusing example from GRAPE is used to design a 200 µs broadband x-phase π pul | 97 |
| `examples/optimal_control/case_studies/Tosner_JMR_2009/coherence_transfer.m` | `coherence_transfer()` | The first optimal control example from A heteronuclear two-spin system (1H–13C) with an scalar and both nuclei set on re | 79 |
| `examples/optimal_control/distortions/distortions_figure_1.m` | `distortions_figure_1()` | Figure 1 from the paper by Rasulov and Kuprov: | 73 |
| `examples/optimal_control/distortions/distortions_figure_2.m` | `distortions_figure_2()` | Figure 2 from the paper by Rasulov and Kuprov: | 39 |
| `examples/optimal_control/distortions/distortions_figure_3.m` | `distortions_figure_3()` | Figure 3 from the paper by Rasulov and Kuprov: | 225 |
| `examples/optimal_control/distortions/distortions_figure_4_bot.m` | `distortions_figure_4_bot()` | Figure 4 (bottom) from the paper by Rasulov and Kuprov: | 132 |
| `examples/optimal_control/distortions/distortions_figure_4_top.m` | `distortions_figure_4_top()` | Figure 4 (top) from the paper by Rasulov and Kuprov: | 138 |
| `examples/optimal_control/distortions/kernel_estimation/kernel_application.m` | `kernel_application()` | HiPER instrument filter function kernel application to a complicated shaped pulse and a comparison with expe- rimental m | 44 |
| `examples/optimal_control/distortions/kernel_estimation/kernel_from_87rb.m` | `kernel_from_87rb()` | Transmitter and probe distortion kernel of a 400 MHz Bruker spectrometer fitted with a 4 mm Phoenix MAS probe, estimated | 130 |
| `examples/optimal_control/distortions/kernel_estimation/kernel_from_antenna.m` | `kernel_from_antenna()` | HiPER instrument filter function kernel estimation from the quadrature components recorded by an antenna placed close to | 61 |
| `examples/optimal_control/distortions/kernel_estimation/kernel_from_transm.m` | `kernel_from_transm()` | Response function extraction from the transmission profile of the HiPER instrument. We are sending a linear chirp via th | 57 |
| `examples/optimal_control/distortions/restrans_test.m` | `restrans_test()` | Resonator transform test. Sends a square pulse into a simple resonator model and plots the time-domain response. | 38 |
| `examples/optimal_control/distortions/rlc_response_1.m` | `rlc_response_1(interp_type)` | An illustration of the effect of the resonator response function on a typical composite pulse in NMR spectroscopy. The a | 80 |
| `examples/optimal_control/distortions/rlc_response_2.m` | `rlc_response_2()` | Probe circuit response effect on the accuracy of the deu- terium pre-phasing pulse designed to set deuterium magne- tisa | 102 |
| `examples/optimal_control/distortions/rlc_response_3.m` | `rlc_response_3()` | Probe circuit response effect on the accuracy of the deu- terium pre-phasing pulse designed to set deuterium magne- tisa | 102 |
| `examples/optimal_control/distortions/scope_readout/pulse_heterodyne.m` | `pulse_heterodyne()` | Digital processing for a recording of the proton optimal control pulse, done on a 1 GHz oscilloscope using the 13C coil  | 75 |
| `examples/optimal_control/features_ampl.m` | `features_ampl()` | An illustration of amplitude profiling in a phase-modulated pulse optimisation. The amplitude profile is supplied by the | 115 |
| `examples/optimal_control/features_bss.m` | `features_bss()` | Optimal control pulse optimisation with Bloch-Siegert shift corrections switched on. A single proton with a Larmor frequ | 90 |
| `examples/optimal_control/features_curv.m` | `features_curv()` | A transfer of coherence from longitudinal magnetization into a two-spin singlet state with a distribution of B1 powers.  | 105 |
| `examples/optimal_control/features_diss_drift.m` | `features_diss_drift()` | Optimal control optimisation of a pulse performing magnetisa- tion transfer from H(N) to C(O) in a typical protein backb | 135 |
| `examples/optimal_control/features_dt_opt.m` | `features_dt_opt()` | Optimisation of slice durations in a composite inversion pulse with specified amplitudes, phases, and a constrain- ed ov | 110 |
| `examples/optimal_control/features_dt_var.m` | `features_dt_var()` | Optimal control pulse optimisation for state-to-state transfer across a scalar coupling in a hydrofluorocarbon fragment  | 94 |
| `examples/optimal_control/features_freeze.m` | `features_freeze()` | Optimal control pulse optimisation for state-to-state transfer across a scalar coupling in a hydrofluorocarbon fragment  | 114 |
| `examples/optimal_control/features_keyhole.m` | `features_keyhole()` | Optimal control pulse optimisation for state-to-state transfer across a scalar coupling in a hydrofluorocarbon fragment  | 102 |
| `examples/optimal_control/features_multitarget.m` | `features_multitarget()` | An example of multi-target optimal control pulse design in the context of singlet state NMR spectroscopy. A pulse is des | 110 |
| `examples/optimal_control/features_newton.m` | `features_newton()` | Optimal control pulse optimisation for state-to-state transfer across two scalar couplings in a hydrofluorocarbon fragme | 106 |
| `examples/optimal_control/features_phase_cycle.m` | `features_phase_cycle()` | Optimal control pulse optimisation for state-to-state transfer across a scalar coupling in a hydrofluorocarbon fragment  | 111 |
| `examples/optimal_control/features_trapezium.m` | `features_trapezium()` | Optimal control pulse optimisation for state-to-state transfer across a scalar coupling in a hydrofluorocarbon fragment  | 113 |
| `examples/optimal_control/features_wave_basis.m` | `features_wave_basis()` | An illustration of basis set coefficient optimisation for a pul- se optimised as a linear combination of user-specified  | 116 |
| `examples/optimal_control/magic_pulse_cart.m` | `magic_pulse_cart()` | A template file for the "magic pulse" optimisations. The term refers to a family of broadband NMR pulses that are tolera | 135 |
| `examples/optimal_control/magic_pulse_phase.m` | `magic_pulse_phase()` | A template file for the "magic pulse" optimisations. The term refers to a family of broadband NMR pulses that are tolera | 135 |
| `examples/optimal_control/mas_powder_control.m` | `mas_powder_control()` | Optimal control pulse starting with Lz and populating the Ly state on 87Rb in a quadrupolar rubidium system under magic  | 106 |
| `examples/optimal_control/pattern_pulse_1.m` | `pattern_pulse_1()` | Nutation frequency selective excitation described in Glaser group paper (https://doi.org/10.1016/j.jmr.2004.12.005). Use | 103 |
| `examples/optimal_control/pattern_pulse_2.m` | `pattern_pulse_2()` | Transmitter offset selective excitation described in Glaser group paper (https://doi.org/10.1016/j.jmr.2004.12.005). Use | 106 |
| `examples/optimal_control/pulse_analysis.m` | `pulse_analysis()` | An example of spectrogram analysis for a quadratic chirp pulse; adapted from Matlab example set. Calculation time: secon | 23 |
| `examples/optimal_control/state_transfer_coop.m` | `state_transfer_coop()` | Optimal control pulse optimisation for state-to-state transfer in a quadrupolar 14N spin at a fixed orientation and powe | 89 |
| `examples/optimal_control/state_transfer_m2s.m` | `state_transfer_m2s()` | A transfer of coherence from longitudinal magnetization into a two-spin singlet state in allyl pyruvate with a distribut | 147 |
| `examples/optimal_control/state_transfer_pro.m` | `state_transfer_pro()` | Optimal control optimisation of a pulse performing magnetisa- tion transfer from H(N) to C(O) in a typical protein backb | 122 |
| `examples/optimal_control/state_transfer_s2m.m` | `state_transfer_s2m()` | A transfer of coherence from a two-proton singlet state to a nearby carbon in a setting typically encountered in parahyd | 100 |
| `examples/optimal_control/state_transfer_wf.m` | `state_transfer_wf()` | A transfer of population from the lowermost energy level in a four-spin system to the uppermost energy level using the w | 92 |
| `examples/optimal_control/static_powder_control.m` | `static_powder_control()` | Optimal control optimisation for a pulse that is designed to set deuterium magnetisation in a -CD3 group of alanine up f | 138 |
| `examples/optimal_control/steady_orbit/solid_effect_chirp.m` | `solid_effect_chirp()` | Panoramic optimisation for stroboscopic steady state DNP with the timing and power settings matching the XiX ex- eriment | 130 |
| `examples/optimal_control/steady_orbit/solid_effect_dq.m` | `solid_effect_dq()` | Panoramic optimisation for stroboscopic steady state DNP with the timing and power settings matching the XiX ex- eriment | 127 |
| `examples/optimal_control/steady_orbit/solid_effect_int.m` | `solid_effect_int()` | Panoramic optimisation for stroboscopic steady state DNP with the timing and power settings matching the XiX ex- eriment | 127 |
| `examples/optimal_control/steady_orbit/solid_effect_xix.m` | `solid_effect_xix()` | Panoramic optimisation for stroboscopic steady state DNP with the timing and power settings matching the XiX ex- eriment | 131 |
| `examples/parahydrogen/altadena_propanal.m` | `altadena_propanal()` | ALTADENA experiment simulation for the parahydrogenation of acrolein into propanal. Simple model of the ALTADENA effect  | 68 |
| `examples/parahydrogen/case_studies/hyperpolarised_deuterium/bubble_pulse_acquire.m` | `bubble_pulse_acquire()` | Simulated PNL (partially negative line) spectrum of ortho-deuterium in the presence of a parahydrogena- tion catalyst. B | 120 |
| `examples/parahydrogen/case_studies/hyperpolarised_deuterium/just_bubbling.m` | `just_bubbling()` | Evolution of state populations under ortho-deuterium bubbling in the presence of a parahydrogenation cata- lyst. No puls | 115 |
| `examples/parahydrogen/case_studies/hyperpolarised_deuterium/kinetic_isotope_effect.m` | `kinetic_isotope_effect()` | Evolution of state populations under ortho-deuterium bubbling in the presence of a parahydrogenation cata- lyst. Bubblin | 133 |
| `examples/parahydrogen/case_studies/spontaneous_singlet_to_z/rlx_trajectory.m` | `rlx_trajectory()` | Time dependence of LzSz and Lz+Sz spin orders in a para- hydrogen molecule coordinated to a nickel cage that cre- ates l | 75 |
| `examples/parahydrogen/ortho_deuterium.m` | `ortho_deuterium()` | Ortho-deuteration simulation for acrylonitrile in Figure 1 of the paper by Natterer, Greve, and Bargon: Simulation time: | 68 |
| `examples/parahydrogen/pasadena_ethylbenzene.m` | `pasadena_ethylbenzene()` | PASADENA experiment simulation for the parahydrogenation of styrene into ethylbenzene. Set to reproduce the top trace of | 70 |
| `examples/parahydrogen/pasadena_propanal.m` | `pasadena_propanal()` | PASADENA experiment simulation for the parahydrogenation of acrolein into propanal. Calculation time: seconds | 64 |
| `examples/parahydrogen/sabre_pyridine.m` | `sabre_pyridine()` | SABRE experiment simulation for Eibe Duecker and Christian Griesinger. Set to reproduce Figure 3b from http://dx.doi.org | 114 |
| `examples/quantum_tech/circuit_qed/cavity_fock_grape_a.m` | `cavity_fock_grape_a()` | GRAPE preparation of a cavity Fock state through a dispersively coupled qubit, using piecewise-constant drives on both t | 127 |
| `examples/quantum_tech/circuit_qed/cavity_fock_grape_b.m` | `cavity_fock_grape_b()` | GRAPE preparation of a cavity Fock state through a dispersively coupled qubit using smooth band-limited drives. The cont | 113 |
| `examples/quantum_tech/circuit_qed/cross_resonance.m` | `cross_resonance()` | Cross-resonance gate mechanism between two fixed-frequency transmons in the laboratory frame. The control transmon is dr | 143 |
| `examples/quantum_tech/circuit_qed/resonator_decay.m` | `resonator_decay()` | Open-system dynamics of a leaky microwave resonator at finite temperature. A Fock state decays as a downward cascade thr | 117 |
| `examples/quantum_tech/circuit_qed/transmon_drag.m` | `transmon_drag()` | DRAG correction of a resonant Gaussian pulse on a three-level Duffing transmon in the laboratory frame. The derivative o | 129 |
| `examples/quantum_tech/circuit_qed/transmon_two_photon.m` | `transmon_two_photon()` | Two-photon transition in a four-level Duffing transmon. The drive carrier is placed halfway between the 0-1 and 1-2 tran | 95 |
| `examples/quantum_tech/diamond_defects/diamond_co_epr_xw.m` | `diamond_co_epr_xw()` | Field-swept powder EPR spectra of a Co centre in diamond at X and W bands. Calculation time: minutes. | 67 |
| `examples/quantum_tech/diamond_defects/diamond_gev0_epr_xw.m` | `diamond_gev0_epr_xw()` | Field-swept powder EPR spectra of GeV0 centre in diamond at X and W bands. Calculation time: seconds. | 67 |
| `examples/quantum_tech/diamond_defects/diamond_n2vm_epr_xw.m` | `diamond_n2vm_epr_xw()` | Field-swept powder EPR spectra of an N2V-centre in diamond at X and W bands. Calculation time: seconds. | 68 |
| `examples/quantum_tech/diamond_defects/diamond_ni_epr_xw.m` | `diamond_ni_epr_xw()` | Field-swept powder EPR spectra of Ni defects in diamond at X and W bands. Calculation time: seconds. | 68 |
| `examples/quantum_tech/diamond_defects/diamond_ninter_epr_xw.m` | `diamond_ninter_epr_xw()` | Field-swept powder EPR spectra of a nitrogen interstitial defect in diamond at X and W bands. Calculation time: seconds. | 68 |
| `examples/quantum_tech/diamond_defects/diamond_nvm_epr_xw.m` | `diamond_nvm_epr_xw()` | Field-swept powder EPR spectra of an NV centre in diamond at X and W bands. Calculation time: seconds. | 67 |
| `examples/quantum_tech/diamond_defects/diamond_p1_13c_epr_xw.m` | `diamond_p1_13c_epr_xw()` | Field-swept powder EPR spectra of a P1 centre in 13C-enriched diamond at X and W bands. Calculation time: minutes. | 71 |
| `examples/quantum_tech/diamond_defects/diamond_p1_epr_xw.m` | `diamond_p1_epr_xw()` | Field-swept powder EPR spectra of a P1 centre in diamond at X and W bands. Calculation time: seconds. | 67 |
| `examples/quantum_tech/diamond_defects/diamond_siv0_epr_xw.m` | `diamond_siv0_epr_xw()` | Field-swept powder EPR spectra of SiV0 centre in diamond at X and W bands. Calculation time: seconds. | 68 |
| `examples/quantum_tech/jaynes_cummings_a.m` | `jaynes_cummings_a()` | Jaynes-Cummings coupling between a spin and an electromagnetic cavity mode with five population numbers included. The av | 70 |
| `examples/quantum_tech/jaynes_cummings_b.m` | `jaynes_cummings_b()` | Jaynes-Cummings coupling between a spin and an electromagnetic cavity mode with five population numbers included. A time | 79 |
| `examples/quantum_tech/jaynes_cummings_c.m` | `jaynes_cummings_c()` | An exchange-coupled two-electron system with the electrons having independent Jaynes-Cummings couplings to the same mode | 74 |
| `examples/quantum_tech/optomechanics/optomech_sideband.m` | `optomech_sideband()` | Optomechanical sideband transfer of a phonon Fock state into a driven cavity. A red-detuned coherent drive on the cavity | 89 |
| `examples/quantum_tech/spin_cavity_purcell_effect.m` | `spin_cavity_purcell_effect()` | Cavity-induced spin relaxation in the EPR Purcell regime. Coherent Jaynes-Cummings exchange is combined with rapid cavit | 124 |
| `examples/quantum_tech/spin_cavity_vacuum_rabi.m` | `spin_cavity_vacuum_rabi()` | Vacuum Rabi oscillation between an electron spin and a micro- wave cavity mode in the Jaynes-Cummings approximation. Thi | 65 |
| `examples/quantum_tech/spin_phonon_avoided_crossing.m` | `spin_phonon_avoided_crossing()` | Avoided crossing between an electron spin transition and a quantised phonon mode in the resonant spin-phonon exchange mo | 78 |
| `examples/quantum_tech/spin_phonon_dephasing.m` | `spin_phonon_dephasing()` | Longitudinal spin-phonon coupling producing spin coherence modulation and spin-conditioned displacement of a quantised v | 72 |
| `examples/quantum_tech/spin_phonon_swap.m` | `spin_phonon_swap()` | Resonant excitation swap between an electron spin and a quantised phonon mode. The model is the spin-phonon Jaynes- Cumm | 64 |
| `examples/quantum_tech/tavis_cummings_splitting.m` | `tavis_cummings_splitting()` | Collective normal-mode splitting in the Tavis-Cummings model for one to four identical electron spins coupled to a commo | 86 |
| `examples/quantum_tech/transmon_cavity_swap.m` | `transmon_cavity_swap()` | Vacuum Rabi swap between a transmon and a microwave cavity mode, both represented by truncated bosonic Weyl algebras. Th | 65 |
| `examples/quantum_tech/transmon_duffing_ladder.m` | `transmon_duffing_ladder()` | Duffing-model energy ladder of a weakly anharmonic transmon, showing how the transition frequencies separate as anharmo- | 62 |
| `examples/quantum_tech/transmon_frog.m` | `transmon_frog()` | Basic implementation of a Frequency Robust Gate (FROG) for a single transmon, based on: Ensemble GRAPE optimisation with | 85 |
| `examples/quantum_tech/transmon_rabi_leakage.m` | `transmon_rabi_leakage()` | Rabi dynamics of a driven four-level transmon in the Duffing approximation, including leakage into the second and third  | 61 |
| `examples/quantum_tech/transmon_ramsey_chevron.m` | `transmon_ramsey_chevron()` | Ramsey chevron of a three-level transmon in the Duffing ap- proximation. A nominal pi/2 pulse prepares a coherence, and  | 86 |
| `examples/quantum_tech/transmon_stirap.m` | `transmon_stirap()` | Basic single transmon system with Duffing model in- teractions, parameters and model from: Ensemble GRAPE optimisation w | 85 |
| `examples/quantum_tech/transmon_transfer.m` | `transmon_transfer()` | Basic two-transmon system with Duffing model interacti- ons and a flip-flop coupling; coherence transfer from transmon 1 | 95 |
| `examples/relaxation_theory/aniso_diff_test_1.m` | `aniso_diff_test_1()` | Relaxation superoperator calculation for a dipole-coupled two-spin system with an anisotropic rotational diffusion tenso | 53 |
| `examples/relaxation_theory/aniso_diff_test_2.m` | `aniso_diff_test_2()` | Relaxation superoperator calculation for an anisotropically shielded two-spin system with an anisotropic rotational diff | 48 |
| `examples/relaxation_theory/cpmg_echo_train.m` | `cpmg_echo_train()` | CPMG echo train in a powder. Calculation time: seconds | 55 |
| `examples/relaxation_theory/csa_antisymm_1.m` | `csa_antisymm_1()` | Longitudinal and transverse relaxation rates in a system with a significant antisymmetry in the shielding tensor. Calcul | 49 |
| `examples/relaxation_theory/csa_csa_xcorr_1.m` | `csa_csa_xcorr_1()` | Complete Bloch-Redfield-Wangsness relaxation superoperator in a system with two anisotropically shielded nuclei. Spinach | 38 |
| `examples/relaxation_theory/csa_csa_xcorr_2.m` | `csa_csa_xcorr_2()` | CSA-CSA cross-correlation in the 103Rh subsystem and its effect on the widths of the three lines of the proton triplet.  | 63 |
| `examples/relaxation_theory/dd_csa_xcorr_1.m` | `dd_csa_xcorr_1()` | Complete Bloch-Redfield-Wangsness relaxation superoperator in a system with two anisotropically shielded nuclei with a d | 44 |
| `examples/relaxation_theory/dd_csa_xcorr_2.m` | `dd_csa_xcorr_2()` | DD-CSA cross-correlation -a reproduction of Fig 5a from the paper by Grace and Kumar (http://dx.doi.org/10.1006/jmra.199 | 103 |
| `examples/relaxation_theory/dd_quad_xcorr_1.m` | `dd_quad_xcorr_1()` | Complete Bloch-Redfield-Wangsness relaxation superoperator in a system with a quadrupolar coupling and a dipole coupling | 60 |
| `examples/relaxation_theory/dd_relaxation_1.m` | `dd_relaxation_1()` | Complete Bloch-Redfield-Wangsness relaxation superoperator in a system with a dipolar coupling between two spins. Dipola | 73 |
| `examples/relaxation_theory/dd_relaxation_2.m` | `dd_relaxation_2()` | Complete Bloch-Redfield-Wangsness relaxation superoperator in a system with dipolar coupling between spins. The dipolar  | 42 |
| `examples/relaxation_theory/dd_relaxation_3.m` | `dd_relaxation_3()` | Extreme narrowing limit case comparison between the dipolar relaxation rates in proton-proton and proton-deuterium syste | 76 |
| `examples/relaxation_theory/from_md/ngce_test.m` | `ngce_test()` | Test of the numerical integral route to the Redfield relaxation superoperator against the analytical results for isotrop | 77 |
| `examples/relaxation_theory/from_md/rwalk_plot.m` | `rwalk_plot()` | A plot of a typical random walk on a sphere. | 33 |
| `examples/relaxation_theory/from_md/sucrose_eight_spins.m` | `sucrose_eight_spins()` | One of the calculations reported in the JMR paper with Jim Prestegard: an eight-spin subsystem from the glucose ring of  | 123 |
| `examples/relaxation_theory/from_md/sucrose_three_spins.m` | `sucrose_three_spins()` | One of the calculations reported in the JMR paper with Jim Prestegard: a three-spin subsystem from the glucose ring of s | 112 |
| `examples/relaxation_theory/hfc_antisymm_1.m` | `hfc_antisymm_1()` | Longitudinal and transverse relaxation rates in a system with a significant antisymmetry in the hyperfine tensor. Calcul | 73 |
| `examples/relaxation_theory/hfc_relaxation_1.m` | `hfc_relaxation_1()` | Computes and prints the full Redfield superoperator for an electron- nucleus system with an anisotropic hyperfine coupli | 74 |
| `examples/relaxation_theory/inv_rec_1.m` | `inv_rec_1()` | A simple inversion-recovery experiment; longitudinal magnetisation is monitored as a function of time. Calculation time: | 63 |
| `examples/relaxation_theory/inv_rec_2.m` | `inv_rec_2()` | An example of inversion recovery experiment simulation for a strychnine spin system. Calculation time: minutes. | 97 |
| `examples/relaxation_theory/maz_noesy_1.m` | `maz_noesy_1()` | 15N-labelled methylaziridine NOESY, including the effects of the scalar relaxation of the first kind, caused by the modu | 137 |
| `examples/relaxation_theory/maz_noesy_2.m` | `maz_noesy_2()` | Methylaziridine NOESY, including the effects of the scalar relaxation of the first kind (caused by the modulation of J-c | 145 |
| `examples/relaxation_theory/maz_pulse_acquire.m` | `maz_pulse_acquire()` | Methylaziridine pulse-acquire, showing the effect of the scalar relaxation of the second kind due to the fast quadrupola | 132 |
| `examples/relaxation_theory/nz_vs_redfield_1.m` | `nz_vs_redfield_1()` | Nakajima-Zwanzig relaxation theory against Redfield theory for a two-spin system with dipolar and CSA cross-correlations | 98 |
| `examples/relaxation_theory/quad_relaxation_1.m` | `quad_relaxation_1()` | 14N quadrupolar relaxation in glycine in liquid state. The numerical output of Spinach is compared to the analytical equ | 51 |
| `examples/relaxation_theory/quad_scalar_1.m` | `quad_scalar_1()` | NMR spectrum of 17O enriched water inside a fullerene cage. A rather exotic combination of quadrupolar relaxation on the | 57 |
| `examples/relaxation_theory/sat_rec_1.m` | `sat_rec_1()` | A simple saturation-recovery experiment. Calculation time: seconds. | 61 |
| `examples/relaxation_theory/scalar_relaxation_1.m` | `scalar_relaxation_1()` | Redfield superoperator for the scalar relaxation of the first kind in a two-proton system with a noisy J-coupling. This  | 43 |
| `examples/relaxation_theory/sle_esr_nitroxide_1.m` | `sle_esr_nitroxide_1()` | Comparison between nitroxide simulation using SLE formalism and Redfield relaxation theory. Calculation time: seconds | 79 |
| `examples/relaxation_theory/sle_esr_nitroxide_2.m` | `sle_esr_nitroxide_2()` | Slow motion regime simulation of an ESR spectrum of a nitroxide radical. Set to reproduce Figure 2 from the paper by Con | 60 |
| `examples/relaxation_theory/sle_nmr_dd_csa.m` | `sle_nmr_dd_csa()` | 15N-1H DD-CSA cross-correlation in a protein amide bond spin system using SLE formalism and Redfield relaxation theory.  | 82 |
| `examples/relaxation_theory/sle_solid_limit.m` | `sle_solid_limit()` | Solid limit of Stochastic Liouville equation formalism. Calculation time: hours | 73 |
| `examples/relaxation_theory/t1t2_strychnine.m` | `t1t2_strychnine()` | Relaxation analysis for strychnine, dipolar processes only. Calculation time: seconds. | 39 |
| `examples/relaxation_theory/trosy_double.m` | `trosy_double()` | Hari Arthanari's Double TROSY effect. Calculation time: seconds. | 85 |
| `examples/relaxation_theory/trosy_fluorine_num.m` | `trosy_fluorine_num()` | Transverse relaxation rate as a function of the applied magnetic field in a 3-fluorotyrosine labelled protein. The fluor | 102 |
| `examples/relaxation_theory/trosy_fluorine_sym.m` | `trosy_fluorine_sym()` | Transverse relaxation rate as a function of the applied magnetic field in a 3-fluorotyrosine labelled protein. The fluor | 74 |
| `examples/relaxation_theory/trosy_methyl.m` | `trosy_methyl()` | Methyl trosy in a rapidly rotating 13CH3 group of a slowly tumbling protein, simulated using the Fokker-Planck formalism | 155 |
| `examples/relaxation_theory/trosy_nh.m` | `trosy_nh()` | Transverse relaxation rate as a function of the applied magnetic field at a typical amide N-H group in a protein. Rotati | 100 |
| `examples/relaxation_theory/trosy_proton.m` | `trosy_proton()` | Transverse relaxation rate as a function of the applied magnetic field at the C-H group in position 3 of the aromatic ri | 101 |
| `examples/relaxation_theory/trosy_selenium.m` | `trosy_selenium()` | Transverse relaxation rate as a function of the applied magnetic field in ethylselenol. The selenium atom and its direct | 102 |
| `examples/shaped_pulses/shaped_pulse_chirp_af.m` | `shaped_pulse_chirp_af()` | Chirp inversion pulse using the Fokker-Planck formalism. Fewer points are required by the amplitude-frequency method tha | 87 |
| `examples/shaped_pulses/shaped_pulse_chirp_xy.m` | `shaped_pulse_chirp_xy()` | Chirped inversion pulse. Calculation time: seconds | 86 |
| `examples/shaped_pulses/shaped_pulse_fp.m` | `shaped_pulse_fp()` | An off-resonance rectangular soft pulse simulated using the Fokker-Planck formalism. Note that the pulse frequency off-  | 79 |
| `examples/shaped_pulses/shaped_pulse_gaussian.m` | `shaped_pulse_gaussian()` | Gaussian 90-degree pulse on a chain of 31 strongly coupled protons. Calculation time: seconds | 91 |
| `examples/shaped_pulses/shaped_pulse_q5.m` | `shaped_pulse_q5()` | 90-degree Q5 pulse on a chain of 31 strongly coupled protons. Calculation time: seconds | 91 |
| `examples/shaped_pulses/shaped_pulse_slr.m` | `shaped_pulse_slr()` | Shinnar-Le Roux band-selective 90-degree pulse on a chain of 31 strongly coupled protons. Calculation time: seconds | 81 |
| `examples/shaped_pulses/shaped_pulse_vg.m` | `shaped_pulse_vg()` | Veshtort-Griffin E1000B 90-degree selective pulse applied to a system of 31 proton spins with nearest-neighbor J co- upl | 83 |
| `examples/singlet_states/carbon_singlet.m` | `carbon_singlet()` | Singlet relaxation rate for the two triple bond carbons in cis-dimethylbut-2-ynedioate. Magnetic parameters com- puted w | 53 |
| `examples/singlet_states/decoherence_benzoquinone.m` | `decoherence_benzoquinone()` | Long-lived spin states in the para-benzoquinone molecule (4 protons, 256-dimensional Liouville space). The relaxation su | 49 |
| `examples/singlet_states/decoherence_bicyclopropylidene.m` | `decoherence_bicyclopropylidene()` | Long-lived spin states in the bicyclopropylidene molecule (8 protons, 65536-dimensional Liouville space). The relaxation | 49 |
| `examples/singlet_states/decoherence_diacetylene.m` | `decoherence_diacetylene()` | Long-lived spin states in the diacetylene molecule (2 protons, 4 carbons, 4096-dimensional Liouville space). The relaxat | 67 |
| `examples/singlet_states/decoherence_naphthalenetetrone.m` | `decoherence_naphthalenetetrone()` | Long-lived spin states in the napthalenetetrone molecule. (4 protons, 256-dimensional Liouville space). The relaxation s | 49 |
| `examples/singlet_states/decoherence_urea.m` | `decoherence_urea()` | A demonstration that the nitrogen singlet state in urea is not long-lived. The relaxation superoperator accounts for eve | 52 |
| `examples/singlet_states/dipolar_singlet.m` | `dipolar_singlet()` | A demonstration that the two-spin singet state is immune to dipolar relaxation. Full Redfield superoperator for dipolar  | 48 |
| `examples/singlet_states/eigenstate_analysis.m` | `eigenstate_analysis()` | Stationary state analysis for the spin system of allyl pyruvate, finding out which component of the singlet state commut | 100 |
| `examples/singlet_states/m2s_example.m` | `m2s_example()` | An example of the M2S sequence for a two-spin system. Calculation time: seconds | 45 |
| `examples/singlet_states/s2m_example.m` | `s2m_example()` | An example of the S2M sequence for a two-spin system. Calculation time: seconds | 45 |
| `examples/singlet_states/singlet_imaging_1.m` | `singlet_imaging_1()` | Singlet imaging in a system with one-dimensional diffusion and flow. Calculation time: minutes | 160 |
| `examples/singlet_states/warren_singlet.m` | `warren_singlet()` | A demonstration that long-lived states exist that are immune not only to dipolar and CSA, but also to quadrupolar relaxa | 45 |
| `examples/spin_chemistry/cidnp_flash_acquire.m` | `cidnp_flash_acquire()` | A model of the CIDNP magnetisation pumping process described in IK's paper: The system is pumped for 0.5 seconds, and th | 100 |
| `examples/spin_chemistry/cidnp_geminate.m` | `cidnp_geminate()` | A basic example of the geminate CIDNP effect simulation. Calculation time: seconds | 59 |
| `examples/spin_chemistry/cidnp_nz_1.m` | `cidnp_nz_1()` | Field dependence of geminate CIDNP from a radical pair in a viscous solvent, computed with Redfield theory and with the  | 124 |
| `examples/spin_chemistry/cidnp_pumping_1.m` | `cidnp_pumping_1()` | A simulation of the matrix in Equation 2 of IK's paper on chemically amplified NOEs (https://doi.org/10.1016/j.jmr.2004. | 69 |
| `examples/spin_chemistry/cidnp_pumping_2.m` | `cidnp_pumping_2()` | A simulation of Figure 2A in IK's paper on chemically amplified NOEs (https://doi.org/10.1016/j.jmr.2004.01.011). Calcul | 84 |
| `examples/spin_chemistry/singlet_yield_1.m` | `singlet_yield_1()` | Liquid state magnetic field effect simulation on a radical pair with four nuclei using exponential recombination ki- net | 56 |
| `examples/spin_chemistry/singlet_yield_2.m` | `singlet_yield_2()` | Liquid state magnetic field effect simulation on a radical pair with six equivalent nuclei using exponential recombi- na | 55 |
| `examples/spin_chemistry/singlet_yield_3.m` | `singlet_yield_3()` | Figure 1 from the paper by Timmel, Till, Brocklehurst, McLauchlan and Hore: Calculation time: seconds | 49 |
| `examples/spin_chemistry/singlet_yield_4.m` | `singlet_yield_4()` | Figure 3 from the paper by Till, Timmel, Brocklehurst and Hore: Note: the original paper only uses electron Zeeman opera | 54 |
| `examples/spin_chemistry/singlet_yield_5.m` | `singlet_yield_5()` | Figure 5 from the paper by Timmel, Till, Brocklehurst, McLauchlan and Hore: Calculation time: seconds | 47 |
| `examples/spin_chemistry/singlet_yield_anisotropy_1.m` | `singlet_yield_anisotropy_1()` | Singlet yield anisotropy calculation for a radical pair using exponential recombination kinetics model. Calculation time | 57 |
| `examples/spin_chemistry/singlet_yield_anisotropy_2.m` | `singlet_yield_anisotropy_2()` | Singlet yield anisotropy calculation for a radical pair using exponential recombination kinetics model. Calculation time | 73 |
| `examples/spin_chemistry/singlet_yield_anisotropy_3.m` | `singlet_yield_anisotropy_3()` | Singlet yield anisotropy calculation for a model radical pair reaction, Haberkorn recombination model. Run time: minutes | 72 |
| `examples/spin_chemistry/singlet_yield_nz_1.m` | `singlet_yield_nz_1()` | Magnetic field effect on a triplet-born benzophenone ketyl / thiyl radical pair in a viscous ionic liquid, computed with | 130 |
| `examples/spin_chemistry/singlet_yield_nz_2.m` | `singlet_yield_nz_2()` | Field dependence of the decay rate of a micelle-confined triplet-born benzophenone ketyl / alkyl radical pair, computed  | 112 |
| `examples/visualisation/cst_peptide_bond.m` | `cst_peptide_bond()` | Example of shielding tensor visualisation for a peptide bond. Gaussian log is parsed. Note: antisymmetric components of  | 33 |
| `examples/visualisation/cst_strychnine.m` | `cst_strychnine()` | Example of carbon shielding tensor visualisation for strychnine molecule. Gaussian log is parsed. Note: antisymmetric co | 28 |
| `examples/visualisation/efg_silicate.m` | `efg_silicate()` | Example of electric field gradient tensor visualisation for an aluminosilicate solid. CASTEP log is parsed. | 25 |
| `examples/visualisation/hfc_porphyrine.m` | `hfc_porphyrine()` | Example of proton hyperfine tensor visualisation for copper porphyrine. ORCA log is parsed. | 26 |
| `examples/visualisation/hfc_pyrene.m` | `hfc_pyrene()` | Example of carbon hyperfine tensor visualisation for pyrene cation radical. Gaussian log is parsed. | 26 |
