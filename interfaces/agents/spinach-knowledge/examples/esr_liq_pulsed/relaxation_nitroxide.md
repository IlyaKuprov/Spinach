# examples/esr_liq_pulsed/relaxation_nitroxide.m

**Call:** <code>relaxation_nitroxide()</code> — no input arguments and no returned outputs.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_liq_pulsed/relaxation_nitroxide.m) · [DFT input](https://github.com/IlyaKuprov/Spinach/blob/main/examples/standard_systems/nitroxide.log)

## What it models

A W-band nitroxide pulse-acquire ESR spectrum. The source imports a DFT-derived system from <code>../standard_systems/nitroxide.log</code> with <code>gparse</code> and <code>g2spinach</code>, maps an electron (<code>E</code>) and <code>14N</code>, and sets <code>options.no_xyz=1</code> because hyperfine couplings are provided. It sets <code>sys.magnet=3.5</code>; the example labels this magnet induction but does not state its unit.

The basis is <code>sphten-liouv</code> with <code>approximation='none'</code>. Relaxation is <code>redfield</code>, with <code>equilibrium='zero'</code>, <code>rlx_keep='secular'</code>, and <code>tau_c={5e-11}</code>. These are the model choices encoded by the example, not a claim of validation against measurement.

## Acquisition and plotted result

The electron channel is observed: <code>spins={'E'}</code>, with <code>rho0</code> and <code>coil</code> both <code>state(...,'L+','E')</code> and no decoupled spins. Acquisition parameters are <code>offset=-2e8</code>, <code>sweep=2e8</code>, <code>npoints=512</code>, and <code>zerofill=1024</code>; the source sets <code>axis_units='GHz-labframe'</code>, <code>derivative=1</code>, and <code>invert_axis=1</code>. The numeric offset, sweep, and magnet values have no unit annotation in this file.

The function calls <code>liquid(spin_system,@acquire,parameters,'esr')</code>, applies <code>{{'none'}}</code> apodisation, computes <code>fftshift(fft(fid,parameters.zerofill))</code>, and plots the real spectrum. It opens a figure; the FID and spectrum are local variables, not returned or saved by this function. The source comments estimate seconds of calculation time.

## Requirements and scope

Run with Spinach and the example's relative input layout available (the relative <code>../standard_systems/nitroxide.log</code> path must resolve). The source uses <code>gparse</code>, <code>g2spinach</code>, <code>create</code>, <code>basis</code>, <code>state</code>, <code>liquid</code>, <code>acquire</code>, <code>apodisation</code>, <code>kfigure</code>, and <code>plot_1d</code>, plus MATLAB FFT routines. The .m file supplies no DOI or bibliography entry; the stable source and input links above are the cited records.
