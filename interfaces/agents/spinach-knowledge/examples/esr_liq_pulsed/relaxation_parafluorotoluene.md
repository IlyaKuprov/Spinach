# examples/esr_liq_pulsed/relaxation_parafluorotoluene.m

**Call:** <code>relaxation_parafluorotoluene()</code> — no input arguments and no returned outputs.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_liq_pulsed/relaxation_parafluorotoluene.m) · [DFT input](https://github.com/IlyaKuprov/Spinach/blob/main/examples/standard_systems/parafluorotoluene.log)

## What it models

An X-band pulse-acquire FFT ESR spectrum of the para-fluorotoluene radical, using explicit time-domain propagation and a Redfield relaxation superoperator. The spin system is imported from <code>../standard_systems/parafluorotoluene.log</code> with <code>gparse</code> and <code>g2spinach</code>; the isotope mapping is <code>E</code>, <code>1H</code>, and <code>19F</code>, and <code>options.no_xyz=1</code> because hyperfine couplings are supplied. The source sets <code>sys.magnet=0.33</code>; its comment calls this a magnet field but does not state a unit.

The basis is <code>sphten-liouv</code> with <code>approximation='none'</code>. Relaxation is Redfield with <code>equilibrium='zero'</code>, <code>rlx_keep='secular'</code>, and <code>tau_c={1e-10}</code>. The .m file does not list the imported tensors explicitly; they come from the linked DFT log through the parser.

## Acquisition and plotted result

The electron channel is observed, with <code>rho0</code> and <code>coil</code> both <code>L+</code> electron states and <code>decouple={}</code>. The script sets <code>offset=0</code>, <code>sweep=3e8</code>, <code>npoints=1024</code>, <code>zerofill=4096</code>, <code>axis_units='GHz-labframe'</code>, <code>derivative=1</code>, and <code>invert_axis=1</code>. The numeric magnet and sweep values are shown as coded; the file does not annotate their units.

The function calls <code>liquid(spin_system,@acquire,parameters,'esr')</code>, applies <code>{{'none'}}</code> apodisation, computes a shifted FFT, and plots its real part. It returns no spectrum or FID and does not save the plotted figure. The source estimates hours of calculation and requires at least 16 GB of RAM per CPU core; plan resources accordingly.

## Requirements and scope

Run with Spinach, MATLAB, and the relative <code>../standard_systems/parafluorotoluene.log</code> input path resolving from the example directory. The source uses <code>gparse</code>, <code>g2spinach</code>, <code>create</code>, <code>basis</code>, <code>state</code>, <code>liquid</code>, <code>acquire</code>, <code>apodisation</code>, <code>kfigure</code>, and <code>plot_1d</code>, as well as MATLAB FFT routines. The example supplies no DOI or bibliography record; the source and DFT-input links above are stable repository paths.
