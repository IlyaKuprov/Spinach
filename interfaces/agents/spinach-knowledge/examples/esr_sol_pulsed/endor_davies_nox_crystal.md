# examples/esr_sol_pulsed/endor_davies_nox_crystal.m

**Call:** <code>endor_davies_nox_crystal()</code> — no input arguments and no returned outputs.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/endor_davies_nox_crystal.m)

## Spin system and model

A single-orientation nitroxide Davies ENDOR calculation, preceded by a pulse-acquire crystal ESR spectrum. The source describes soft-pulse treatment with the Fokker–Planck formalism. It defines an <code>E</code>/<code>14N</code> system at <code>sys.magnet=3.35</code>, with electron Zeeman matrix <code>[2.01045 0 0; 0 2.00641 0; 0 0 2.00211]</code> and electron–nitrogen coupling matrix <code>1e7*[1.2356 0 0.6322; 0 1.1266 0; 0.6322 0 8.2230]</code>.

The file does not annotate a unit for the magnet or matrix entries. The basis is <code>sphten-liouv</code> without approximation. Relaxation is <code>t1_t2</code>, with <code>rlx_keep='diagonal'</code>, <code>equilibrium='zero'</code>, and both <code>r1_rates</code> and <code>r2_rates</code> set to <code>{20e6 0.5e6}</code>.

## Stage 1: single-crystal ESR

The source sets <code>spins={'E'}</code>, electron <code>L+</code> initial state and coil, no decoupling, <code>offset=0</code>, <code>sweep=5e8</code>, <code>npoints=128</code>, <code>zerofill=512</code>, <code>axis_units='MHz'</code>, <code>orientation=[0 0 0]</code>, <code>derivative=0</code>, and <code>invert_axis=1</code>. It runs <code>crystal(spin_system,@acquire,parameters,'esr')</code>, applies no apodisation, FFTs the FID, and plots the real spectrum as the top-left panel. The source names the plot “Single crystal ESR spectrum.”

## Stage 2: Davies ENDOR

After clearing the first-stage parameters, the function uses <code>spins={'E','14N'}</code>, <code>rho0=state(...,'Lz','E')</code>, an electron <code>L+</code> coil, orientation <code>[0 0 0]</code>, offsets <code>[0 0]</code>, and <code>method='expm'</code>. Electron-pulse settings are <code>e_rnk=2</code>, <code>e_phi=0</code>, <code>e_dur=10e-9</code>, <code>e_pwr=2*pi*16.5e7</code>; nuclear-pulse settings are <code>n_rnk=3</code>, <code>n_phi=0</code>, <code>n_frq=linspace(-200e6,200e6,200)</code>, <code>n_dur=1e-7</code>, and <code>n_pwr=pi*1e7</code>. The source does not label units for these parameter literals. It evaluates <code>crystal(...,@endor_davies,parameters,'esr')</code> at <code>e_frq=+93e6</code>, <code>+10e6</code>, and <code>-73e6</code>.

The output is a four-panel figure: ESR plus three ENDOR traces. ENDOR curves plot <code>real(answer)</code> against <code>n_frq/1e6</code>, with the horizontal axis labelled nuclear frequency in MHz and vertical axis <code>(RF on)/(RF off)</code>. The function returns no arrays and writes no data file. The source estimates minutes of calculation time.

## Requirements and scope

Requires MATLAB and Spinach's <code>create</code>, <code>basis</code>, <code>state</code>, <code>crystal</code>, <code>acquire</code>, <code>endor_davies</code>, <code>apodisation</code>, <code>kfigure</code>, <code>scale_figure</code>, and plotting helpers. The result is for the one orientation explicitly set to <code>[0 0 0]</code>, not a powder average. The .m file supplies no DOI or bibliography entry; the stable source link above is provided instead.
