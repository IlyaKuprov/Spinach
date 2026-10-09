# examples/esr_liq_pulsed/relaxation_parafluoronitrobenzene.m

**Call:** <code>relaxation_parafluoronitrobenzene()</code> — no input arguments and no returned outputs.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_liq_pulsed/relaxation_parafluoronitrobenzene.m)

## What it models

A pulse-acquire FFT ESR simulation of para-fluoronitrobenzene, described in the source as an EasySpin test-file version; the source acknowledges Stefan Stoll. The spin system is specified directly in this file: <code>E</code>, <code>14N</code>, <code>19F</code>, and four <code>1H</code> spins. It sets <code>sys.magnet=0.33898</code> (the source labels it a magnet field but gives no unit). The electron Zeeman principal values are <code>[2.0032 2.0012 2.0097]</code> with Euler angles <code>[0 0 0]</code>.

The electron couplings are encoded as follows (the factors of <code>1e6</code> are part of the source expressions): <code>14N</code>, <code>(40.40+[24 -12 -12])*1e6</code>; <code>19F</code>, <code>(22.51+[34.9 -19.8 -15])*1e6</code>; protons 4 and 5, <code>[9.69 9.69 9.69]*1e6</code>; protons 6 and 7, <code>[3.16 3.16 3.16]*1e6</code>. All listed coupling Euler angles are <code>[0 0 0]</code>. The basis is <code>sphten-liouv</code> with no approximation, zero-quantum terms for <code>1H</code>, and <code>S2</code> symmetry on proton pairs <code>[4 5]</code> and <code>[6 7]</code>.

Relaxation is secular Redfield (<code>inter.relaxation={'redfield'}</code>, <code>rlx_keep='secular'</code>) with <code>equilibrium='zero'</code> and <code>tau_c={160e-12}</code>. These are the model choices in the script, not an assertion of agreement with experiment.

## Acquisition and plotted result

The electron channel is observed with <code>rho0</code> and <code>coil</code> both <code>L+</code> electron states and <code>decouple={}</code>. The source sets <code>offset=-1e7</code>, <code>sweep=2e8</code>, <code>npoints=1024</code>, <code>zerofill=4096</code>, <code>axis_units='GHz-labframe'</code>, <code>derivative=1</code>, and <code>invert_axis=1</code>. These numeric field/frequency parameters are transcribed as coded; the file does not attach units to the magnet, offset, or sweep values.

It calls <code>liquid(spin_system,@acquire,parameters,'esr')</code>, applies <code>{{'none'}}</code> apodisation, takes the shifted FFT, and plots the real spectrum. The function returns no data and saves no figure or spectrum; the source estimates calculation time in minutes.

## Requirements and scope

Requires MATLAB and Spinach routines used by the script (<code>create</code>, <code>basis</code>, <code>state</code>, <code>liquid</code>, <code>acquire</code>, <code>apodisation</code>, <code>kfigure</code>, <code>plot_1d</code>) and MATLAB FFT routines. Unlike the DFT-import examples, its spin tensors are entered directly, so it has no auxiliary system-file dependency. The source describes it as an EasySpin test-file version but supplies no DOI or external citation URL; the stable source link above is provided instead.
