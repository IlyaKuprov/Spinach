# examples/nmr_overtone/cpmas_glycine_simple.m

- MATLAB implementation: [examples/nmr_overtone/cpmas_glycine_simple.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_overtone/cpmas_glycine_simple.m)

- Signature: `cpmas_glycine_simple()`

## What it models

This is a single-condition Spinach calculation of proton-to-`14N` overtone cross-polarisation in glycine under MAS, according to the source header. It constructs a spin system from parameters in the script and calls `singlerot(spin_system,@overtone_cp,parameters,'qnmr')`; it does not load or fit a measured spectrum. The sequence implementation is delegated to `@overtone_cp`, which is not defined in this wrapper, so pulse shapes, phase cycling, and other sequence internals are not specified here. No gradient list or laboratory pulse-acquire schedule is defined in this wrapper.

The header says the glycine quadrupolar tensor data come from O'Dell and Ratcliffe, [DOI 10.1016/j.cplett.2011.08.030](https://doi.org/10.1016/j.cplett.2011.08.030). The input is `C_q=1.18e6` Hz (1.18 MHz), `eta_q=0.53`, and `I=1`, converted with `eeqq2nqi`. This paper is the named source for a model tensor; the example contains no experimental data file.

## Spin system and relaxation

The isotopes are `14N` and `1H`; the field input is `sys.magnet=14.1`. The wrapper sets `inter.zeeman.scalar={32.4,0}`, places the nuclei at `[0,0,0]` and `[0,0,1.00]`, and leaves the homonuclear `inter.coupling.matrix{2,2}` empty. No additional N–H tensor is assigned by hand: `create` uses these 1 Å coordinates to generate the point-dipolar N–H coupling automatically. Adding a second dipolar tensor would double-count it. Relaxation uses the damping option, diagonal retained terms, zero equilibrium, and `damp_rate=300` (no unit is given inline).

The basis is spherical-tensor Liouville space with no approximation. The wrapper disables `krylov` and `trajlevel`, sets `max_rank=7`, and uses the rough powder grid `rep_2ang_6400pts_sph`. Its rotor-rate input is `-19840`; the wrapper does not attach a unit to that literal.

## Fixed RF condition and simulated spectrum

The spectral sweep is 44-52 kHz with 256 points and 256-point zero filling. The axis is identified as kHz. The initial state is an oriented `1H` state, while the receiver and overtone channel operators are built for `14N`. The RF frequency input is `48e3` (48 kHz); the contact-duration input is `1e-4` s. The two-channel RF power input is `2*pi*[55.0e3,35.1e3]/sin(theta)`, where `theta=atan(sqrt(2))` is the magic angle.

The program runs one `singlerot` calculation, plots the real part of the simulated spectrum, and does not report any measured signal or fitted parameter. “Hours” in the source header is a runtime estimate only.
