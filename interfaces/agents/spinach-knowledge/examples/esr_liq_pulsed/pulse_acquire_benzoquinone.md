# examples/esr_liq_pulsed/pulse_acquire_benzoquinone.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/pulse_acquire_benzoquinone.m`
- Signature: `pulse_acquire_benzoquinone()`
- Total lines: 75

## Purpose

Pulse-acquire FFT ESR on 2-methoxy-1,4-benzoquinone radical in liquid state. Set to reproduce Figure 1 in Simple common linewidth is used as a relaxation model. Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: Magnet induction; implemented by `sys.magnet=0.33`.
- Lines 18-19: Isotope list; implemented by `sys.isotopes={'E','1H','1H','1H','1H','1H','1H'}`.
- Lines 21-22: Zeeman interactions and couplings; implemented by `inter.zeeman.scalar={2.004577 0 0 0 0 0 0}`.
- Lines 31-32: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 37-38: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 45-46: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 49-50: Experiment parameters; implemented by `parameters.spins={'E'}`.
- Lines 62-63: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'esr')`.
- Lines 65-66: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'none'}})`.
- Lines 68-69: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 71-72: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 19: computes `sys.isotopes` using `sys.isotopes={'E','1H','1H','1H','1H','1H','1H'}`.
- Lines 22: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.004577 0 0 0 0 0 0}`.
- Lines 23: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(7,7)`.
- Lines 24: computes `inter.coupling.scalar{2,1}` using `inter.coupling.scalar{2,1}=mt2hz(0.08)`.
- Lines 25: computes `inter.coupling.scalar{3,1}` using `inter.coupling.scalar{3,1}=mt2hz(0.08)`.
- Lines 26: computes `inter.coupling.scalar{4,1}` using `inter.coupling.scalar{4,1}=mt2hz(0.08)`.
- Lines 27: computes `inter.coupling.scalar{5,1}` using `inter.coupling.scalar{5,1}=mt2hz(-0.059)`.
- Lines 28: computes `inter.coupling.scalar{6,1}` using `inter.coupling.scalar{6,1}=mt2hz(-0.364)`.
- Lines 29: computes `inter.coupling.scalar{7,1}` using `inter.coupling.scalar{7,1}=mt2hz(-0.204)`.
- Lines 32: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 33: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 34: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 35: computes `inter.damp_rate` using `inter.damp_rate=1e6`.
- Lines 38: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 39: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 40: computes `bas.longitudinals` using `bas.longitudinals={'1H'}`.
- Lines 41: computes `bas.projections` using `bas.projections=+1`.

## Implementation structure

- Pulse-acquire FFT ESR on 2-methoxy-1,4-benzoquinone radical in
- liquid state. Set to reproduce Figure 1 in
- Simple common linewidth is used as a relaxation model.
- Calculation time: seconds
- Magnet induction
- Isotope list
- Zeeman interactions and couplings
- Relaxation theory
- Basis set
- Spinach housekeeping
- Experiment parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `mt2hz()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
