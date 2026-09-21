# examples/esr_sol_pulsed/mas_diamond_p1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/mas_diamond_p1.m`
- Signature: `mas_diamond_p1()`
- Total lines: 95

## Purpose

Two-pulse echo-detected frequency-swept EPR spectra of the P1 substitutional nitrogen defect in diamond, static and under magic angle spinning, after Figure 1a of Khamrui et al., J. Phys. Chem. Lett. 2026, <https://doi.org/10.1021/acs.jpclett.6c02108>: 400 ns long pulses with 416 kHz nutation frequency, 300 ns interpulse delay, at 6.9 T, static and at 10, 25, and 37 kHz MAS. The carrier is stepped across the spectrum and the integrated echo is recorded at each frequency, with the spectra normalised to the static one. Calculation time: hours on a 256-core node.

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The 14N hyperfine coupling of the P1 centre (dipolar part 10.9 MHz) makes the two outer lines dephase under spinning as their resonance frequencies move during the sequence, whereas the central line survives.
- All P1 centre tensors are axial and coaxial, so a two-angle powder grid is sufficient; relaxation is omitted because it scales the four spectra by the same factor, which the normalisation removes.

## Known limitations against Figure 1a

- Line positions and the collapse of the outer lines under spinning match the paper, but the intensities do not reproduce the paper's own simulation: normalised to the static central peak, the static outer perpendicular edges come out at 0.20 and 0.18 (paper about 0.4), and the central line keeps 0.98, 0.90, and 0.82 of its echo at 10, 25, and 37 kHz (paper about 0.87, 0.52, and 0.33). The paper's simulated outer lines are broad humps where this example gives the perpendicular-edge singularity convolved with the roughly 1 MHz pulse response.
- Zeroing the 14N quadrupole leaves the 37 kHz central survival unchanged, so it is not the source. Untested candidates: the paper's model keeps only the secular hyperfine term with no nuclear Zeeman interaction, its echo integration window is not stated (this example integrates the complex echo over 1.0 us after the second pulse), and its carrier step is not stated. Relaxation (T1=100 us, T2=4 us in the paper) is omitted here because it scales all four spectra by the same factor.

## Numerical / algorithmic content

- Time propagation is explicit in Hilbert space (`zeeman-hilb`): the pulse sequence `experiments/echo_sweep.m` steps through the Hamiltonian rotor stack that `singlerot()` supplies, applying a propagator per 5 ns time step to the density matrix.
- The rotor stack index at each time step follows from `parameters.rate`, so `rate=0` is the static case; the rotor rank (2700) is set so that the rotor phase resolution of the stack matches the 5 ns time step at 37 kHz, the fastest rate used; at slower rates consecutive steps reuse a stack element.
- The electron coherence pathway (-1 after the first pulse, +1 after the second) is selected with `coherence()` in place of the phase cycle; the sequence averages over 100 rotor phases at the start of the sequence and sweeps the carrier offset inside the sequence, reusing the free-evolution propagators because the rotor stack commutes with the electron `Lz`.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: P1 centre parameters; implemented by `p1_params.orientation='111'`.
- Lines 36-37: Build the spin system; implemented by `[sys,inter]=diamond_p1(p1_params)`.
- Lines 39-40: Magnet field, central line at 193.797 GHz; implemented by `sys.magnet=6.9156`.
- Lines 42-43: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 46-47: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 50-51: Rotor parameters; implemented by `parameters.axis=[1 1 1]`.
- Lines 54-55: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 72-73: Spinning rates; implemented by `rates=[0 10e3 25e3 37e3]`.
- Lines 75-76: Simulation; implemented by `spectra=zeros(parameters.npoints,numel(rates))`.
- Lines 82-83: Normalisation to the static spectrum; implemented by `spectra=spectra/max(spectra(:,1))`.
- Lines 85-86: Plotting; implemented by `kfigure(); hold on`.

### Control flow inferred from the code

- Line 77: `for` loop over `n=1:numel(rates)`.
- Line 87: `for` loop over `n=1:numel(rates)`.

### Key state/data transformations

- Lines 33: computes `p1_params.orientation` using `p1_params.orientation='111'`.
- Lines 34: computes `p1_params.nitrogen` using `p1_params.nitrogen='14N'`.
- Lines 37: computes `[sys,inter]` using `[sys,inter]=diamond_p1(p1_params)`.
- Lines 40: computes `sys.magnet` using `sys.magnet=6.9156`.
- Lines 43: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 44: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 47: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 51: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 52: computes `parameters.max_rank` using `parameters.max_rank=2700`.
- Lines 55: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 56: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','E')`.
- Lines 57: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 58: computes `parameters.pulse_dur` using `parameters.pulse_dur=400e-9`.
- Lines 59: computes `parameters.pulse_frq` using `parameters.pulse_frq=416e3`.
- Lines 60: computes `parameters.tau` using `parameters.tau=300e-9`.
- Lines 61: computes `parameters.echo_win` using `parameters.echo_win=1.0e-6`.
- Lines 62: computes `parameters.timestep` using `parameters.timestep=5e-9`.
- Lines 63: computes `parameters.nphases` using `parameters.nphases=100`.
- Lines 64-70: computes `parameters.offset`, `parameters.sweep`, `parameters.npoints`, `parameters.zerofill`, `parameters.grid`, `parameters.axis_units`, and `parameters.verbose` using `parameters.offset=0`, `parameters.sweep=3e8`, `parameters.npoints=601`, `parameters.zerofill=601`, `parameters.grid='rep_2ang_400pts_sph'`, `parameters.axis_units='GHz-labframe'`, and `parameters.verbose=0`.
- Lines 78-79: computes `spectra(:,n)` using `spectra(:,n)=abs(singlerot(spin_system,@echo_sweep,parameters,'esr'))` after `parameters.rate=rates(n)`.
- Lines 83: computes `spectra` using `spectra=spectra/max(spectra(:,1))`.

## Implementation structure

- Two-pulse echo-detected frequency-swept EPR spectra of the P1 sub-
- stitutional nitrogen defect in diamond, static and under magic ang-
- le spinning, after Figure 1a of Khamrui et al., J. Phys. Chem. Lett.
- 2026, <https://doi.org/10.1021/acs.jpclett.6c02108>: 400 ns long
- pulses with 416 kHz nutation frequency, 300 ns interpulse delay, at
- 6.9 T, static and at 10, 25, and 37 kHz MAS. The carrier is stepped
- across the spectrum and the integrated echo is recorded at each fre-
- quency, with the spectra normalised to the static one.
- The 14N hyperfine coupling (dipolar part 10.9 MHz) makes the two outer
- lines dephase under spinning as their resonance frequencies move during
- the sequence, whereas the central line survives. Compared to the paper's
- own simulation, the static outer edges here are weaker (0.2 of the cent-
- ral peak against 0.4) and the central line dephases less (0.8 at 37 kHz).
- Calculation time: hours on a 256-core node.
- P1 centre parameters
- Build the spin system
- Magnet field, central line at 193.797 GHz
- Basis set
- Spinach housekeeping
- Rotor parameters
- Sequence parameters
- Spinning rates
- Simulation
- Normalisation to the static spectrum
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `diamond_p1()`, `create()`, `basis()`, `state()`, `singlerot()` with the `echo_sweep()` handle, `kfigure()`, `plot_1d()`, `klegend()`, `kylabel()`.
