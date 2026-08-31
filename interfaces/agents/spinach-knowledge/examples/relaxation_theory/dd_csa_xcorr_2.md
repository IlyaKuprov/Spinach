# examples/relaxation_theory/dd_csa_xcorr_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/dd_csa_xcorr_2.m`
- Signature: `dd_csa_xcorr_2()`
- Total lines: 103

## Purpose

DD-CSA cross-correlation -a reproduction of Fig 5a from the paper by Grace and Kumar (http://dx.doi.org/10.1006/jmra.1995.1151). Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-12: Read the spin system parameters (vacuum DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/fdnb.log'), {{'H','1H'},{'F','19F'}},[32.0 270.0],[])`.
- Lines 14-15: Set up the calculation; implemented by `sys.magnet=9.4`.
- Lines 26-27: Proximity cut-off; implemented by `sys.tols.prox_cutoff=4.0`.
- Lines 29-30: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 33-34: Set simulation parameters; implemented by `parameters.offset=-521`.
- Lines 40-41: Set the assumptions to high-field NMR; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 43-44: Get the Hamiltonian superoperator; implemented by `L=hamiltonian(spin_system)`.
- Lines 46-47: Add Redfield superoperator,; implemented by `L=L+1i*relaxation(spin_system)`.
- Lines 49-50: Apply the offset; implemented by `L=frqoffset(spin_system,L,parameters)`.
- Lines 52-53: Get pulse operators on fluorine; implemented by `Lx=operator(spin_system,'Lx','19F')`.
- Lines 56-57: Set detection state to L+; implemented by `coil=state(spin_system,'L+','19F')`.
- Lines 59-60: Calculate the time step of the simulation; implemented by `timestep=1/parameters.sweep`.
- Lines 62-63: Compute the axis; implemented by `axis_hz=sweep2ticks(parameters.offset,parameters.sweep,parameters.zerofill)`.
- Lines 65-66: Set the initial state to thermal equilibrium; implemented by `rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'))`.
- Lines 68-69: Get the figure going; implemented by `kfigure()`.
- Lines 71-72: Loop over mixing times; implemented by `for t_mix=[0.1 1.4 1.6 1.8 2.0 2.2 2.4 10]`.
- Lines 74-75: Apply the inversion pulse; implemented by `rho=step(spin_system,Lx,rho_eq,pi)`.
- Lines 77-78: Run the mixing time; implemented by `rho=evolution(spin_system,L,coil,rho,t_mix,1,'final')`.

### Control flow inferred from the code

- Line 72: `for` loop over `t_mix=[0.1 1.4 1.6 1.8 2.0 2.2 2.4 10]`.

### Key state/data transformations

- Lines 11-12: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/fdnb.log'), {{'H','1H'},{'F','19F'}},[32.0 270.0],[])`.
- Lines 15: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 16: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=5`.
- Lines 17: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 18: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 19: computes `inter.equilibrium` using `inter.equilibrium='dibari'`.
- Lines 20: computes `inter.temperature` using `inter.temperature=298`.
- Lines 21: computes `inter.tau_c` using `inter.tau_c={9.6e-12}`.
- Lines 23: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 34: computes `parameters.offset` using `parameters.offset=-521`.
- Lines 35: computes `parameters.sweep` using `parameters.sweep=50`.
- Lines 36: computes `parameters.npoints` using `parameters.npoints=128`.
- Lines 37: computes `parameters.zerofill` using `parameters.zerofill=512`.
- Lines 38: computes `parameters.spins` using `parameters.spins={'19F'}`.
- Lines 44: computes `L` using `L=hamiltonian(spin_system)`.
- Lines 53: computes `Lx` using `Lx=operator(spin_system,'Lx','19F')`.

## Implementation structure

- DD-CSA cross-correlation -a reproduction of Fig 5a from the paper
- by Grace and Kumar (http://dx.doi.org/10.1006/jmra.1995.1151).
- Calculation time: seconds
- Read the spin system parameters (vacuum DFT calculation)
- Set up the calculation
- Proximity cut-off
- Run Spinach housekeeping
- Set simulation parameters
- Set the assumptions to high-field NMR
- Get the Hamiltonian superoperator
- Add Redfield superoperator,
- Apply the offset

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `relaxation()`, `frqoffset()`, `operator()`, `state()`, `sweep2ticks()`, `equilibrium()`, `kfigure()`, `step()`, `evolution()`, `fftshift()`.
