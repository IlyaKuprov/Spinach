# examples/parahydrogen/case_studies/hyperpolarised_deuterium/bubble_pulse_acquire.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/parahydrogen/case_studies/hyperpolarised_deuterium/bubble_pulse_acquire.m`
- Signature: `bubble_pulse_acquire()`
- Total lines: 120

## Purpose

Simulated PNL (partially negative line) spectrum of ortho-deuterium in the presence of a parahydrogena- tion catalyst. Bubbling is followed by a 45-degree pulse. Paper link to follow in due course.

## Physical / mathematical content

- Parahydrogen examples. The physical motif is highly non-Boltzmann singlet order imported from para-H2 and converted into observable nuclear magnetisation through hydrogenation, exchange, or catalytic transfer processes.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: Spin system; implemented by `sys.isotopes={'2H','2H','2H','2H'}`.
- Lines 18-19: Experimental chemical shifts; implemented by `inter.zeeman.scalar={4.55 4.55 -13.5 -16.5}`.
- Lines 21-22: Experimental J-couplings; implemented by `inter.coupling.scalar=cell(4,4)`.
- Lines 26-27: NQI tensors from a DFT calculation; implemented by `inter.coupling.matrix{3,3}=1e3*[ 108.2 0.1 28.1`.
- Lines 34-35: Cartesian coordinates, DFT calculation; implemented by `inter.coordinates={[]; []`.
- Lines 39-40: Kinetics; implemented by `inter.chem.parts={[1 2],[3 4]}`.
- Lines 45-46: Magnet field; implemented by `sys.magnet=7.05`.
- Lines 48-49: Simulation formalsim; implemented by `bas.formalism='sphten-liouv'`.
- Lines 52-53: Relaxation theory; implemented by `inter.relaxation={'redfield','t1_t2'}`.
- Lines 60-61: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 64-65: Relaxation analysis; implemented by `relaxan(spin_system)`.
- Lines 67-68: Spin states; implemented by `[S,~,Q]=deut_pair(spin_system,1,2)`.
- Lines 70-71: Initial spin state: nothing; implemented by `rho0=unit_state(spin_system)`.
- Lines 73-74: Assumptions; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 76-77: Evolution generators; implemented by `H=hamiltonian(spin_system)`.
- Lines 80-81: Free kinetics; implemented by `KF=kinetics(spin_system)`.
- Lines 83-84: Kinetics with bubbling (a guess, needs proper rate); implemented by `pumped_state=S+Q{1}+Q{2}+Q{3}+Q{4}+Q{5}; pumped_state(1)=0`.
- Lines 87-88: Run the bubbling for 7 seconds; implemented by `rho=evolution(spin_system,H+1i*R+1i*KB,[],rho0,7,1,'final')`.

### Key state/data transformations

- Lines 16: computes `sys.isotopes` using `sys.isotopes={'2H','2H','2H','2H'}`.
- Lines 19: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={4.55 4.55 -13.5 -16.5}`.
- Lines 22: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(4,4)`.
- Lines 23: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2} = 12.0`.
- Lines 24: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4} = 0.24`.
- Lines 27: computes `inter.coupling.matrix{3,3}` using `inter.coupling.matrix{3,3}=1e3*[ 108.2 0.1 28.1`.
- Lines 30: computes `inter.coupling.matrix{4,4}` using `inter.coupling.matrix{4,4}=1e3*[ -55.2 8.1 -14.5`.
- Lines 35: computes `inter.coordinates` using `inter.coordinates={[]; []`.
- Lines 40: computes `inter.chem.parts` using `inter.chem.parts={[1 2],[3 4]}`.
- Lines 41-42: computes `inter.chem.rates` using `inter.chem.rates=[-1 5000; 1 -5000]`.
- Lines 43: computes `inter.chem.concs` using `inter.chem.concs=[1 0]`.
- Lines 46: computes `sys.magnet` using `sys.magnet=7.05`.
- Lines 49: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 50: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 53: computes `inter.relaxation` using `inter.relaxation={'redfield','t1_t2'}`.
- Lines 54: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 55: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 56: computes `inter.tau_c` using `inter.tau_c={1e-12, 400e-12}`.

## Implementation structure

- Simulated PNL (partially negative line) spectrum of
- ortho-deuterium in the presence of a parahydrogena-
- tion catalyst. Bubbling is followed by a 45-degree
- pulse. Paper link to follow in due course.
- Spin system
- Experimental chemical shifts
- Experimental J-couplings
- NQI tensors from a DFT calculation
- Cartesian coordinates, DFT calculation
- Kinetics
- Magnet field
- Simulation formalsim

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxan()`, `deut_pair()`, `unit_state()`, `assume()`, `hamiltonian()`, `relaxation()`, `kinetics()`, `pumped_state()`, `magpump()`, `evolution()`, `state()`, `operator()`, `liquid()`, `apodisation()`.
