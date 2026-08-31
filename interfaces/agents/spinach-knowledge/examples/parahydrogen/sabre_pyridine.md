# examples/parahydrogen/sabre_pyridine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/parahydrogen/sabre_pyridine.m`
- Signature: `sabre_pyridine()`
- Total lines: 114

## Purpose

SABRE experiment simulation for Eibe Duecker and Christian Griesinger. Set to reproduce Figure 3b from http://dx.doi.org/10.1021/ja903601p Calculation time: minutes

## Physical / mathematical content

- Parahydrogen examples. The physical motif is highly non-Boltzmann singlet order imported from para-H2 and converted into observable nuclear magnetisation through hydrogenation, exchange, or catalytic transfer processes.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Spin system; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 13-14: Chemical shifts; implemented by `inter.zeeman.scalar={8.54 7.44 7.86 7.44 8.54 -23.5 -23.5}`.
- Lines 16-17: Couplings inside pyridine; implemented by `inter.coupling.scalar=cell(7)`.
- Lines 29-30: Couplings of the hydride group; implemented by `inter.coupling.scalar{1,6}=1.12`.
- Lines 34-35: Magnetic fields; implemented by `polarization_field=25e-3`.
- Lines 38-39: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 42-43: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 45-46: Do the housekeeping; implemented by `sys.magnet=polarization_field`.
- Lines 50-51: Get the Hamiltonian superoperator; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 53-54: Start in a singlet state; implemented by `rho=singlet(spin_system,6,7)`.
- Lines 56-57: Evolve the system for 2.5 seconds; implemented by `rho=evolution(spin_system,H,[],rho,2.5,1,'final')`.
- Lines 59-60: Disconnect the parahydrogen; implemented by `[H,rho]=decouple(spin_system,H,rho,[6 7])`.
- Lines 62-63: Evolve the system for a further 2.5 seconds; implemented by `rho=evolution(spin_system,H,[],rho,2.5,1,'final')`.
- Lines 65-66: Split the Hamiltonian superoperator; implemented by `Hz=hamiltonian(assume(spin_system,'nmr','zeeman'))`.
- Lines 71-72: Lift the field exponentially in 1024 steps over 5 seconds; implemented by `for field=exp(linspace(log(polarization_field),log(nmr_field),1024))`.
- Lines 76-77: Set the field to NMR field; implemented by `H=Hc+nmr_field*Hz`.
- Lines 79-80: Evolve the system for a further second in high field; implemented by `rho=evolution(spin_system,H,[],rho,1.0,1,'final')`.
- Lines 82-83: Set NMR experiment parameters; implemented by `parameters.offset=2400`.

### Control flow inferred from the code

- Line 72: `for` loop over `field=exp(linspace(log(polarization_field),log(nmr_field),1024))`.

### Key state/data transformations

- Lines 11: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 14: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={8.54 7.44 7.86 7.44 8.54 -23.5 -23.5}`.
- Lines 17: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(7)`.
- Lines 18: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}= 4.88`.
- Lines 19: computes `inter.coupling.scalar{4,5}` using `inter.coupling.scalar{4,5}= 4.88`.
- Lines 20: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}= 1.00`.
- Lines 21: computes `inter.coupling.scalar{2,5}` using `inter.coupling.scalar{2,5}= 1.00`.
- Lines 22: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}= 1.84`.
- Lines 23: computes `inter.coupling.scalar{3,5}` using `inter.coupling.scalar{3,5}= 1.84`.
- Lines 24: computes `inter.coupling.scalar{1,5}` using `inter.coupling.scalar{1,5}=-0.13`.
- Lines 25: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}= 7.67`.
- Lines 26: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}= 7.67`.
- Lines 27: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}= 1.37`.
- Lines 30: computes `inter.coupling.scalar{1,6}` using `inter.coupling.scalar{1,6}=1.12`.
- Lines 31: computes `inter.coupling.scalar{1,7}` using `inter.coupling.scalar{1,7}=1.02`.
- Lines 32: computes `inter.coupling.scalar{6,7}` using `inter.coupling.scalar{6,7}=7.00`.
- Lines 35: computes `polarization_field` using `polarization_field=25e-3`.
- Lines 36: computes `nmr_field` using `nmr_field=7.05`.

## Implementation structure

- SABRE experiment simulation for Eibe Duecker and Christian Griesinger.
- Set to reproduce Figure 3b from http://dx.doi.org/10.1021/ja903601p
- Calculation time: minutes
- Spin system
- Chemical shifts
- Couplings inside pyridine
- Couplings of the hydride group
- Magnetic fields
- Basis set
- Algorithmic options
- Do the housekeeping
- Get the Hamiltonian superoperator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `singlet()`, `evolution()`, `decouple()`, `step()`, `state()`, `operator()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
