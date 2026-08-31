# examples/nmr_liquids/noe_two_spin_hom.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/noe_two_spin_hom.m`
- Signature: `noe_two_spin_hom()`
- Total lines: 56

## Purpose

Nuclear overhauser effect in a homonuclear two-spin system in the long correlation time case. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Set the spin system; implemented by `sys.isotopes={'1H','1H'}`.
- Lines 16-17: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 23-24: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 30-31: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 34-35: Build the relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 37-38: Get thermal equilibrium state; implemented by `rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'))`.
- Lines 40-41: Start in a state with one spin inverted; implemented by `Lz1=state(spin_system,{'Lz'},{1})`.
- Lines 44-46: Compute the evolution trajectory; implemented by `coil=[state(spin_system,{'Lz'},{1}) state(spin_system,{'Lz'},{2})]`.
- Lines 49-50: Plot the longitudinal magnetization; implemented by `kfigure(); plot(linspace(0,10,1001),answer)`.

### Key state/data transformations

- Lines 11: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 12: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0 0.0}`.
- Lines 13: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 17: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 24: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 25: computes `inter.equilibrium` using `inter.equilibrium='dibari'`.
- Lines 26: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 27: computes `inter.temperature` using `inter.temperature=298`.
- Lines 28: computes `inter.tau_c` using `inter.tau_c={1e-9}`.
- Lines 31: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 35: computes `R` using `R=relaxation(spin_system)`.
- Lines 38: computes `rho_eq` using `rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'))`.
- Lines 41: computes `Lz1` using `Lz1=state(spin_system,{'Lz'},{1})`.
- Lines 42: computes `rho` using `rho=rho_eq-2*Lz1*(Lz1'*rho_eq)/norm(Lz1)^2`.
- Lines 45-46: computes `coil` using `coil=[state(spin_system,{'Lz'},{1}) state(spin_system,{'Lz'},{2})]`.
- Lines 47: computes `answer` using `answer=evolution(spin_system,1i*R,coil,rho,1e-2,1000,'multichannel')`.

## Implementation structure

- Nuclear overhauser effect in a homonuclear two-spin system in
- the long correlation time case.
- Calculation time: seconds
- Set the spin system
- Magnet field
- Basis set
- Relaxation theory parameters
- Spinach housekeeping
- Build the relaxation superoperator
- Get thermal equilibrium state
- Start in a state with one spin inverted
- Compute the evolution trajectory

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `equilibrium()`, `hamiltonian()`, `assume()`, `state()`, `evolution()`, `kfigure()`, `kxlabel()`, `kylabel()`, `klegend()`.
