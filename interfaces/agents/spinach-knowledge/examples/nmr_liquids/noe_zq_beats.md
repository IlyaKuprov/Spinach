# examples/nmr_liquids/noe_zq_beats.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/noe_zq_beats.m`
- Signature: `noe_zq_beats()`
- Total lines: 60

## Purpose

Zero-quantum beats in the Overhauser effect in a strongly coupled two-spin system. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Set the spin system; implemented by `sys.isotopes={'1H','1H'}`.
- Lines 17-18: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 31-32: Proximity cut-off; implemented by `sys.tols.prox_cutoff=4.0`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Build the Liouvillian; implemented by `L=hamiltonian(assume(spin_system,'nmr'))+1i*relaxation(spin_system)`.
- Lines 41-42: Get thermal equilibrium state; implemented by `rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'))`.
- Lines 44-45: Start in a state with one spin inverted; implemented by `Lz1=state(spin_system,{'Lz'},{1})`.
- Lines 48-49: Compute the evolution trajectory; implemented by `coil=[state(spin_system,{'Lz'},{1}) state(spin_system,{'Lz'},{2})]`.
- Lines 52-53: Plot the longitudinal magnetization; implemented by `kfigure()`.

### Key state/data transformations

- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 13: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0 0.01}`.
- Lines 14: computes `inter.coupling.scalar` using `inter.coupling.scalar={0.0 3.0; 0.0 0.0}`.
- Lines 15: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 18: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 26: computes `inter.equilibrium` using `inter.equilibrium='dibari'`.
- Lines 27: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 28: computes `inter.temperature` using `inter.temperature=298`.
- Lines 29: computes `inter.tau_c` using `inter.tau_c={1e-9}`.
- Lines 32: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 39: computes `L` using `L=hamiltonian(assume(spin_system,'nmr'))+1i*relaxation(spin_system)`.
- Lines 42: computes `rho_eq` using `rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'))`.
- Lines 45: computes `Lz1` using `Lz1=state(spin_system,{'Lz'},{1})`.
- Lines 46: computes `rho` using `rho=rho_eq-2*Lz1*(Lz1'*rho_eq)/norm(Lz1)^2`.

## Implementation structure

- Zero-quantum beats in the Overhauser effect in a strongly
- coupled two-spin system.
- Calculation time: seconds
- Set the spin system
- Magnet field
- Basis set
- Relaxation theory parameters
- Proximity cut-off
- Spinach housekeeping
- Build the Liouvillian
- Get thermal equilibrium state
- Start in a state with one spin inverted

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `relaxation()`, `equilibrium()`, `state()`, `evolution()`, `kfigure()`, `kxlabel()`, `kylabel()`, `klegend()`.
