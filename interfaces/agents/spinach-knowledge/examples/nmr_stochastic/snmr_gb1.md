# examples/nmr_stochastic/snmr_gb1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_stochastic/snmr_gb1.m`
- Signature: `snmr_gb1()`
- Total lines: 165

## Purpose

A Primas-style stochastic NMR experiment on GB1 protein. The calculation requires a terabyte of RAM and NVidia A100 GPU. Calculation time: hours

## Physical / mathematical content

- Stochastic NMR examples. These scripts model random processes, trajectories, or stochastic Liouville dynamics and connect fluctuating Hamiltonians or transport processes to ensemble-averaged observables.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Protein data import; implemented by `options.pdb_mol=1`.
- Lines 16-17: Magnet field; implemented by `sys.magnet=18.79`.
- Lines 19-20: Tolerances; implemented by `sys.tols.inter_cutoff=2.0`.
- Lines 23-24: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 29-30: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 36-40: Use GPU arithmetic sys.enable={'gpu'};; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 43-44: Get the Hamiltonian; implemented by `H0=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 46-47: Get the relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 49-50: Save the workspace; implemented by `save('gb1_workspace.mat','-v7.3','-nocompression')`.
- Lines 52-53: Control operators; implemented by `Hx=operator(spin_system,'Lx','1H')`.
- Lines 60-63: Observables; implemented by `coils=[state(spin_system,'Lx','1H') state(spin_system,'Ly','1H') state(spin_system,'Lx','13C') state(spin_system,'Ly','13C') state(spin_system,'Lx','15N') state(spin_sys…`.
- Lines 65-66: Isotropic thermal equilibrium; implemented by `rho=equilibrium(spin_system)`.
- Lines 68-69: Re-save the workspace; implemented by `save('gb1_workspace.mat','-v7.3','-nocompression')`.
- Lines 71-72: Stochastic process parameters; implemented by `dt=1e-5`.
- Lines 76-77: Control noise track; implemented by `cHx=sig_omega*randn(nsteps,1)`.
- Lines 84-85: GPU uploads; implemented by `L0=gpuArray(H0+1i*R)`.
- Lines 91-92: Trajectory calculation; implemented by `report(spin_system,'computing trajectory ')`.

### Control flow inferred from the code

- Line 94: `for` loop over `n=1:nsteps`.
- Line 108: conditional branch on `mod(n,1000)==0`.

### Key state/data transformations

- Lines 11: computes `options.pdb_mol` using `options.pdb_mol=1`.
- Lines 12: computes `options.noshift` using `options.noshift='delete'`.
- Lines 13: computes `options.select` using `options.select='all'`.
- Lines 14: computes `[sys,inter]` using `[sys,inter]=protein('2N9K.pdb','2N9K.bmrb',options)`.
- Lines 17: computes `sys.magnet` using `sys.magnet=18.79`.
- Lines 20: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=2.0`.
- Lines 21: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 24: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 25: computes `bas.approximation` using `bas.approximation='IK-1'`.
- Lines 26: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 27: computes `bas.level` using `bas.level=4; bas.space_level=3`.
- Lines 30: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 31: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 32: computes `inter.equilibrium` using `inter.equilibrium='IME'`.
- Lines 33: computes `inter.tau_c` using `inter.tau_c={5e-9}`.
- Lines 34: computes `inter.temperature` using `inter.temperature=298`.
- Lines 40: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 44: computes `H0` using `H0=hamiltonian(assume(spin_system,'nmr'))`.

## Implementation structure

- A Primas-style stochastic NMR experiment on GB1 protein. The
- calculation requires a terabyte of RAM and NVidia A100 GPU.
- Calculation time: hours
- Protein data import
- Magnet field
- Tolerances
- Basis set
- Relaxation theory
- Use GPU arithmetic
- sys.enable={'gpu'};
- Spinach housekeeping
- Get the Hamiltonian

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `protein()`, `create()`, `basis()`, `hamiltonian()`, `assume()`, `relaxation()`, `save()`, `operator()`, `state()`, `equilibrium()`, `gpuArray()`, `report()`, `fids()`, `gather()`, `cHx()`, `cHy()`.
