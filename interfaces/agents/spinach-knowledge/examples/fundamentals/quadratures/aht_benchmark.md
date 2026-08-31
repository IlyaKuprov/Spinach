# examples/fundamentals/quadratures/aht_benchmark.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/quadratures/aht_benchmark.m`
- Signature: `aht_benchmark()`
- Total lines: 135

## Purpose

Accuracy benchmark for the Hamiltonian period propagator caclulation using Lie group integrators. Bacause the mo- dulation is purely sinusoidal, 2nd and 4th order integra- tors show the same apparent accuracy. Calculation time: seconds

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: System specification; implemented by `sys.magnet=14.1; sys.isotopes={'14N'}`.
- Lines 18-19: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 22-23: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 28-29: Algorithmic options; implemented by `sys.disable={'krylov','trajlevel'}`.
- Lines 31-32: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 35-36: Magic angle; implemented by `theta=atan(sqrt(2))`.
- Lines 38-39: Spectrum setup; implemented by `parameters.max_rank=6`.
- Lines 49-50: Pull the Hamiltonian out of the wrapper; implemented by `stuff=singlerot(spin_system,@impound,parameters,'qnmr')`.
- Lines 53-54: Get the overtone frequency; implemented by `ovt_frq=-2*spin('14N')*spin_system.inter.magnet/(2*pi)`.
- Lines 56-57: Get the pulse frequency; implemented by `omega=2*pi*ovt_frq-2*pi*parameters.rf_frq`.
- Lines 59-60: Project pulse operators; implemented by `Lx=kron(speye(parameters.spc_dim),parameters.Lx)`.
- Lines 62-63: Get the pulse Hamiltonian; implemented by `pulseop=parameters.rf_pwr*Lx; Hp=pulseop/2; Hm=pulseop/2`.
- Lines 65-66: Reference calculation using 4th order Lie quadrature; implemented by `spin_system.tols.prop_chop=eps()`.
- Lines 75-76: Pre-allocate error array; implemented by `nslices=2:32; error=zeros(numel(nslices),4)`.
- Lines 78-79: Loop over point count; implemented by `for k=1:numel(nslices)`.
- Lines 81-82: Discretise the period of the rotating frame; implemented by `slice_dur=2*pi/(omega*nslices(k))`.
- Lines 84-85: Left edge quadrature; implemented by `PL=eye(size(H0))`.
- Lines 92-93: Midpoint quadrature; implemented by `PM=eye(size(H0))`.

### Control flow inferred from the code

- Line 68: `for` loop over `n=1:nslices`.
- Line 79: `for` loop over `k=1:numel(nslices)`.
- Line 86: `parfor` loop over `n=1:nslices(k)`.
- Line 94: `parfor` loop over `n=1:nslices(k)`.
- Line 102: `parfor` loop over `n=1:nslices(k)`.
- Line 111: `parfor` loop over `n=1:nslices(k)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=14.1; sys.isotopes={'14N'}`.
- Lines 15: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0])`.
- Lines 16: computes `inter.zeeman.scalar{1}` using `inter.zeeman.scalar{1}=32.4`.
- Lines 19: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 20: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 23: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 24: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 25: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 26: computes `inter.damp_rate` using `inter.damp_rate=300`.
- Lines 29: computes `sys.disable` using `sys.disable={'krylov','trajlevel'}`.
- Lines 32: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 36: computes `theta` using `theta=atan(sqrt(2))`.
- Lines 39: computes `parameters.max_rank` using `parameters.max_rank=6`.
- Lines 40: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 41: computes `parameters.rate` using `parameters.rate=-19840`.
- Lines 42: computes `parameters.grid` using `parameters.grid='single_crystal'`.
- Lines 43-44: computes `parameters.Lx` using `parameters.Lx=cos(theta)*operator(spin_system,'Lz','14N')+ sin(theta)*operator(spin_system,'Lx','14N')`.
- Lines 45: computes `parameters.spins` using `parameters.spins={'14N'}`.

## Implementation structure

- Accuracy benchmark for the Hamiltonian period propagator
- caclulation using Lie group integrators. Bacause the mo-
- dulation is purely sinusoidal, 2nd and 4th order integra-
- tors show the same apparent accuracy.
- Calculation time: seconds
- System specification
- Basis set
- Relaxation theory
- Algorithmic options
- Spinach housekeeping
- Magic angle
- Spectrum setup

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `atan()`, `operator()`, `singlerot()`, `spin()`, `speye()`, `eps()`, `propagator()`, `isergen()`, `nslices()`, `kfigure()`, `set()`, `kxlabel()`, `kylabel()`.
