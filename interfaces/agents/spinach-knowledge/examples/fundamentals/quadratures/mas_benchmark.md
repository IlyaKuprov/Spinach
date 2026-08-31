# examples/fundamentals/quadratures/mas_benchmark.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/quadratures/mas_benchmark.m`
- Signature: `mas_benchmark()`
- Total lines: 129

## Purpose

Integrating the Lioville -von Neumann equation through one period of the MAS rotor using the piecewise-constant Hamiltonian approximation as well as the more accurate Lie group integrators. Calculation time: seconds

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: System specification; implemented by `sys.magnet=14.1`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 29-30: Get Hamiltonian; implemented by `[H,Q]=hamiltonian(spin_system)`.
- Lines 32-33: Spinning speed, Hz; implemented by `rate=50000; T=1/rate`.
- Lines 35-36: Magic angle; implemented by `mag=atan(sqrt(2))`.
- Lines 38-39: Initial condition; implemented by `rho_init=state(spin_system,'L+','1H')`.
- Lines 41-42: Discretise the period; implemented by `np_ref=2^13+1`.
- Lines 46-47: Pre-allocate H(t); implemented by `HF=cell(1,numel(t_ref))`.
- Lines 49-50: Pre-compute Euler angles and H(t); implemented by `for n=1:numel(t_ref)`.
- Lines 55-56: Run reference simulation; implemented by `rho_ref=rho_init`.
- Lines 63-64: Pre-allocate error array; implemented by `np=2.^(4:12)+1; error=zeros(numel(np),4)`.
- Lines 66-67: Benchmarking loop -maybe parfor; implemented by `for k=1:numel(np)`.
- Lines 69-70: Discretise the period; implemented by `t=linspace(0,T,np(k))`.
- Lines 73-74: Pre-allocate H(t); implemented by `HF=cell(1,numel(t))`.
- Lines 76-77: Pre-compute Euler angles and H(t); implemented by `for n=1:numel(t)`.
- Lines 81-82: Piecewise constant, left point; implemented by `rho_one=rho_init`.
- Lines 88-89: Piecewise constant, mid point; implemented by `rho_one=rho_init`.

### Control flow inferred from the code

- Line 50: `for` loop over `n=1:numel(t_ref)`.
- Line 57: `for` loop over `n=1:2:(numel(t_ref)-2)`.
- Line 67: `for` loop over `k=1:numel(np)`.
- Line 77: `for` loop over `n=1:numel(t)`.
- Line 83: `for` loop over `n=1:2:(numel(t)-1)`.
- Line 90: `for` loop over `n=1:2:(numel(t)-1)`.
- Line 97: `for` loop over `n=1:2:(numel(t)-2)`.
- Line 105: `for` loop over `n=1:2:(numel(t)-2)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 16: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={5.0,-2.0}`.
- Lines 17: computes `inter.coordinates` using `inter.coordinates={[0 0 0]; [0 3.9 0.1]}`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 22: computes `bas.projections` using `bas.projections=+1`.
- Lines 25: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 30: computes `[H,Q]` using `[H,Q]=hamiltonian(spin_system)`.
- Lines 33: computes `rate` using `rate=50000; T=1/rate`.
- Lines 36: computes `mag` using `mag=atan(sqrt(2))`.
- Lines 39: computes `rho_init` using `rho_init=state(spin_system,'L+','1H')`.
- Lines 42: computes `np_ref` using `np_ref=2^13+1`.
- Lines 43: computes `t_ref` using `t_ref=linspace(0,T,np_ref)`.
- Lines 44: computes `dt` using `dt=T/(numel(t_ref)-1)`.
- Lines 47: computes `HF` using `HF=cell(1,numel(t_ref))`.
- Lines 51: computes `HF{n}` using `HF{n}=H+orientation(Q,[0 mag 2*pi*(n-1)/(numel(t_ref)-1)])`.
- Lines 56: computes `rho_ref` using `rho_ref=rho_init`.

## Implementation structure

- Integrating the Lioville -von Neumann equation through
- one period of the MAS rotor using the piecewise-constant
- Hamiltonian approximation as well as the more accurate
- Lie group integrators.
- Calculation time: seconds
- System specification
- Basis set
- Spinach housekeeping
- Get Hamiltonian
- Spinning speed, Hz
- Magic angle
- Initial condition

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `hamiltonian()`, `atan()`, `state()`, `orientation()`, `step()`, `kfigure()`, `set()`, `kxlabel()`, `kylabel()`, `klegend()`, `ylim()`, `scale_figure()`.
