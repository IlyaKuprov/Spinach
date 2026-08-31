# examples/fundamentals/nutation_dist_test.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/nutation_dist_test.m`
- Signature: `nutation_dist_test()`
- Total lines: 98

## Purpose

Recovery of an RF field distribution from a nutation curve measured with the same coil used for excitation and detection. A proton en- semble is driven on-resonance by a bimodal distribution of RF field amplitudes; following the reciprocity principle, the detected sig- nal of every ensemble member is weighted by its own RF field ampli- tude. The transverse magnetisation components are combined into a complex nutation

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: Isotopes; implemented by `sys.isotopes={'1H'}`.
- Lines 21-22: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 24-25: Chemical shift; implemented by `inter.zeeman.scalar={0.0}`.
- Lines 27-28: No console output; implemented by `sys.output='hush'`.
- Lines 30-31: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Initial and observable states; implemented by `rho=state(spin_system,'Lz','1H')`.
- Lines 43-44: RF field operator; implemented by `Lx=operator(spin_system,'Lx','1H')`.
- Lines 46-47: Drift Hamiltonian; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 49-50: Bimodal RF field distribution in rad/s; implemented by `b1_freq=2*pi*linspace(25e3,65e3,201).'`.
- Lines 55-56: Probability masses on the RF field grid; implemented by `b1_mass=b1_dist*(b1_freq(2)-b1_freq(1))`.
- Lines 58-59: Nutation curve time grid; implemented by `dt=2e-6; npts=256`.
- Lines 61-62: Nutation curve accumulator; implemented by `curve=zeros(npts,1)`.
- Lines 64-65: Loop over the RF field ensemble; implemented by `for n=1:numel(b1_freq)`.
- Lines 67-69: Magnetisation trajectories for the current RF field; implemented by `traj=evolution(spin_system,H+b1_freq(n)*Lx,[coil_x coil_y], rho,dt,npts-1,'multichannel')`.
- Lines 71-72: Reciprocity-weighted accumulation of the nutation curve; implemented by `curve=curve+b1_mass(n)*b1_freq(n)*(traj(1,:)+1i*traj(2,:)).'`.
- Lines 76-77: Receiver phase and gain; implemented by `curve=exp(1i*1.9)*curve/max(abs(curve))`.
- Lines 79-80: Complex measurement noise; implemented by `rng(1); curve=curve+2e-3*(randn(npts,1)+1i*randn(npts,1))`.

### Control flow inferred from the code

- Line 65: `for` loop over `n=1:numel(b1_freq)`.

### Key state/data transformations

- Lines 19: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 22: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 25: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0}`.
- Lines 28: computes `sys.output` using `sys.output='hush'`.
- Lines 31: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 39: computes `rho` using `rho=state(spin_system,'Lz','1H')`.
- Lines 40: computes `coil_x` using `coil_x=state(spin_system,'Lx','1H')`.
- Lines 41: computes `coil_y` using `coil_y=state(spin_system,'Ly','1H')`.
- Lines 44: computes `Lx` using `Lx=operator(spin_system,'Lx','1H')`.
- Lines 47: computes `H` using `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 50: computes `b1_freq` using `b1_freq=2*pi*linspace(25e3,65e3,201).'`.
- Lines 51-52: computes `b1_dist` using `b1_dist=0.8*exp(-(b1_freq-2*pi*50e3).^2/(2*(2*pi*3.0e3)^2))+ 0.2*exp(-(b1_freq-2*pi*38e3).^2/(2*(2*pi*2.5e3)^2))`.
- Lines 56: computes `b1_mass` using `b1_mass=b1_dist*(b1_freq(2)-b1_freq(1))`.
- Lines 59: computes `dt` using `dt=2e-6; npts=256`.
- Lines 62: computes `curve` using `curve=zeros(npts,1)`.
- Lines 68-69: computes `traj` using `traj=evolution(spin_system,H+b1_freq(n)*Lx,[coil_x coil_y], rho,dt,npts-1,'multichannel')`.

## Implementation structure

- Recovery of an RF field distribution from a nutation curve measured
- with the same coil used for excitation and detection. A proton en-
- semble is driven on-resonance by a bimodal distribution of RF field
- amplitudes; following the reciprocity principle, the detected sig-
- nal of every ensemble member is weighted by its own RF field ampli-
- tude. The transverse magnetisation components are combined into a
- complex nutation curve, a receiver phase is applied, noise is add-
- ed, and nutation_dist.m is called with a user-specified Tikhonov
- regularisation parameter to recover the true nutation frequency
- distribution with the reception weight divided out.
- Calculation time: seconds
- Isotopes

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `trapz()`, `b1_freq()`, `evolution()`, `b1_mass()`, `traj()`, `rng()`, `nutation_dist()`, `kfigure()`, `kxlabel()`, `kgrid()`.
