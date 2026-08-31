# examples/fundamentals/quadratures/product_quadratures_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/quadratures/product_quadratures_2.m`
- Signature: `product_quadratures_2()`
- Total lines: 140

## Purpose

A test of Lie-group product quadratures on a chirped frequency oscillator with radiation damping that has a state-dependent and time-dependent evolution generator.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Set system parameters; implemented by `cr=2*pi*400`.
- Lines 15-16: Bootstrap the object; implemented by `spin_system=bootstrap()`.
- Lines 18-20: Make Bloch-Maxwell generator; implemented by `G=@(t,mu)(-1i*[r2 cr*t 0; -cr*t r2 0; 0 0 r1] -1i*rrd*[mu(3) 0 0; 0 mu(3) 0; mu(1) mu(2) 0])`.
- Lines 22-23: Set initial magnetisation; implemented by `mu0=euler2dcm(0,pi*178/180,0)*[0 0 1]'`.
- Lines 25-26: Run reference RKMK-DP8 simulation and keep the trajectory; implemented by `np_ref=2^12; dt_ref=0.5/(np_ref-1); mu_traj=zeros(3,np_ref)`.
- Lines 34-35: Benchmark arrays; implemented by `np=ceil(2.^linspace(8,11,10)); bench=zeros(numel(np),7)`.
- Lines 37-38: Benchmarking loop; implemented by `for k=1:numel(np)`.
- Lines 40-41: Half a second; implemented by `dt=0.5/(np(k)-1)`.
- Lines 43-44: Piecewise-constant, left edge; implemented by `mu=mu0`.
- Lines 51-52: Lie group, order 2; implemented by `mu=mu0`.
- Lines 59-60: Lie group, order 4; implemented by `mu=mu0`.
- Lines 67-68: Lie group, order 4 (from Casas-Iserles Appendix A1); implemented by `mu=mu0`.
- Lines 75-76: Fourth-order RKMK; implemented by `mu=mu0`.
- Lines 83-84: Fifth-order Dormand-Prince RKMK; implemented by `mu=mu0`.
- Lines 91-92: Eighth-order Dormand-Prince RKMK; implemented by `mu=mu0`.
- Lines 101-102: Plotting; implemented by `kfigure(); scale_figure([1.8 1.0])`.
- Lines 123-124: Compute empirical convergence orders; implemented by `idx=floor(2*numel(np)/3):numel(np)`.
- Lines 128-129: Print empirical convergence orders; implemented by `fprintf('\nEmpirical convergence orders (slope of log-log error curve):\n')`.

### Control flow inferred from the code

- Line 28: `for` loop over `n=2:np_ref`.
- Line 38: `for` loop over `k=1:numel(np)`.
- Line 45: `for` loop over `n=2:np(k)`.
- Line 53: `for` loop over `n=2:np(k)`.
- Line 61: `for` loop over `n=2:np(k)`.
- Line 69: `for` loop over `n=2:np(k)`.
- Line 77: `for` loop over `n=2:np(k)`.
- Line 85: `for` loop over `n=2:np(k)`.
- Line 93: `for` loop over `n=2:np(k)`.
- Line 131: `for` loop over `n=1:size(bench,2)`.

### Key state/data transformations

- Lines 11: computes `cr` using `cr=2*pi*400`.
- Lines 12: computes `r1` using `r1=10; r2=10`.
- Lines 13: computes `rrd` using `rrd=40`.
- Lines 16: computes `spin_system` using `spin_system=bootstrap()`.
- Lines 19-20: computes `G` using `G=@(t,mu)(-1i*[r2 cr*t 0; -cr*t r2 0; 0 0 r1] -1i*rrd*[mu(3) 0 0; 0 mu(3) 0; mu(1) mu(2) 0])`.
- Lines 23: computes `mu0` using `mu0=euler2dcm(0,pi*178/180,0)*[0 0 1]'`.
- Lines 26: computes `np_ref` using `np_ref=2^12; dt_ref=0.5/(np_ref-1); mu_traj=zeros(3,np_ref)`.
- Lines 27: computes `mu_traj(:,1)` using `mu_traj(:,1)=mu0`.
- Lines 29: computes `t_curr` using `t_curr=(n-2)*dt_ref`.
- Lines 30: computes `mu_traj(:,n)` using `mu_traj(:,n)=iserstep(spin_system,{G,t_curr,'RKMK-DP8'},mu_traj(:,n-1),dt_ref)`.
- Lines 32: computes `mu_ref` using `mu_ref=mu_traj(:,end)`.
- Lines 35: computes `np` using `np=ceil(2.^linspace(8,11,10)); bench=zeros(numel(np),7)`.
- Lines 41: computes `dt` using `dt=0.5/(np(k)-1)`.
- Lines 44: computes `mu` using `mu=mu0`.
- Lines 49: computes `bench(k,1)` using `bench(k,1)=norm(mu-mu_ref)/norm(mu_ref)`.
- Lines 57: computes `bench(k,2)` using `bench(k,2)=norm(mu-mu_ref)/norm(mu_ref)`.
- Lines 65: computes `bench(k,3)` using `bench(k,3)=norm(mu-mu_ref)/norm(mu_ref)`.
- Lines 73: computes `bench(k,4)` using `bench(k,4)=norm(mu-mu_ref)/norm(mu_ref)`.

## Implementation structure

- A test of Lie-group product quadratures on a chirped frequency
- oscillator with radiation damping that has a state-dependent
- and time-dependent evolution generator.
- Set system parameters
- Bootstrap the object
- Make Bloch-Maxwell generator
- Set initial magnetisation
- Run reference RKMK-DP8 simulation and keep the trajectory
- Benchmark arrays
- Benchmarking loop
- Half a second
- Piecewise-constant, left edge

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `bootstrap()`, `euler2dcm()`, `mu_traj()`, `iserstep()`, `bench()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `klegend()`, `ylim()`, `set()`, `orders()`, `polyfit()`.
