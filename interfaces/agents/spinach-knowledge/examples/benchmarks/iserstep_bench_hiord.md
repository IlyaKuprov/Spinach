# examples/benchmarks/iserstep_bench_hiord.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/benchmarks/iserstep_bench_hiord.m`
- Signature: `iserstep_bench_hiord()`
- Total lines: 191

## Purpose

Benchmarks iserstep higher-order methods on a chirped-frequency oscillator with radiation damping, that has a state-dependent, and a time-dependent evolution generator. Syntax: iserstep_bench_hiord()

## Physical / mathematical content

- Benchmark examples. These files stress-test Spinach performance, scaling, and numerical throughput on representative spin-dynamics workloads, so runtime, memory pressure, and solver/pathway choices are part of the intended content.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file also defines local helper function(s): `run_method()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: Set the chirp rate; implemented by `cr=2*pi*400`.
- Lines 21-22: Set the relaxation rates; implemented by `r1=10; r2=10`.
- Lines 24-25: Set the radiation damping rate; implemented by `rrd=40`.
- Lines 27-28: Bootstrap the object; implemented by `spin_system=bootstrap()`.
- Lines 30-32: Make Bloch-Maxwell generator (Liouvillian, including -1i factors); implemented by `G=@(t,mu)(-1i*[r2 cr*t 0; -cr*t r2 0; 0 0 r1] -1i*rrd*[mu(3) 0 0; 0 mu(3) 0; mu(1) mu(2) 0])`.
- Lines 34-35: Set the initial magnetisation; implemented by `mu0=euler2dcm(0,pi*178/180,0)*[0 0 1]'`.
- Lines 37-38: Set the reference grid size; implemented by `np_ref=2^12`.
- Lines 40-41: Set the reference time step; implemented by `dt_ref=0.5/(np_ref-1)`.
- Lines 43-44: Preallocate the reference trajectory array; implemented by `mu_traj=zeros(3,np_ref)`.
- Lines 46-47: Set the initial state; implemented by `mu_traj(:,1)=mu0`.
- Lines 49-50: Compute the reference trajectory using RKMK-DP8; implemented by `for n=2:np_ref`.
- Lines 55-56: Extract the reference final state; implemented by `mu_ref=mu_traj(:,end)`.
- Lines 58-59: Set the benchmark grid sizes; implemented by `np=ceil(2.^linspace(8,10.5,10))`.
- Lines 61-62: Preallocate the accuracy array; implemented by `bench=zeros(numel(np),7)`.
- Lines 64-65: Benchmarking loop; implemented by `for k=1:numel(np)`.
- Lines 67-68: Set the time step; implemented by `dt=0.5/(np(k)-1)`.
- Lines 70-71: Run PWCL propagation and time it; implemented by `mu=run_method(spin_system,G,mu0,dt,np(k),'PWCL',false)`.
- Lines 75-76: Run LG2 propagation and time it; implemented by `mu=run_method(spin_system,G,mu0,dt,np(k),'LG2',false)`.

### Control flow inferred from the code

- Line 50: `for` loop over `n=2:np_ref`.
- Line 65: `for` loop over `k=1:numel(np)`.
- Line 157: `for` loop over `n=1:size(bench,2)`.

### Key state/data transformations

- Lines 19: computes `cr` using `cr=2*pi*400`.
- Lines 22: computes `r1` using `r1=10; r2=10`.
- Lines 25: computes `rrd` using `rrd=40`.
- Lines 28: computes `spin_system` using `spin_system=bootstrap()`.
- Lines 31-32: computes `G` using `G=@(t,mu)(-1i*[r2 cr*t 0; -cr*t r2 0; 0 0 r1] -1i*rrd*[mu(3) 0 0; 0 mu(3) 0; mu(1) mu(2) 0])`.
- Lines 35: computes `mu0` using `mu0=euler2dcm(0,pi*178/180,0)*[0 0 1]'`.
- Lines 38: computes `np_ref` using `np_ref=2^12`.
- Lines 41: computes `dt_ref` using `dt_ref=0.5/(np_ref-1)`.
- Lines 44: computes `mu_traj` using `mu_traj=zeros(3,np_ref)`.
- Lines 47: computes `mu_traj(:,1)` using `mu_traj(:,1)=mu0`.
- Lines 51: computes `t_curr` using `t_curr=(n-2)*dt_ref`.
- Lines 52: computes `mu_traj(:,n)` using `mu_traj(:,n)=step(spin_system,{G,t_curr,'RKMK-DP8'},mu_traj(:,n-1),dt_ref)`.
- Lines 56: computes `mu_ref` using `mu_ref=mu_traj(:,end)`.
- Lines 59: computes `np` using `np=ceil(2.^linspace(8,10.5,10))`.
- Lines 62: computes `bench` using `bench=zeros(numel(np),7)`.
- Lines 68: computes `dt` using `dt=0.5/(np(k)-1)`.
- Lines 71: computes `mu` using `mu=run_method(spin_system,G,mu0,dt,np(k),'PWCL',false)`.
- Lines 72: computes `bench(k,1)` using `bench(k,1)=norm(mu-mu_ref)/norm(mu_ref)`.

### Local helper functions

- Line 179: `run_method()` — `function mu=run_method(spin_system,G,mu0,dt,np,method,use_step)`.
  - Representative operation: `mu=mu0`.
  - Representative operation: `for n=2:np`.

## Outputs

- (none) -produces a set of diagnostic plots, and prints the empirical
- convergence orders for each method

## Implementation structure

- Benchmarks iserstep higher-order methods on a chirped-frequency oscillator
- with radiation damping, that has a state-dependent, and a time-dependent
- evolution generator. Syntax:
- iserstep_bench_hiord()
- (none) -produces a set of diagnostic plots, and prints the empirical
- convergence orders for each method
- Set the chirp rate
- Set the relaxation rates
- Set the radiation damping rate
- Bootstrap the object
- Make Bloch-Maxwell generator (Liouvillian, including -1i factors)
- Set the initial magnetisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `bootstrap()`, `euler2dcm()`, `mu_traj()`, `step()`, `run_method()`, `bench()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `klegend()`, `ylim()`, `set()`, `orders()`, `polyfit()`.
