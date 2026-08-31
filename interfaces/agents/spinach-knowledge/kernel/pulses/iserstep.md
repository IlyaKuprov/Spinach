# kernel/pulses/iserstep.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/iserstep.m`
- Signature: `rho_b=iserstep(spin_system,LTM,rho_a,dt)`
- Total lines: 517

## Purpose

Lie-group and Runge-Kutta-Munthe-Kaas solvers for the Lie equa- tion. LG methods are implementations of Equation A.1, with mi- nor typos fixed, from The key difference from step() function is that the Liouvillian can depend on the density matrix. Syntax: rho_b=iserstep(spin_system,{L,t,method},rho_a,dt)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `RKMK()`, `dexpinv()`, `bernoulli_B()`, `rk_tableau()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 48-49: Check consistency; implemented by `grumble(LTM,rho_a,dt)`.
- Lines 51-52: Dereference input pointers; implemented by `L=LTM{1}; t=LTM{2}; method=LTM{3}`.
- Lines 54-55: Decide method; implemented by `switch method`.
- Lines 59-60: Left generator from left state; implemented by `LL=L(t,rho_a)`.
- Lines 62-63: Assume the generator stays constant; implemented by `rho_b=step(spin_system,LL,rho_a,dt)`.
- Lines 70-71: Estimate midpoint state and generator; implemented by `rho_mid=step(spin_system,LL,rho_a,dt/2)`.
- Lines 74-75: Assume the generator stays constant; implemented by `rho_b=step(spin_system,LM,rho_a,dt)`.
- Lines 86-87: Estimate endpoint using midpoint generator; implemented by `rho_b=step(spin_system,LM,rho_a,dt)`.
- Lines 89-90: Take the step using fourth-order Lie method; implemented by `rho_b=step(spin_system,{LL,LM,L(t+dt,rho_b)},rho_a,dt)`.
- Lines 94-95: Return the input state for zero time step; implemented by `if dt==0`.
- Lines 100-101: Match the paper notation; implemented by `a_fun=@(t,rho)(-1i*L(t,rho))`.
- Lines 103-104: Convert Lie algebra elements into generators for step(); implemented by `h=dt; t_n=t; prop_scale=1i/dt`.
- Lines 106-107: Stage 1; implemented by `k1=h*a_fun(t_n,rho_a)`.
- Lines 110-111: Stage 2; implemented by `u2=0.5*Q1`.
- Lines 116-117: Stage 3; implemented by `u3=0.5*Q1+0.25*Q2`.
- Lines 122-123: Stage 4; implemented by `u4=Q1+Q2`.
- Lines 128-129: Stage 5 (Note: [A, B] is the commutator A*B -B*A); implemented by `u5=0.5*Q1+0.25*Q2+(1/3)*Q3-(1/24)*Q4-(1/48)*(Q1*Q2-Q2*Q1)`.
- Lines 134-135: Stage 6; implemented by `u6=Q1+Q2+(2/3)*Q3+(1/6)*Q4-(1/6)*(Q1*Q2-Q2*Q1)`.

### Control flow inferred from the code

- Line 55: dispatches on `method`; cases `'PWCL'`, `'LG2'`, `'LG4'`, `'LG4A'`, `'RKMK4'`, `'RKMK-DP5'`, `'RKMK-DP8'`, `'RKMK-RKF45'`.
- Line 95: conditional branch on `dt==0`.

### Key state/data transformations

- Lines 52: computes `L` using `L=LTM{1}; t=LTM{2}; method=LTM{3}`.
- Lines 60: computes `LL` using `LL=L(t,rho_a)`.
- Lines 63: computes `rho_b` using `rho_b=step(spin_system,LL,rho_a,dt)`.
- Lines 71: computes `rho_mid` using `rho_mid=step(spin_system,LL,rho_a,dt/2)`.
- Lines 72: computes `LM` using `LM=L(t+0.5*dt,rho_mid)`.
- Lines 101: computes `a_fun` using `a_fun=@(t,rho)(-1i*L(t,rho))`.
- Lines 104: computes `h` using `h=dt; t_n=t; prop_scale=1i/dt`.
- Lines 107: computes `k1` using `k1=h*a_fun(t_n,rho_a)`.
- Lines 108: computes `Q1` using `Q1=k1`.
- Lines 111: computes `u2` using `u2=0.5*Q1`.
- Lines 112: computes `rho_u2` using `rho_u2=step(spin_system,prop_scale*u2,rho_a,dt)`.
- Lines 113: computes `k2` using `k2=h*a_fun(t_n+h/2,rho_u2)`.
- Lines 114: computes `Q2` using `Q2=k2-k1`.
- Lines 117: computes `u3` using `u3=0.5*Q1+0.25*Q2`.
- Lines 118: computes `rho_u3` using `rho_u3=step(spin_system,prop_scale*u3,rho_a,dt)`.
- Lines 119: computes `k3` using `k3=h*a_fun(t_n+h/2,rho_u3)`.
- Lines 120: computes `Q3` using `Q3=k3-k2`.
- Lines 123: computes `u4` using `u4=Q1+Q2`.

### Local helper functions

- Line 194: `RKMK()` — `function rho_b=RKMK(order_or_name,spin_system,L,t,rho_a,dt,varargin)`. Process built-in method names
  - Representative operation: `if ischar(order_or_name)||isstring(order_or_name)`.
  - Representative operation: `method_name=char(order_or_name)`.
- Line 299: `dexpinv()` — `function w=dexpinv(u,v,q_dexp,comm_fun,bern_coeffs)`. Initialise the series and first adjoint
  - Representative operation: `w=v`.
  - Representative operation: `ad_term=v`.
- Line 317: `bernoulli_B()` — `function B=bernoulli_B(n)`. Set available Bernoulli numbers
  - Representative operation: `B_all=[1,-1/2,1/6,0,-1/30,0,1/42,0,-1/30,0,5/66]`.
  - Representative operation: `if n>10`.
- Line 333: `rk_tableau()` — `function [A,b,c,q]=rk_tableau(name)`. Standardise method name
  - Representative operation: `name=upper(strtrim(name))`.
  - Representative operation: `switch name`.
- Line 484: `grumble()` — `function grumble(LTM,rho_a,dt)`.
  - Representative operation: `if ~iscell(LTM)`.
  - Representative operation: `error('second input must be a cell array.')`.

## Parameters / inputs

- spin_system -Spinach data structure from create.m
- and basis.m constructors
- L -a handle to a function L(t,rho) that must take
- time and state vector, and return the evolution
- generator (in rad/s) of the Lie equation:
- d_rho/d_t = -i*L(t,rho)*rho
- rho_a -state vector at the start of the evolution
- period
- t -time at the start of the evolution, seconds
- dt -evolution time step, seconds
- method -'PWCL', 'PWCM', 'RKMK4', 'RKMK-DP5',
- 'RKMK-DP8', or 'LG4'; the latter one
- has a good balance of efficiency and
- numerical accuracy

## Outputs

- rho_b -state vector at the end of the evolution
- time step

## Implementation structure

- Lie-group and Runge-Kutta-Munthe-Kaas solvers for the Lie equa-
- tion. LG methods are implementations of Equation A.1, with mi-
- nor typos fixed, from
- The key difference from step() function is that the Liouvillian
- can depend on the density matrix. Syntax:
- rho_b=iserstep(spin_system,{L,t,method},rho_a,dt)
- spin_system -Spinach data structure from create.m
- and basis.m constructors
- L -a handle to a function L(t,rho) that must take
- time and state vector, and return the evolution
- generator (in rad/s) of the Lie equation:
- d_rho/d_t = -i*L(t,rho)*rho

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `step()`, `a_fun()`, `comm()`, `RKMK()`, `ischar()`, `isstring()`, `char()`, `rk_tableau()`, `lower()`, `bernoulli_B()`, `isequal()`, `dexpinv()`, `comm_fun()`, `bern_coeffs()`, `factorial()`.
