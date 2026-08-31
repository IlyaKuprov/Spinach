# kernel/includes/redfield_integral_serial.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/includes/redfield_integral_serial.m`
- Signature: `(script file)`
- Total lines: 116

## Purpose

Bloch-Wangsness-Redfield and Nakajima-Zwanzig integral evaluati- on, the serial path. This include is called from within the rel- axation.m theory blocks and follows the notation used in IK's paper: with the difference that the numerical quadrature method propo- sed there has been superceded by the much faster auxiliary mat- rix method described in: The calling theory block must set rlx_onshell (true selects the back

## Physical / mathematical content

- Include scripts and shared setup fragments. These files implement tightly scoped runtime setup, parallel profiling, resource guards, or shared kernels included by other Spinach routines.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 1-29: Bloch-Wangsness-Redfield and Nakajima-Zwanzig integral evaluati- on, the serial path. This include is called from within the rel- axation.m theory blocks and follows the notation used in IK's paper: with the difference that the numerical quadrature method propo- sed there has been superceded by the much faster auxiliary mat- rix method described in: The calling theory block must set rlx_onshell (true selects the back-rotated kernel that reduces to Redfield theory at zero shift, false the resolvent kernel of Nakajima-Zwanzig theory) and rlx_shift (the Laplace evaluation point, Hz); Redfield the- ory is the on-shell form at zero shift. This include is called when relaxation.m is not at the top of the parallelisation call stack. When it is, the asynchronous include is called instead.; implemented by `for n=1:numel(Q)`.
- Lines 28-29: Sum over spherical ranks; implemented by `for n=1:numel(Q)`.
- Lines 31-32: Sum over first projection index; implemented by `for k=1:(2*n+1)`.
- Lines 34-35: Sum over second projection index; implemented by `for m=1:(2*n+1)`.
- Lines 37-38: Only proceed if necessary; implemented by `if cheap_norm(Q{n}{k,m})>0`.
- Lines 40-41: Sum over third projection index; implemented by `for p=1:(2*n+1)`.
- Lines 43-44: Sum over fourth projection index; implemented by `for q=1:(2*n+1)`.
- Lines 46-47: Only proceed if necessary; implemented by `if cheap_norm(Q{n}{p,q})>0`.
- Lines 49-50: Get decay weights and rates for correlation function; implemented by `[weights,rates,states]=corrfun(spin_system,n,k,m,p,q)`.
- Lines 52-53: Loop over chemical species; implemented by `for s=1:numel(states)`.
- Lines 55-56: Loop over correlation function exponentials; implemented by `for j=1:numel(weights{s})`.
- Lines 58-59: Only proceed if necessary; implemented by `if abs(weights{s}(j))>0`.
- Lines 61-62: Set the upper integration limit according to the accuracy goal; implemented by `upper_limit=-1.5*(1/rates{s}(j))*log(1/spin_system.tols.rlx_integration)`.
- Lines 64-67: Inform the user; implemented by `report(spin_system,['integrating L=' num2str(n) ', k=' num2str(n+1-k,'%+d') ', m=' num2str(n+1-m,'%+d') ', p=' num2str(n+1-p,'%+d') ', q=' num2str(n+1-q,'%+d') ', chemic…`.
- Lines 69-70: Kill the terms in L0 that are irrelevant on the time scale of the integration; implemented by `B=clean_up(spin_system,L0,spin_system.tols.rlx_integration/abs(upper_limit))`.
- Lines 72-73: Prepare the coupling matrices; implemented by `A=Q{n}{k,m}; C=Q{n}{p,q}'`.
- Lines 75-76: Kernel form and evaluation point set by the calling theory; implemented by `if rlx_onshell`.
- Lines 82-83: Obliterate irrelevant elements; implemented by `A(~states{s},~states{s})=0; B(~states{s},~states{s})=0`.

### Control flow inferred from the code

- Line 29: `for` loop over `n=1:numel(Q)`.
- Line 32: `for` loop over `k=1:(2*n+1)`.
- Line 35: `for` loop over `m=1:(2*n+1)`.
- Line 38: conditional branch on `cheap_norm(Q{n}{k,m})>0`.
- Line 41: `for` loop over `p=1:(2*n+1)`.
- Line 44: `for` loop over `q=1:(2*n+1)`.
- Line 47: conditional branch on `cheap_norm(Q{n}{p,q})>0`.
- Line 53: `for` loop over `s=1:numel(states)`.
- Line 56: `for` loop over `j=1:numel(weights{s})`.
- Line 59: conditional branch on `abs(weights{s}(j))>0`.
- Line 76: conditional branch on `rlx_onshell`.

### Key state/data transformations

- Lines 50: computes `[weights,rates,states]` using `[weights,rates,states]=corrfun(spin_system,n,k,m,p,q)`.
- Lines 62: computes `upper_limit` using `upper_limit=-1.5*(1/rates{s}(j))*log(1/spin_system.tols.rlx_integration)`.
- Lines 65-67: computes `report(spin_system,['integrating L` using `report(spin_system,['integrating L=' num2str(n) ', k=' num2str(n+1-k,'%+d') ', m=' num2str(n+1-m,'%+d') ', p=' num2str(n+1-p,'%+d') ', q=' num2str(n+1-q,'%+d') ', chemic…`.
- Lines 66-67: computes `', m` using `', m=' num2str(n+1-m,'%+d') ', p=' num2str(n+1-p,'%+d') ', q=' num2str(n+1-q,'%+d') ', chemical species ' num2str(s) ' '])`.
- Lines 67: computes `', q` using `', q=' num2str(n+1-q,'%+d') ', chemical species ' num2str(s) ' '])`.
- Lines 70: computes `B` using `B=clean_up(spin_system,L0,spin_system.tols.rlx_integration/abs(upper_limit))`.
- Lines 73: computes `A` using `A=Q{n}{k,m}; C=Q{n}{p,q}'`.
- Lines 77: computes `D` using `D=B-1i*(rates{s}(j)-rlx_shift)*speye(size(B))`.
- Lines 83: computes `A(~states{s},~states{s})` using `A(~states{s},~states{s})=0; B(~states{s},~states{s})=0`.
- Lines 84: computes `C(~states{s},~states{s})` using `C(~states{s},~states{s})=0; D(~states{s},~states{s})=0`.
- Lines 87: computes `R_int` using `R_int=-weights{s}(j)*A*expmint(spin_system,B,C,D,upper_limit)`.
- Lines 91: computes `R` using `R=R+R_int`.

## Implementation structure

- Bloch-Wangsness-Redfield and Nakajima-Zwanzig integral evaluati-
- on, the serial path. This include is called from within the rel-
- axation.m theory blocks and follows the notation used in IK's
- paper:
- with the difference that the numerical quadrature method propo-
- sed there has been superceded by the much faster auxiliary mat-
- rix method described in:
- The calling theory block must set rlx_onshell (true selects the
- back-rotated kernel that reduces to Redfield theory at zero
- shift, false the resolvent kernel of Nakajima-Zwanzig theory)
- and rlx_shift (the Laplace evaluation point, Hz); Redfield the-
- ory is the on-shell form at zero shift.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cheap_norm()`, `corrfun()`, `report()`, `num2str()`, `clean_up()`, `speye()`, `expmint()`, `clear()`.
