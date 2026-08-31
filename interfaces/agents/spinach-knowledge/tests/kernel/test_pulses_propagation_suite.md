# tests/kernel/test_pulses_propagation_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_pulses_propagation_suite.m`
- Signature: `result=test_pulses_propagation_suite()`
- Total lines: 128

## Purpose

Tests pulse-coordinate and propagation helpers. Syntax: result=test_pulses_propagation_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Pulse coordinate and propagation helpers\n')`.
- Lines 20-23: State the propagation-helper target of the test; implemented by `result=new_test_result('kernel/pulses_propagation_suite', 'Pulse coordinate and propagation helpers', 'pulse propagation helpers must preserve analytic coordinate and pr…`.
- Lines 25-26: Choose non-zero amplitudes away from the polar singularity; implemented by `r=[1.2 2.3 3.4]`.
- Lines 29-30: Use f=sum(r.^2)=sum(x.^2+y.^2), whose Cartesian Hessian is exactly 2I; implemented by `Dr=2*r`.
- Lines 37-38: Convert polar coordinates, gradients, and Hessians to Cartesian and back; implemented by `[x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy]=polar2cartesian(r,p,Dr,Dp,Drr,Drp,Dpr,Dpp)`.
- Lines 70-71: Check Iserles second-order and fourth-order product quadrature formulae; implemented by `HL=[0 1;1 0]`.
- Lines 82-83: Build a one-proton Hilbert-space spin system for Lie-step checks; implemented by `sys.magnet=0`.
- Lines 95-96: Constant generators should reduce all low-order Lie methods to the exact exponential step; implemented by `methods={'PWCL','LG2','LG4','RKMK4','LG4A'}`.
- Lines 103-104: Check R-sequence phase construction and nutation calibration; implemented by `[phases,pulse_amp,pulse_dur]=rsequence(1,4,1,1,1000,'180_pulse','homo_double_quantum_nucycle')`.
- Lines 114-115: Check R-sequence compiler indexing when the RF amplitude is zero; implemented by `S=pauli(2)`.

### Control flow inferred from the code

- Line 97: `for` loop over `n=1:numel(methods)`.
- Line 121: `for` loop over `n=1:numel(P)`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/pulses_propagation_suite', 'Pulse coordinate and propagation helpers', 'pulse propagation helpers must preserve analytic coordinate and pr…`.
- Lines 26: computes `r` using `r=[1.2 2.3 3.4]`.
- Lines 27: computes `p` using `p=[0.2 -0.7 1.1]`.
- Lines 30: computes `Dr` using `Dr=2*r`.
- Lines 31: computes `Dp` using `Dp=zeros(size(r))`.
- Lines 32: computes `Drr` using `Drr=2*eye(numel(r))`.
- Lines 33: computes `Drp` using `Drp=zeros(numel(r),numel(r))`.
- Lines 34: computes `Dpr` using `Dpr=zeros(numel(r),numel(r))`.
- Lines 35: computes `Dpp` using `Dpp=zeros(numel(r),numel(r))`.
- Lines 38: computes `[x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy]` using `[x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy]=polar2cartesian(r,p,Dr,Dp,Drr,Drp,Dpr,Dpp)`.
- Lines 71: computes `HL` using `HL=[0 1;1 0]`.
- Lines 72: computes `HM` using `HM=[0 -1i;1i 0]`.
- Lines 73: computes `HR` using `HR=[1 0;0 -1]`.
- Lines 74: computes `dt` using `dt=0.125`.
- Lines 75: computes `H_ref` using `H_ref=(HL+HR)/2+(1i*dt/6)*(HL*HR-HR*HL)`.
- Lines 83: computes `sys.magnet` using `sys.magnet=0`.
- Lines 84: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 85: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.

## Outputs

- result -regression test result with explanatory messages
- The test checks RF coordinate Hessian round-trips, Iserles generators,
- Lie-step methods on a constant generator, and R-sequence phase/compiler
- invariants.

## Implementation structure

- Tests pulse-coordinate and propagation helpers. Syntax:
- result=test_pulses_propagation_suite()
- result -regression test result with explanatory messages
- The test checks RF coordinate Hessian round-trips, Iserles generators,
- Lie-step methods on a constant generator, and R-sequence phase/compiler
- invariants.
- Announce the test target
- State the propagation-helper target of the test
- Choose non-zero amplitudes away from the polar singularity
- Use f=sum(r.^2)=sum(x.^2+y.^2), whose Cartesian Hessian is exactly 2I
- Convert polar coordinates, gradients, and Hessians to Cartesian and back
- Check Iserles second-order and fourth-order product quadrature formulae

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `polar2cartesian()`, `cartesian2polar()`, `test_close()`, `isergen()`, `test_spin_system()`, `operator()`, `state()`, `step()`, `iserstep()`, `rsequence()`, `pauli()`, `rseq_compiler()`, `int2str()`, `speye()`.
