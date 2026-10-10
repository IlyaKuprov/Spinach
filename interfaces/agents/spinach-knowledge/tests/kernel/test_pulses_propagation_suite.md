# tests/kernel/test_pulses_propagation_suite.m

## Purpose

Regression test suite for pulse-coordinate and propagation helper functions in Spinach. It verifies RF coordinate Hessian round-trips, Iserles generators, Lie-step methods on a constant generator, and R-sequence phase/compiler invariants.

## Behaviour

- Announces the test target with `fprintf` and initialises a regression test result object via `new_test_result` for `kernel/pulses_propagation_suite`.
- **Polar-Cartesian round-trip**: Uses non-zero amplitudes `r=[1.2 2.3 3.4]` and phases `p=[0.2 -0.7 1.1]` (away from the polar singularity). Defines `f=sum(r.^2)=sum(x.^2+y.^2)` whose Cartesian Hessian is exactly `2I`, with gradient `Dr=2*r` and all other derivative blocks zero. Converts polar coordinates, gradients, and Hessians to Cartesian via `polar2cartesian` and back via `cartesian2polar`, then checks:
  - Amplitude and phase Hessian round-trips (tolerance `1e-12`).
  - Cartesian gradients `Dx=2*x`, `Dy=2*y` (tolerance `1e-12`).
  - Hessian blocks `Dxx=2I`, `Dyy=2I`, `Dxy=0`, `Dyx=0` (tolerance `1e-12`).
  - Gradient round-trips for amplitude and phase (tolerance `1e-10`).
  - Hessian round-trips `Drr`, `Drp`, `Dpr`, `Dpp` (tolerance `1e-9`).
- **Iserles generators**: With Pauli-like matrices `HL=[0 1;1 0]`, `HM=[0 -1i;1i 0]`, `HR=[1 0;0 -1]` and `dt=0.125`:
  - Second order: `isergen(HL,[],HR,dt)` must equal `(HL+HR)/2+(1i*dt/6)*(HL*HR-HR*HL)` (tolerance `1e-15`).
  - Fourth order: `isergen(HL,HM,HR,dt)` must equal `(HL+4*HM+HR)/6+(1i*dt/12)*(HL*HR-HR*HL)` (tolerance `1e-15`).
- **Lie-step methods**: Builds a one-proton Hilbert-space spin system (`sys.magnet=0`, `sys.isotopes={'1H'}`, `inter.zeeman.scalar={0}`, `bas.formalism='zeeman-hilb'`, `bas.approximation='none'`) via `test_spin_system`. Uses `L=operator(spin_system,'Lx',1)`, `rho=state(spin_system,'Lz',1)`, `dt=0.25`, and the exact step `step(spin_system,L,rho,dt)`. For each method in `{'PWCL','LG2','LG4','RKMK4','LG4A'}`, `iserstep(spin_system,{Lfun,0,method},rho,dt)` with a constant generator `Lfun=@(~,~)L` must match the exact exponential step (tolerance `1e-12`).
- **R-sequence phases**: Calls `rsequence(1,4,1,1,1000,'180_pulse','homo_double_quantum_nucycle')`. Checks:
  - Phases equal `[base;-base]` with `base=[pi/4;-pi/4;pi/4;-pi/4]` (tolerance `1e-15`), i.e., R4_1 phase alternation is `+/-pi*nu/N` with the nucycle appending the inverted block.
  - Pulse amplitude equals `4000*pi` (tolerance `1e-12`), i.e., a 180-degree pulse over one quarter rotor period has nutation `pi/duration`.
  - Pulse duration equals `1/(4*1000)` (tolerance `1e-15`), i.e., the R element duration is n rotor periods divided by N blocks.
- **R-sequence compiler**: With Pauli matrices `S=pauli(2)`, `Sx=S.x`, `Sy=S.y`, calls `rseq_compiler(spin_system,zeros(2),Sx,Sy,[0;pi;0],0,0.1,'180_pulse')`. Checks:
  - Index map `T` equals `[1;2;1]` (tolerance `1e-15`), i.e., unique phases are compiled once and reused by their index map.
  - Each compiled propagator `P{n}` equals `speye(2)` (tolerance `1e-15`), i.e., with zero RF amplitude and zero drift each compiled propagator is the identity.

## Inputs and outputs

**Syntax**

```matlab
result = test_pulses_propagation_suite()
```

**Outputs**

- `result` — regression test result object with explanatory messages, accumulated through repeated `test_close` calls.

## References

- Source: [tests/kernel/test_pulses_propagation_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_pulses_propagation_suite.m)
