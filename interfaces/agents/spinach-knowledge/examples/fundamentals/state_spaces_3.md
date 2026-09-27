# examples/fundamentals/state_spaces_3.m

- Signature: `state_spaces_3()`

## Purpose

Simulate transverse magnetisation in the fatty-acid system produced by `fatty_acid(15)` under scalar-coupling evolution and repeated 180° pulses, without relaxation, and inspect how the state spreads through correlation orders. The source estimates minutes of calculation time.

## Method

The example uses the `sphten-liouv` formalism, `IK-2` basis approximation, proximal level 1, scalar-coupling connectivity, and a 14.1 T field. It starts from proton Lx magnetisation and observes with proton L+. After 50 trajectory points at `4e-5` s intervals, it applies eight pi-rotation pulses about Lx; each pulse is followed by 100 further points at the same interval. The combined trajectory is analysed with `trajan(...,'correlation_order')`.

The system enables `greedy` and `prop_cache`. No relaxation operator or relaxation evolution is included.
