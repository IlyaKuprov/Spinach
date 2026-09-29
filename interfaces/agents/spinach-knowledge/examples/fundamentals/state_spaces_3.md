# examples/fundamentals/state_spaces_3.m

- Signature: `state_spaces_3()`
- Source: [`examples/fundamentals/state_spaces_3.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/state_spaces_3.m)

## Model and question

This example follows transverse magnetisation in the fatty-acid system returned by `fatty_acid(15)`, under strong scalar coupling and repeated refocusing pulses. Its analysis asks how the trajectory occupies successive spin-correlation orders; it does not test angular or powder quadrature. The source states that relaxation is absent and estimates a runtime of minutes.

The field is 14.1 T. The basis uses `sphten-liouv`, `IK-2`, proximity level 1, and scalar-coupling connectivity; `greedy` and `prop_cache` are enabled. The initial state is proton `Lx`; the source also constructs a proton `L+` coil state, but does not pass that coil into the trajectory calculation. The Hamiltonian is built under the NMR assumption, without a relaxation term.

## Sequence and output

Trajectory-mode `evolution` first propagates 50 points at 4e-5 s per point. Eight loop iterations then apply a `pi` rotation about proton `Lx`, each followed by 100 more trajectory points at the same time step. Thus the source specifies 850 sampled evolution points and nominally 0.034 s of free evolution, apart from the instantaneous pulse steps. `trajan(...,'correlation_order')` displays the resulting correlation-order trajectory.

The source defines no numerical pass/fail threshold, convergence test, or reference trace; any visual axis scaling is not an acceptance criterion. The page therefore describes the sequence and diagnostic output without claiming that the trajectory was run or that a particular correlation order is sufficient.
