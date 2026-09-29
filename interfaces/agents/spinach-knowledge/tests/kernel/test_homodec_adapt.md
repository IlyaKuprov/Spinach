# tests/kernel/test_homodec_adapt.m

[View source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_homodec_adapt.m)

## Purpose

Regression test for the admission of acquisition-stage irradiation specified in Hilbert space into Liouville-space simulation. The test verifies that a complex-phase soft pulse with transverse acquisition irradiation produces the same FID after Hilbert admission as under native Liouville execution, and that zero-power, absent-irradiation, and no-op formalism cases are handled correctly. Optional operator caching is not enabled.

## Behaviour

- Builds a coupled `1H`–`13C` spin system with unequal Zeeman offsets (`inter.zeeman.scalar={1 3}`) and scalar coupling (`inter.coupling.scalar={0 150;0 0}`), at 9.4 T, with `bas.approximation='none'`.
- Loops over the formalisms `zeeman-hilb`, `zeeman-liouv`, `sphten-liouv`, and `zeeman-wavef`, building the Hamiltonian and setting `parameters.homodec_oper` to `cos(0.4)*Lx(1H) + sin(0.4)*Ly(1H)` with `parameters.homodec_pwr=100`.
- Checks that the drift Hamiltonian and the irradiation operator genuinely do not commute (Frobenius norm of the commutator greater than 1).
- Calls `sim2liouv` on each formalism. For `zeeman-hilb`, verifies that the admitted `homodec_oper` equals `hilb2liouv(parameters.homodec_oper,'comm')` to absolute and relative tolerances of `1e-14`. For the other formalisms, verifies a no-op: the converted system and admitted parameters are returned unchanged (`isequaln`).
- For the first two formalisms only, runs `liquid` with `@sp_acquire` using a soft pulse (`pulse_frq=50`, `pulse_phi=0.3`, `pulse_pwr=2*pi*1000`, `pulse_dur=0.00025`, `pulse_rnk=2`, `method='expm'`, `dead_time=0.00013`) and detection with `coil=state(spin_system,'L+','1H')`, `rho0=state(spin_system,'Lz','1H')`, `spins={'1H'}`, `offset=100`, `sweep=2000`, `npoints=8`.
- Compares zero irradiation (`homodec_pwr=0`) against complete removal of the `homodec_oper` and `homodec_pwr` fields, requiring agreement to `1e-10` and confirming that admission does not invent the optional `homodec_oper` field when it is absent.
- After the loop, asserts a nonzero signal (`norm(fids{2,1})>0.1`), a measurable acquisition-stage irradiation effect (`norm(fids{2,1}-fids{2,3})>0.01`), and that the Hilbert-admitted and native Liouville FIDs agree to `1e-7` for all three irradiation variants.

## Inputs and outputs

- **Inputs**: none. The function takes no arguments.
- **Outputs**: `result` — a regression test record created with `new_test_result('kernel/homodec_adapt', ...)`, accumulated through `test_true` and `test_close` assertions covering the noncommuting irradiation check, commutation operator conversion, no-op formalism handling, absent operator handling, zero-power equivalence, nonzero signal, irradiation effect, and Hilbert/native FID agreement.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_homodec_adapt.m)
