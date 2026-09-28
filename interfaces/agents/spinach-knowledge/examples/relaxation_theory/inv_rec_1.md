# examples/relaxation_theory/inv_rec_1.m

- Signature: `inv_rec_1()`

## Purpose

Simulates a simple inversion-recovery experiment and plots the proton longitudinal signal during recovery. The source reports a calculation time of seconds.

## Physical / mathematical content

The model is a single proton at 14.1 T with a scalar Zeeman interaction of 1.5. It uses the `t1_t2` relaxation model with `r1_rates={5.0}`, `r2_rates={5.0}`, secular retention, `dibari` equilibrium, and temperature 298. The initial state is the thermal equilibrium density operator, inverted by a π pulse about `Lx`; longitudinal magnetization is then detected.

## Numerical / algorithmic content

The script uses the complete `sphten-liouv` basis, forms the NMR Hamiltonian and relaxation superoperator, and evolves the post-pulse state under `L + 1i*R`. `evolution` returns the detected observable for 1000 steps of 1 ms each, covering one second; the plotted signal is its real part.

## Implementation structure

After defining the spin system and relaxation parameters, the code builds the basis, obtains `equilibrium`, defines the `Lz` proton coil, and constructs the static Liouvillian and `Lx` pulse operator. A π `step` inverts the state, then the evolution result is plotted against `linspace(0,1,1001)` with the longitudinal-signal axis label.
