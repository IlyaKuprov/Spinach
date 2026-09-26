# experiments/hyperpol/masdnp.m

- Signature: `dnp=masdnp(spin_system,parameters)`

## Purpose

Magic angle spinning DNP simulation, returning the rotor period averaged steady state magnetization. This function takes a lot of inspiration from the code donated by Fred Mentink, please ci- te Fred's papers if you are using it. Syntax: dnp=masdnp(spin_system,parameters)

## Physical / mathematical content

- Hyperpolarisation experiment implementations. They propagate driven electron-nuclear systems under microwave irradiation, MAS, relaxation, and repetition until transient or steady-state observables are assembled.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- parameters.spins -the spins to microwave
- parameters.rate -spinning rate, Hz
- parameters.axis -spinning axis direction vector.
- parameters.max_rank -rotor discretization grid rank
- parameters.mw_pwr -microwave power, rad/s
- parameters.mw_frq -microwave frequency, Hz
- parameters.mw_time -microwave irradiation duration
- before the average magnetistion
- is computed, seconds
- parameters.grid -the name of the spherical avera-
- ging grid
- parameters.coil -detection state
- parameters.verbose -set this to 1 to enable diag-
- nostic output

## Outputs

- dnp -enhancement of the user-specified state relative to
- the thermal equilibrium
- Note: increase the rotor rank and the spherical grid size until
- the answer stops changing. You will likely need huge values
- for both parameters.
- Note: this function must be called directly, without a context
- wrapper.

## Implementation structure

- Magic angle spinning DNP simulation, returning the rotor period
- averaged steady state magnetization. This function takes a lot
- of inspiration from the code donated by Fred Mentink, please ci-
- te Fred's papers if you are using it. Syntax:
- dnp=masdnp(spin_system,parameters)
- parameters.spins - the spins to microwave
- parameters.rate - spinning rate, Hz
- parameters.axis - spinning axis direction vector.
- parameters.max_rank - rotor discretization grid rank
- parameters.mw_pwr - microwave power, rad/s
- parameters.mw_frq - microwave frequency, Hz
- parameters.mw_time - microwave irradiation duration
