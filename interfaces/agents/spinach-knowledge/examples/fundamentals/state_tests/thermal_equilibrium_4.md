# examples/fundamentals/state_tests/thermal_equilibrium_4.m

- Signature: `thermal_equilibrium_4()`
- Source: [`examples/fundamentals/state_tests/thermal_equilibrium_4.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/state_tests/thermal_equilibrium_4.m)

## Purpose

Checks whether a computed thermal-equilibrium state is stationary under two thermalised relaxation superoperators: the Dibari–Levitt method and the inhomogeneous master equation (IME) route.

## Model and relaxation settings

The model has four `19F` spins, a magnet-field parameter of `9.4`, and scalar Zeeman parameters `-120.5380`, `-133.9429`, `-129.3169`, and `-129.5320`. Its six scalar-coupling entries are: pairs 1–2 and 3–4, `271.2924` each; 1–3, `0.5401`; 1–4, `-25.9884`; 2–3, `9.9625`; and 2–4, `-40.7675`. The source does not state units for these parameters.

Relaxation is configured as `damp`, with temperature parameter `40`, damping rate `5.0`, and `rlx_keep='labframe'`; the equilibrium setting starts at `zero` before the explicit thermalisation calls. The basis uses `approximation='none'` for both `sphten-liouv` and `zeeman-liouv`. Full lab-frame retention is required for the `zeeman-liouv` branch: `relaxation` explicitly rejects `rlx_keep='diagonal'` in that formalism.

## Procedure and check

For each formalism the example builds `rho_eq` and the relaxation superoperator `R`. It constructs the lab-frame left-action Hamiltonian for the Dibari–Levitt call to `thermalize`, then also obtains an IME-thermalised superoperator. After each method, it checks `norm(Rt*rho_eq,2)` against `1e-9` and raises an error above that threshold.

## Output and limits

The observable is stationarity of the supplied equilibrium state under each thermalised operator; this example does not produce an acquisition or spectrum. The source includes a success message, but no execution result is claimed here. The numerical settings are reproduced as written because the source does not declare their units.
