# examples/fundamentals/state_tests/thermal_equilibrium_5.m

- Signature: `thermal_equilibrium_5()`
- Source: [`examples/fundamentals/state_tests/thermal_equilibrium_5.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/state_tests/thermal_equilibrium_5.m)

## Purpose

Compares thermodynamic-equilibrium recovery trajectories across the spherical-tensor and Zeeman Liouville formalisms, using both Dibari–Levitt and IME thermalisation settings.

## Model and settings

The model is the same four-`19F` spin system specified in the adjacent equilibrium-superoperator test: magnet-field parameter `9.4`, Zeeman scalar entries `-120.5380`, `-133.9429`, `-129.3169`, `-129.5320`, and scalar couplings `271.2924` (1–2 and 3–4), `0.5401` (1–3), `-25.9884` (1–4), `9.9625` (2–3), and `-40.7675` (2–4). Relaxation is `damp`, with temperature parameter `40`, damping rate `5.0`, and `inter.rlx_keep='labframe'`. Full lab-frame retention is required for the `zeeman-liouv` branch: `relaxation` explicitly rejects `rlx_keep='diagonal'` in that formalism. Units for these numerical parameters are not specified in this source.

The comparison crosses `sphten-liouv` and `zeeman-liouv` with the `dibari` and `IME` equilibrium methods; each basis uses `approximation='none'`. For each combination, it obtains the equilibrium state, applies a π-radian `19F` `Lx` rotation, and calculates an `Lz` observable recovery trajectory using the NMR Hamiltonian-plus-relaxation generator. The evolution call uses the source parameters `1e-3` and `1000`.

## Comparison and output

Four relative norm differences compare the trajectories between formalisms for each method and between the two methods for each formalism. The test's threshold is `1e-3`; exceeding it raises an error. Its output is this numerical cross-check, not a plotted spectrum. The source contains a success message, but no successful execution is claimed here.
