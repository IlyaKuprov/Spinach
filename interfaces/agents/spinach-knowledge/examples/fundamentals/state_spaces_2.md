# examples/fundamentals/state_spaces_2.m

- Signature: `state_spaces_2()`
- Source: [`examples/fundamentals/state_spaces_2.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/state_spaces_2.m)

## Model and question

This pulse-acquire proton NMR example displays the density-operator content by spin-correlation order for anti-3,5-difluoroheptane. It concerns correlation-order state-space representation, not an angular or powder-quadrature test. The source comment calls the molecule a 16-spin example; the isotope array has 23 sites: seven spin-zero `12C`, fourteen `1H`, and two `19F`, so the listed spin-bearing nuclei number 16.

The field is 11.7464 T. Chemical shifts and scalar couplings are specified explicitly in the source. The two `19F` shifts at sites 10 and 18 are set to zero, with source comments giving -184.1865 and explaining that zero is used because those values do not matter for this calculation and are faster. The basis uses `sphten-liouv`, `IK-0`, inter-level 1, manually populated projections at levels 1-3, two `S3` symmetry groups, longitudinal `19F` states, and projection `{1}`. Automatic `zte` state dropout is disabled. GPU enablement is commented out.

## Propagation and output

The initial state is proton `L+`; the NMR-assumption Hamiltonian is propagated in trajectory mode with a 1 ms step for 1000 steps. The source does not apply an explicit RF pulse with `step`; it starts from the transverse initial condition. `trajan(...,'correlation_order')` plots the contributions, with a logarithmic y axis and display range 1e-6 to 3. Those plot limits are not pass/fail criteria.

The source comments say nine- and ten-spin contributions appear in the lower part of the figure and suggest that correlations through order eight suffice for practical simulation without relaxation. This is a source-stated interpretation, not an independently verified bound or a general truncation rule. The code defines no acceptance tolerance or quadrature comparison. Runtime is estimated in the source as minutes, with a GPU noted as faster.
