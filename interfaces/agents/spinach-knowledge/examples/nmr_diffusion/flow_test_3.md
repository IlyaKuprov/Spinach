# examples/nmr_diffusion/flow_test_3.m

- MATLAB implementation: [examples/nmr_diffusion/flow_test_3.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_diffusion/flow_test_3.m)

- Signature: `flow_test_3()`

## Purpose

This example transports a signed concentration-like field through a three-dimensional circular flow while it diffuses. The Spinach system is deliberately a ghost spin with empty Zeeman and coupling matrices, so the calculation is a transport demonstration rather than an NMR signal simulation. The source estimates minutes, with faster execution on a GPU.

## Spin system and transport model

The system uses isotope G, a field setting of 5.9, spherical-tensor Liouville formalism, and no basis approximation. Although a field is configured, the empty interaction matrices leave no spin evolution in this example. A cubic 0.02 m box is sampled on a 50 × 50 × 50 grid. The derivative setting is `{'period',7}`; the source does not attach a unit to 7.

At grid coordinates X, Y, and Z, the coded circular-flow field is `u=-1000*Y`, `v=1000*X`, and `w=0`. The diffusion tensor is spatially uniform and isotropic: its three diagonal entries are 8×10⁻⁶ and its off-diagonal entries are zero. The grid coordinates are in metres; the source does not annotate a unit for the diffusion entries or flow coefficient.

## Initial field and propagation

The initial field is the sum of three Gaussian-shaped peaks, with signs positive, negative, positive and centres at (−0.003, −0.003, −0.003), (0.001, 0.001, 0.001), and (0.003, −0.003, 0.003) m. The source sets `sigma=2e-6` in the Gaussian expression, without stating its unit. The Fokker–Planck generator is made with `v2fplanck` and inflated before `evolution` propagates the flattened field with the coded step parameter 5×10⁻⁵ for 200 steps.

## Output and scope

The routine opens a figure and sends each propagated three-dimensional field to `volplot`, updating the view through the trajectory. It does not calculate a receiver signal, spectrum, or spin observable, and the script itself supplies no saved numerical result or convergence study; the page describes the configured example, not a measured run.
