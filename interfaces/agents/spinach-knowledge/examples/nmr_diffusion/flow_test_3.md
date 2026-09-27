# examples/nmr_diffusion/flow_test_3.m

- Signature: `flow_test_3()`

## Purpose

Shows a three-dimensional circular-flow advection–diffusion calculation with no active spin interactions. The source describes a minutes-long calculation, faster on GPU.

## Physical and numerical content

The system uses a ghost spin with empty Zeeman and coupling matrices. A 50 × 50 × 50 grid spans a 0.02 m cube; the velocity field is u = −1000y, v = 1000x, w = 0, and the diffusion tensor is isotropic with diagonal entries 8×10⁻⁶ and zero off-diagonal entries. The derivative settings are {period, 7}. The initial state is a signed combination of three Gaussian peaks centered at the coded coordinates, with sigma = 2×10⁻⁶.

The example builds the Fokker–Planck generator with `v2fplanck`, inflates it, and obtains a 200-point trajectory using `evolution` with a step parameter of 5×10⁻⁵. It plots each three-dimensional state with `volplot`.
