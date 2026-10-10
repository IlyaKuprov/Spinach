# examples/spin_chemistry/cidnp_flash_acquire.m

Source: [examples/spin_chemistry/cidnp_flash_acquire.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/spin_chemistry/cidnp_flash_acquire.m)

## Purpose

Simulate CIDNP magnetisation pumping during illumination, followed by dark relaxation. The source identifies the model with Ilya Kuprov's paper at https://doi.org/10.1016/j.jmr.2004.01.011.

## Spin system and relaxation

The model uses 1H and 19F at 14.1 T, zero isotropic shifts, and a 50 Hz scalar coupling. The fluorine chemical-shift-anisotropy principal values are entered as [-47, -16, 63], with zero Euler angles; the two coordinates are [0, 0, 0] and [0, 2.60, 0]. Redfield relaxation uses secular retention, zero equilibrium, and a 110 ps correlation time.

The script normalises the 1H longitudinal, 19F longitudinal, and two-spin longitudinal-product states. It obtains the NMR Hamiltonian and relaxation matrix, then adjusts the matrix's identity-state column to thermalise to the sum of the two longitudinal states. Two magpump terms add light-driven contributions, with the source parameters 1.3 for 1H and 34.0 for 19F; the source does not label units for these parameters. The initial state is the unit state plus both longitudinal magnetisations.

## Pump and acquisition schedule

The illuminated trajectory uses 50 steps of 0.01 s, giving a 0.5 s pump interval with both pumping terms active. The subsequent trajectory starts from its endpoint and uses the ordinary relaxation matrix for 500 steps of 0.01 s, a further 5 s. The concatenated time axis spans 0 to 5.5 s in 0.01 s increments.

## Observables

The plotted channels are the 19F longitudinal magnetisation, the 1H longitudinal magnetisation, and minus twice the longitudinal product state. These are simulated state projections plotted against time; the source-only description does not claim measured polarisation or a validated fit.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
