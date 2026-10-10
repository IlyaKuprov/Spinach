# examples/spin_chemistry/cidnp_transport.m

A small radical-pair model with two electrons and a nucleus demonstrates transport of nuclear polarisation into a diamagnetic singlet-channel product. A triplet-channel spin-free sink tracks the other recombination channel. Matching carries the nucleus while tracing out the electrons. Effective electron offsets and hyperfine coupling are specified directly in rad/s, with recombination rates in inverse seconds; this is not a fitted molecular system.

The initial electron singlet is scaled to unit concentration. Constant-generator evolution returns species populations and an unweighted nuclear-product signal. The example reports concentration balance and the final product polarisation and plots both population and polarisation trajectories.

Detection and reference operator vectors explicitly use the `exact` method of the four-argument `coil_state` primitive.
