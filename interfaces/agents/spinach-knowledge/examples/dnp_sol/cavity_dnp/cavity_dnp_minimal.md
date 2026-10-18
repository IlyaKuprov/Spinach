# cavity_dnp_minimal

`result=cavity_dnp_minimal(nlevels,max_rank,ncycles)` uses a cavity Fock dimension, pulse Fokker–Planck rank, and RF sweep count. The usual call is `cavity_dnp_minimal(3,2,300)`; abbreviated calls retain rank 2 and 300 sweeps.

The calculation retains all twenty Bi:Si donor levels and the cavity in `sphten-liouv`. It illustrates Fig. 1 and Appendix D of Lorriaux et al., arXiv:2610.19180v1, rather than fitting an ensemble measurement.

Thermal preparation at 25 mK is followed by RF saturation of the eight lower-manifold lines, a selective microwave pulse, a selective RF pulse, and a 20 s cavity reset. `reduce` removes only the disconnected cavity spectator during the pulses. `shaped_pulse_af` retains the donor drift, and `coherence` averages donor phase without removing electron–nuclear zero-quantum coherences. Vacuum attachment, spin-cavity exchange, photon loss, and propagation use Spinach's bosonic state, interaction, relaxation, and evolution APIs.

RF sweeps replace experimental band-limited noise; rectangular pulses replace WURST shapes. Other approximations are a disconnected cavity during pulses, a zero-temperature photon bath, illustrative 100 Hz single-donor coupling, and omission of intrinsic spin relaxation and direct nuclear magnetic drive.

The result contains nuclear Iz, photon number, trace, and Fock-edge trajectories, plus preparation diagnostics and objects used by `validate_cavity_dnp`. Observable duals use left-product operators; the trace uses `state(spin_system,'E',1)`, not a product of identities. The nuclear polarisation change is not a claim about twofold |6> population enhancement.
