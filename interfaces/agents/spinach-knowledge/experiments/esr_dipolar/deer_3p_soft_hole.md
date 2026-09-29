# experiments/esr_dipolar/deer_3p_soft_hole.m

This is a pulse diagnostic for the three-pulse DEER/PELDOR experiment, not the DEER echo-stack calculation. The source describes a hypothetical test in which a selected soft pulse is followed by an ideal hard `pi/2` pulse and time-domain acquisition.

## Preparation and acquisition

Each of the three shaped pulses is applied independently to `parameters.rho0`; the three pulse responses are not composed sequentially. The unpulsed state is retained as a reference. The code assembles these four states, applies a hard `pi/2` rotation about the constructed `Ey` operator, and calls `acquire`. That operator is built for the first entry of `parameters.spins`; the header describes the hard pulse as acting on all spins. The acquisition uses the receiver offset, sweep, and point count supplied in `parameters`.

Required fields are `parameters.pulse_frq`, `parameters.pulse_pwr`, `parameters.pulse_dur`, `parameters.pulse_phi`, and `parameters.pulse_rnk` (three pulse values each), plus `parameters.offset`, `parameters.sweep`, `parameters.npoints`, `parameters.spins`, `parameters.rho0`, `parameters.coil`, and `parameters.method`. Pulse frequencies are Hz, powers rad/s, durations seconds, phases radians, and ranks integer Fokker-Planck ranks. Offset and sweep are Hz. The method is `expm`, `expv`, or `evolution`; context matrices `H`, `R`, and `K` must be same-sized. The spin list normally identifies electron spins.

The function returns `fids` from `acquire`. Its implementation passes four prepared states (reference plus three pulse-specific states), and the diagnostic wrapper consumes four FID columns. The header output note instead says three FIDs. The returned `fids` matrix has `parameters.npoints` time-sample rows and four state columns (the reference and three pulse responses); the source does not state numeric signal units. It describes the acquisition as infinite-bandwidth while also requiring a sweep value; the implementation passes that value to `acquire`, so the source does not resolve the wording further.

Source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/deer_3p_soft_hole.m
