# experiments/imaging/press_1d.m

Signature: `fid=press_1d(spin_system,parameters,H,R,K,G,F)`

## Contract and sequence

This is an imaging callback, normally invoked as imaging(spin_system,@press_1d,parameters). The framework supplies H, R, K, F, gradient operators, and the spatially prepared initial state. The source accepts sphten-liouv or zeeman-liouv formalisms; H, R, K, and F must be same-sized matrices, and G must contain at least one gradient operator.

The code forms L = H + F + iR + iK and builds transverse RF controls Lx and Ly from L+ for parameters.spins{1}, replicated over the spatial grid. It applies the shaped pulse with +ss_grad_amp G{1}, then evolves under -ss_grad_amp G{1} for half the sum of the RF pulse durations to rephase the slice-select gradient. It then calls the common acquire routine with the selected state and coil. This implementation contains one slice-select/rephase stage before acquisition; the source does not apply a second refocusing pulse or a crusher gradient here.

## Parameters, units, and returned signal

The sequence requires spins (a nonempty cell array of spin labels), scalar positive-integer npts, numeric rho0, scalar ss_grad_amp, rf_frq_list, rf_amp_list, rf_dur_list, rf_phi, and positive-integer max_rank. The three RF lists must have the same number of entries. The source documents RF frequency in Hz, RF amplitude in rad/s, and pulse duration in seconds; ss_grad_amp is in T/m. max_rank controls the Fokker-Planck pulse operator; the source comment says 2 is usually enough. The downstream acquire call also needs sweep, npoints, coil, and decouple in the imaging parameter context. sweep is in Hz and acquire samples at 1/sweep-second spacing.

The returned fid is the one-dimensional complex FID from acquire, with `parameters.npoints` samples (a `parameters.npoints×1` column vector), not a spatial image. The source does not report measured data.

## Source-backed configuration example

examples/imaging/press_1d_example.m configures six 1H spins on a 100-point, 0.30 m sample, with ss_grad_amp 30 mT/m, rf_frq_list -100 kHz, rf_amp_list 2 pi times 5 kHz, rf_dur_list 50 microseconds, rf_phi pi/2, and max_rank 3. The acquisition settings are sweep 1000 Hz and npoints 128. The example also defines spatial longitudinal-state phantoms and an L+ coil; these are configuration choices, not a reported spectrum or executed result.

## References

- [Spinach Wiki: press_1d.m](https://spindynamics.org/wiki/index.php?title=press_1d.m)
- [Canonical MATLAB source: experiments/imaging/press_1d.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/press_1d.m)
- [Example configuration: press_1d_example.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/press_1d_example.m)
