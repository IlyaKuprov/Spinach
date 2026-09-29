# experiments/imaging/udd_dec.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/udd_dec.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=udd_dec.m)

## Purpose

`udd_dec` applies an Uhrig dynamic-decoupling (UDD) pulse train to an imaging phantom, then projects a user-selected spin state to form a sample image. It is a parameterised pulse-sequence simulation, not a measured phantom result. The source uses `H`, `R`, `K`, and `F` from the `imaging()` context to form `B=H+F+1i*R+1i*K`.

## Inputs and units

- `parameters.dec_time`: total sequence duration, in the simulation time unit (seconds in the Spinach imaging convention).
- `parameters.npulses`: number of UDD refocusing pulses, excluding the initial 90° pulse.
- `parameters.spins`: cell array of spin-name strings; the first entry selects the pulse operators (the source example is `{'1H'}`).
- `parameters.rho0`: initial state, supplied by the imaging setup.
- `parameters.npts`: imaging sample-grid dimensions.
- `parameters.coil_st`: cell array whose first entry is the user-specified detection state.

The sequence starts with a 90° y pulse, obtains intervals from `uhrig_times(dec_time,npulses)`, alternates each delay with a 180° x pulse, and applies the final delay. The pulse operators are built from `L+` for the selected first spin. The observed state is explicitly chosen through `coil_st{1}`; this function ignores the coil phantom rather than detecting through the imaging coil. It does not apply an explicit coherence-order filter.

## Detection and return

`mri` is the projected amplitude of the selected state at the sample points, produced by `fpl2phan(rho,parameters.coil_st{1},parameters.npts)`. Its spatial axes/grid are those encoded by `npts`; the function returns the phantom array without a separate coordinate vector. The timing and pulse count define a simulated design; no run-verified image or experimental measurement is claimed.
