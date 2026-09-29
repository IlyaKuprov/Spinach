# experiments/cp_contact_hard.m

- Signature: `contact_curve=cp_contact_hard(spin_system,parameters,H,R,K)`

## Purpose

A general rotating-frame cross-polarisation contact calculation. It starts from `parameters.rho0`, applies an ideal 90-degree excitation assembled from user-supplied operators, then propagates under the context Liouvillian plus the channel-specific spin-lock fields. Unlike the soft-pulse contact routine, the excitation pulse is idealised rather than given a finite duration.

## Inputs and timing

- `parameters.exc_opers`: cell array of excitation spin operators; the implementation sums them and applies the same ideal flip angle to all channels.
- `parameters.irr_opers`: cell array of spin-lock operators, one per RF channel.
- `parameters.irr_powers`: nutation frequencies in Hz, with one row per channel and one column per contact time slice; row count must match `irr_opers`.
- `parameters.time_steps`: vector of slice durations in seconds. The fields determine a piecewise-constant RF schedule; no uniform spacing is assumed by the function.
- `parameters.rho0`: initial state; `parameters.coil`: detection state. `parameters.spc_dim` is also required and checked as a positive integer for operator embedding.
- `H`, `R`, and `K`: same-sized Hamiltonian, relaxation, and kinetics matrices supplied by the context function; the routine forms `H+1i*R+1i*K`.

## Output

`contact_curve` contains the coil expectation at the initial state and after each RF slice. Its source-defined dimensions are `size(parameters.coil,2)` by `numel(parameters.time_steps)+1`; the first column is the pre-evolution value. During each slice, the Hamiltonian receives the sum of the channel terms `2*pi*irr_powers(channel,slice)*irr_opers{channel}`.

## Source limits

The sequence is operator-driven and does not assign isotope identities, pulse phases, or an experimental CP calibration; those choices are represented by the supplied operators and powers.

Source implementation: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/cp_contact_hard.m
