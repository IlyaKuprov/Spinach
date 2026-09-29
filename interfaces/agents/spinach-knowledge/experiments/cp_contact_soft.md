# experiments/cp_contact_soft.m

- Signature: `contact_curve=cp_contact_soft(spin_system,parameters,H,R,K)`

## Purpose

A two-channel rotating-frame CP contact with a finite-duration high-gamma excitation pulse and observable detection through the contact. The routine wipes the low-gamma component of `parameters.rho0`, excites the high-gamma spin along +X, and applies high-gamma -Y and low-gamma +X irradiation during the contact.

## Inputs and timing

- `parameters.spins`: two isotope names, high-gamma first and low-gamma second (for example, `{'1H','13C'}`).
- `parameters.rho0`: initial state; the low-gamma spin state is wiped before the pulse.
- `parameters.hi_pwr`: high-gamma excitation nutation frequency in Hz; the +X 90-degree pulse duration is `1/(4*hi_pwr)` seconds.
- `parameters.cp_pwr`: two CP-channel nutation frequencies in Hz, ordered by `parameters.spins`.
- `parameters.timestep`: CP integration step in seconds; `parameters.nsteps`: number of steps, so the specified contact integration spans `timestep*nsteps` seconds.
- `parameters.coil`: detection state. `H`, `R`, and `K` are same-sized Hamiltonian, relaxation, and kinetics matrices supplied by the context function; the routine forms `H+1i*R+1i*K`.

The implementation also reads `parameters.spc_dim` when embedding control operators. The parameter header does not document this field's meaning or units.

## Output

Returns `contact_curve`, the observable-mode result detected on `parameters.coil`: `nsteps+1` time rows, beginning with the post-excitation state before contact and continuing with one row after each contact step. The step spacing is `parameters.timestep`; multiple detection states, if supplied, occupy separate columns.

## Source limits

This routine computes only the CP contact and its detection. It does not append an FID acquisition or specify instrument-specific RF calibration or a matching-condition search.

Source implementation: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/cp_contact_soft.m
