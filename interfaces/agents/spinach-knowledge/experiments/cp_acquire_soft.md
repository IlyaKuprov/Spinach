# experiments/cp_acquire_soft.m

- Signature: `fid=cp_acquire_soft(spin_system,parameters,H,R,K)`

## Purpose

A rotating-frame cross-polarisation transfer followed by time-domain FID acquisition. The source describes the first spin as high-gamma and the second as low-gamma (for example, proton then carbon): it wipes the low-gamma component of the starting state, excites the high-gamma channel, applies a two-channel CP contact, then wipes the high-gamma state and acquires while decoupling that channel.

## Inputs

- `parameters.spins`: two isotope names, high-gamma first and low-gamma second (for example, `{'1H','13C'}`).
- `parameters.rho0`: initial state; the low-gamma spin state is wiped before the transfer.
- `parameters.hi_pwr`: high-gamma excitation nutation frequency in Hz. The routine applies a +X 90-degree excitation of duration `1/(4*hi_pwr)` seconds.
- `parameters.cp_pwr`: two nutation frequencies in Hz, ordered by `parameters.spins`, for the CP contact. During contact, the high-gamma channel is irradiated along -Y and the low-gamma channel along +X.
- `parameters.cp_dur`: contact duration in seconds.
- `parameters.coil`: detection state; `parameters.sweep` is FID sweep width in Hz and `parameters.npoints` is its positive integer point count.
- `H`, `R`, and `K`: same-sized Hamiltonian, relaxation, and kinetics matrices supplied by the context function; the routine combines them as `H+1i*R+1i*K`.

The implementation also reads `parameters.spc_dim` to embed the channel operators in the full space. This field is not listed in the function's parameter header, which does not define its meaning or units.

## Output

Returns `fid`, the observable-mode FID from evolution with the high-gamma channel decoupled after CP. The call uses the reciprocal sweep width as the sampling interval and `npoints-1` evolution steps; the source does not prescribe a reshaping or row/column convention for the returned array.

## Source limits

The function fixes the two irradiation axes and uses nutation frequencies as RF amplitudes in its rotating-frame generator; it does not define instrument-specific RF calibration or a CP matching-condition search.

Source implementation: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/cp_acquire_soft.m
