# experiments/respiration.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/respiration.m
Spinach Wiki: https://spindynamics.org/wiki/index.php?title=respiration.m

## Purpose and model

`respiration` implements the RESPIRATION cross-polarisation method described by the Aarhus-group paper, DOI [10.1021/jz3000905](http://dx.doi.org/10.1021/jz3000905). It accepts an initial state `parameters.rho0`, applies the coded loop sequence, then returns an FID detected with `parameters.coil`. The initial state is caller-supplied; this function does not prescribe its preparation.

The propagation generator is `L=H+1i*R+1i*K`. `H`, `R` and `K` must be numeric matrices of the same dimensions. The two entries in `parameters.spins` select the operator species. The source example is `{'1H','13C'}`; both labels must occur in the system. It constructs `L+` operators for those entries and their X components.

## Loop and acquisition

For each of `parameters.nloops` iterations, the source propagates under `L+2*pi*2*rate*Hx` and `L-2*pi*2*rate*Hx`, each for `1/(2*rate)`. Here `parameters.rate` is a positive pulse-train rate in Hz. It then applies the source's ideal `theta` step under `Cx+Hx`. The parameter documentation calls `theta` an angle; the code passes it directly as the final argument to `step` and specifies no conversion or independent unit for it.

After the loops, the code decouples the first selected spin and acquires with `evolution` using the supplied coil state, dwell `1/sweep` (with `sweep` in Hz), and `npoints-1` evolution steps. The result is one FID trace of `npoints` samples; the source does not explicitly set row-versus-column orientation.

## Required parameters and output

- `sweep`: positive real scalar, Hz; `npoints`: positive integer.
- `rho0`: initial state; `coil`: detection state.
- `nloops`: positive integer; `theta`: real scalar (documented as the ideal pulse angle).
- `rate`: positive real scalar, Hz; `spins`: two-element cell array of isotope names present in the system (for example, `{'1H','13C'}`).
- `fid`: the acquired FID as seen by `coil`.

The loop description follows the implementation; no measured transfer efficiency or experimental result is claimed.

## References

- Aarhus-group paper: [DOI 10.1021/jz3000905](http://dx.doi.org/10.1021/jz3000905)
- [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=respiration.m)
