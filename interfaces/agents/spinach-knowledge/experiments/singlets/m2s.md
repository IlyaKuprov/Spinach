# experiments/singlets/m2s.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/singlets/m2s.m
Spinach Wiki: https://spindynamics.org/wiki/index.php?title=m2s.m

## Purpose and state model

`m2s` implements the M2S sequence of Pileio and Levitt. It receives the initial state vector `rho` from the caller and returns the final state vector; the source does not construct a particular magnetisation or coherence input. The steps below describe the programmed transformation, not a measured or guaranteed singlet yield.

The caller supplies the background Liouvillian `L`, X and Y spin operators `Hx` and `Hy`, and scalar `J` and `delta_v`. The source labels both `J` (J-coupling) and `delta_v` (Zeeman frequency difference) in Hz; each must be finite, real, scalar, and non-zero. `L`, `Hx` and `Hy` must be same-size matrices, and `rho` must have a row dimension matching `L`'s column dimension.

## Programmed sequence

The code sets `t=1/(4*sqrt(J^2+delta_v^2))` and computes `m1=floor(pi*abs(J)/(2*abs(delta_v)))`; if `m1` is odd it increments it to the next even integer. No additional unit conversion is performed.

1. Apply `pi/2` about `Hy` to the supplied `rho`.
2. Repeat `m1` times: evolve under `L` for `t`, apply `pi` about `Hx`, then evolve under `L` for another `t`.
3. Apply an `Hx` pulse of `sign(J)*pi/2` (90° in magnitude), then evolve under `L` for `t`.
4. Repeat `m1/2` times the same `L`-`Hx`-`L` block from step 2.

There is no acquisition axis in this routine: its output is the final state vector `rho`. The source does not define an input preparation stage or separately assert the state composition at intermediate steps.

## Reference

- [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=m2s.m)
