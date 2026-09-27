# etc/textbook/r2csa2tauc.m

- Signature: `tauc=r2csa2tauc(R2,del_sq,B0,isotope)`

## Purpose

Calculates rotational correlation-time solutions consistent with a transverse CSA relaxation rate under the model implemented by this function. The Zeeman angular frequency is obtained from `B0` and the isotope's magnetogyric ratio.

## Method

The transverse CSA relation is rearranged as a cubic in the correlation time, and the function evaluates its three closed-form roots. It returns the roots in the `tauc` vector; entries may be complex. If all three are non-real, the function reports that there are no real solutions under this model.

## Inputs

- `R2` — transverse relaxation rate in Hz; positive real scalar.
- `del_sq` — positive real scalar, the second-rank CSA invariant (see `blinv`).
- `B0` — real scalar magnetic field in tesla.
- `isotope` — character array identifying the isotope (for example, `'1H'`).

## Output

- `tauc` — vector of three algebraic correlation-time solutions in seconds.

## Reference

See the [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=r2csa2tauc.m).
