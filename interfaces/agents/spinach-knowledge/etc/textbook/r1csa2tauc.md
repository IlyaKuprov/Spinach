# etc/textbook/r1csa2tauc.m

- Signature: `tauc=r1csa2tauc(R1,del_sq,B0,isotope)`

## Purpose

Finds the rotational correlation times consistent with a longitudinal CSA relaxation rate under the model implemented by this function. The Zeeman angular frequency is obtained from `B0` and the isotope's magnetogyric ratio.

## Method

The CSA relaxation relation is rearranged as a quadratic in the correlation time. The function evaluates the larger root first to reduce cancellation, then obtains the other root from the product of the roots. It returns both candidate times in the `tauc` vector, with the larger candidate first. If neither solution is real, it reports that the input case is physically impossible under this model.

## Inputs

- `R1` — longitudinal relaxation rate in Hz; positive real scalar.
- `del_sq` — positive real scalar, the second-rank CSA invariant (see `blinv`).
- `B0` — real scalar magnetic field in tesla.
- `isotope` — character array identifying the isotope (for example, `'1H'`).

## Output

- `tauc` — two candidate rotational correlation times in seconds, ordered with the larger first.

## Reference

See the [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=r1csa2tauc.m).
