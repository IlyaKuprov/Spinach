# kernel/utilities/wigner_3j.m

## Purpose

Calculates Wigner 3j-symbols of the form

```
/ j1 j2 j3 \
\ m1 m2 m3 /
```

using the Spinach library's Clebsch-Gordan coefficient routine. Source: [kernel/utilities/wigner_3j.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/wigner_3j.m).

## Behaviour

- Syntax: `w = wigner_3j(j1,m1,j2,m2,j3,m3)`.
- The function first runs an internal consistency check (`grumble`) on all six arguments, then converts the requested 3j-symbol to a Clebsch-Gordan coefficient via

```
w = ((-1)^(-m3+j1+j2))/sqrt(2*j3+1) * clebsch_gordan(j3,-m3,j1,m1,j2,m2)
```

- If physically inadmissible indices are supplied, a zero is returned (per the header documentation).
- The consistency check raises errors when: any argument is non-numeric; any argument is not real; any argument does not have exactly one element; or any argument is not integer or half-integer (tested via `mod(2*x+1,1)~=0`).

## Inputs and outputs

**Inputs**

- `j1`, `j2`, `j3` - integers (or half-integers, per the validation logic) arranged as the top row of the 3j-symbol.
- `m1`, `m2`, `m3` - integers (or half-integers, per the validation logic) arranged as the bottom row of the 3j-symbol.

**Outputs**

- `w` - the resulting 3j-symbol.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=wigner_3j.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/wigner_3j.m>
