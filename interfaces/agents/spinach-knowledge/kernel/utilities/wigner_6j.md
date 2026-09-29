# kernel/utilities/wigner_6j.m

## Purpose

Computes the Wigner 6j-symbol

```
     / j1 j2 j3 \
     \ j4 j5 j6 /
```

for integer or half-integer angular momentum arguments, using a direct summation over products of four Wigner 3j-symbols.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/wigner_6j.m>

## Behaviour

- Syntax: `w = wigner_6j(j1,j2,j3,j4,j5,j6)`.
- The function first calls an internal consistency check (`grumble`) on all six arguments.
- The result is accumulated starting from `w = 0` over all combinations of magnetic quantum numbers `m1` through `m6`, each running from `-j` to `+j` in unit steps.
- For each combination, a sign factor `(-1)^power_of_minus` is computed, where `power_of_minus = (j1-m1)+(j2-m2)+(j3-m3)+(j4-m4)+(j5-m5)+(j6-m6)`.
- A screening condition restricts contributions to terms satisfying all four relations simultaneously:
  - `m1 + m2 - m3 == 0`
  - `-m1 + m5 + m6 == 0`
  - `m4 - m5 + m3 == 0`
  - `-m4 - m2 - m6 == 0`
- Screened terms add `(-1)^power_of_minus` times the product of four Wigner 3j-symbols:
  - `wigner_3j(j1, m1, j2, m2, j3, -m3)`
  - `wigner_3j(j1, -m1, j5, m5, j6, m6)`
  - `wigner_3j(j4, m4, j5, -m5, j3, m3)`
  - `wigner_3j(j4, -m4, j2, -m2, j6, -m6)`
- If physically inadmissible indices are supplied, a zero is returned.

## Inputs and outputs

**Inputs**

- `j1` ... `j6` — integers or half-integers arranged in the order

```
     / j1 j2 j3 \
     \ j4 j5 j6 /
```

Each argument must be numeric, real, and a scalar (exactly one element). Additionally, `2*j + 1` must be an integer for every argument, i.e. each must be integer or half-integer. Violations raise errors:

- `'all arguments must be numeric.'`
- `'all arguments must be real.'`
- `'all arguments must have one element.'`
- `'all arguments must be integer or half-integer.'`

**Outputs**

- `w` — the resulting 6j-symbol.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=wigner_6j.m>
