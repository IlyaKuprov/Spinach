# etc/textbook/quad_shift.m

- Signature: `delta=quad_shift(Cq,eta,v0,S,m)`

## Purpose

Evaluates the second-order quadrupolar shift of the powder-pattern centre of gravity for the `|S,m⟩ → |S,m−1⟩` NMR transition of a quadrupolar nucleus.

## Expression

The implementation uses Samoson's expression (Equation 3 in the cited paper):

```matlab
delta = -1e6*(3/40)*(Cq/v0)^2*(1+eta^2/3) ...
        *(S*(S+1)-9*m*(m-1)-3)/(S^2*(2*S-1)^2);
```

The result is in ppm. The source notes that this expression was checked against numerical calculations.

## Inputs

- `Cq` — quadrupolar coupling constant in Hz; real scalar.
- `eta` — quadrupolar asymmetry parameter; real scalar.
- `v0` — nuclear Larmor frequency in Hz; real scalar.
- `S` — nuclear spin; integer or half-integer greater than 1/2.
- `m` — magnetic quantum number for an existing `|S,m⟩ → |S,m−1⟩` transition.

## Output

- `delta` — second-order quadrupolar shift in ppm.

## Reference

Samoson, Equation 3: [original article](https://doi.org/10.1016/0009-2614(85)85414-2). See also the [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=quad_shift.m).
