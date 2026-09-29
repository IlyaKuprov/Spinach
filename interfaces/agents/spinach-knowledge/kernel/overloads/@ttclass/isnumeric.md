# kernel/overloads/@ttclass/isnumeric.m

[Mapped MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/isnumeric.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ttclass/isnumeric.m)

## Signature

`answer=isnumeric(tt)`

## Behaviour

Returns logical true exactly when `tt` is a `ttclass` object and its `cores` field is non-empty. Otherwise it returns logical false. The test is on the core container itself; it does not inspect core values, the coefficient, or whether the represented tensor is nonzero.

## Input and output

- `tt` — tensor-train object.
- `answer` — logical predicate result.

This method does not change the train, its shape, ranks, or coefficient and does not materialise the represented tensor.
