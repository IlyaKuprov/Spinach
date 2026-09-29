# kernel/overloads/@ttclass/isreal.m

[Mapped MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/isreal.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ttclass/isreal.m)

## Signature

`answer=isreal(tt)`

## Behaviour

For a `ttclass` input, the method first evaluates `all(isreal(tt.coeff))`. Only if that is true does it visit every stored core `tt.cores{k,n}` for `n=1:tt.ntrains` and `k=1:tt.ncores`; it returns early when a core is not real. Thus the predicate covers the stored coefficients and core entries, rather than forming or inspecting a materialised tensor. It applies no conjugation and changes no train data.

A non-`ttclass` input raises the error `input is not a ttclass.`.

## Input and output

- `tt` — tensor-train object.
- `answer` — logical result of the coefficient and core checks.
