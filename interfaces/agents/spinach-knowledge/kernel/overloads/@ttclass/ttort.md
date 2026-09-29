# kernel/overloads/@ttclass/ttort.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/ttort.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=ttclass/ttort.m)

## Signature

`[tt,lognrm]=ttort(tt,direct)`

## Purpose and core ordering

Orthogonalises each buffered train in `tt` independently. Core row `k` of column `n` is the `k`th core of train `n`, with mode sizes `sz(k,1)` and `sz(k,2)` and bond ranks `r(k)` and `r(k+1)`. The direction is selected by `direct`: `+1` sweeps from core 1 to core `d`; `-1` sweeps from core `d` back to core 1. Other values raise an error.

## QR sweeps

For `direct=+1`, each core `k=1:d-1` is reshaped to `[r(k)*sz(k,1)*sz(k,2),r(k+1)]` and thin-QR factorised. The `Q` factor is reshaped back into core `k`, and `R` is multiplied into the next core, which is reshaped as `[r(k+1),sz(k+1,1)*sz(k+1,2)*r(k+2)]`. The next bond rank becomes the number of columns of `Q`.

For `direct=-1`, core `k=d:-1:2` is reshaped to `[r(k),sz(k,1)*sz(k,2)*r(k+1)]`; the preceding core is reshaped as `[r(k-1)*sz(k-1,1)*sz(k-1,2),r(k)]`. Thin QR is applied to the non-conjugate transpose of core `k`; the transposed `Q` becomes core `k`, and `R.'` is absorbed into core `k-1`. The rank at bond `k` is updated to the number of columns of `Q`.

At each step the code measures `norm(R,2)` and, when nonzero, normalises `R` by that value before passing it to the adjacent core. It also normalises the final boundary core for each train. With one output, each removed norm factor is multiplied into that train's coefficient. With two outputs, `tt.coeff` is set to ones, `lognrm` is initialised from `log(tt.coeff)` before that reset, and each `log(norm(R,2))` and final-core log norm is accumulated into the corresponding entry. The log-output form is intended for tensor norms that may exceed the preserved numeric example `realmax()=1.7977e+308`.

## Inputs and outputs

- `tt` — a tensor-train object, possibly with multiple buffered trains.
- `direct` — `+1` for left-to-right or `-1` for right-to-left orthogonalisation.
- `tt` — the orthogonalised tensor-train object, with core and rank arrays updated in the selected direction.
- `lognrm` — when requested, one accumulated log value per buffered train.
