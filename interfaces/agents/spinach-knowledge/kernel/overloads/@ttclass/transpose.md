# kernel/overloads/@ttclass/transpose.m

## Signature

`ttrain=transpose(ttrain)`

## Behaviour

For every train and core, the method permutes the four core dimensions in the order `[1 3 2 4]`. This leaves the left and right bond dimensions and the core sequence unchanged while swapping the two physical dimensions (row and column modes). The resulting tensor train represents the non-conjugating matrix transpose: complex entries are not conjugated, and the represented matrix dimensions are exchanged.

The core ranks, train coefficients, and tolerance metadata are unchanged. This is a dimension permutation only; no rank reduction, truncation, or rounding is performed.

## References

- [Source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/transpose.m)
- [Spin Dynamics Wiki: `ttclass/transpose.m`](https://spindynamics.org/wiki/index.php?title=ttclass/transpose.m)
