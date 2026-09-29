# kernel/overloads/@ttclass/sum.m

## Signature

`answer=sum(ttrain,dim)`

## Behaviour

The method sums a tensor-train matrix over one matrix dimension: `dim=1` sums rows and `dim=2` sums columns, as specified by the function interface. If `dim` is omitted, it follows MATLAB's first-nonsingleton-dimension convention for the two matrix dimensions: dimension 1 when the row size is nonsingleton, otherwise dimension 2. An input whose matrix dimensions are both singleton is materialised immediately as a scalar.

For each core of each train, the selected physical mode is summed locally, then reshaped back into a core with that physical dimension set to one: `[left-rank,1,column-mode,right-rank]` for `dim=1`, or `[left-rank,row-mode,1,right-rank]` for `dim=2`. Core order, bond ranks, train coefficients, and the unsummed physical modes are retained. The output is the corresponding row- or column-shaped tensor-train matrix; if all its physical modes are singleton, it is materialised as a scalar instead.

No rank reduction, truncation, or rounding is performed. The output tolerance field is explicitly reset to zero for every train, rather than propagated from the input; this routine uses no tolerance-based approximation. The supported summation dimensions are 1 and 2.

## References

- [Source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/sum.m)
- [Spin Dynamics Wiki: `ttclass/sum.m`](https://spindynamics.org/wiki/index.php?title=ttclass/sum.m)
