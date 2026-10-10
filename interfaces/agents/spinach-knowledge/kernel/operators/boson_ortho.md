# kernel/operators/boson_ortho.m

Direct source: [kernel/operators/boson_ortho.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/boson_ortho.m)

- Signature: `B=boson_ortho(nlevels)`

## Purpose

Return the `boson_mono(nlevels)` monomials after sequential Gram–Schmidt subtraction with respect to the Frobenius inner product implemented by `hdot`. The source processes cells in their inherited `boson_mono` order. For each current operator `B{n}` and each preceding `B{k}`, it subtracts `B{k}*hdot(B{k},B{n})/hdot(B{k},B{k})`. Here `hdot(X,Y)=sum(conj(X).*Y,'all')`.

## Ordering, dimensions, and normalisation

The result retains the input cell ordering and shape: an `nlevels^2-by-1` cell array of `nlevels-by-nlevels` matrices. The subtractions make each operator orthogonal to the preceding processed operators under `hdot`. There is no unit-norm rescaling: the routine does not divide each resulting matrix by its own norm, and the first item is unchanged.

## Operator action

Each item remains an operator matrix. The function performs no left/right action, Liouville-space lifting, commutator construction, or propagation.

## Reference

- [Spin Dynamics documentation for `boson_ortho.m`](https://spindynamics.org/wiki/index.php?title=boson_ortho.m)
- Related entry: [`boson_mono.m`](boson_mono.md)
