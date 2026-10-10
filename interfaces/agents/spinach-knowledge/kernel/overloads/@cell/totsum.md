# kernel/overloads/@cell/totsum.m

- Signature: `S=totsum(A)`

## Behaviour

`totsum` reduces a cell array by adding all its numeric entries into one matrix; unlike the binary cell operators, this does combine entries. The code first checks that `isnumeric` is true for every cell member. It does not explicitly check that the cell is nonempty or that all member sizes match.

If every member is sparse, the implementation gathers each member's nonzero row indices, column indices, and values, concatenates those triplets, and constructs a sparse result sized to `size(A{1})`. Repeated row/column positions are combined by MATLAB's sparse constructor. There is no orientation weighting or tensor contraction: this is a sum over cell entries, with their matrix coordinates identifying contributions. The constructed dimensions come from the first entry, so all contributed indices must fit those dimensions.

If at least one member is not sparse, the code initialises `S` as zeros with the first member's size and class, then adds entries in sequence using MATLAB addition. Compatible array sizes are required by those additions, subject to MATLAB's elementwise expansion rules; the wrapper does not reconcile shapes. This branch does not use the sparse-triplet construction. It neither implements `opium`/`polyadic` methods nor defines their semantics.

## Inputs and validation

The intended input is a nonempty cell array whose members all pass `isnumeric`. The source checks member numeric status, but has no explicit nonempty or common-size check; it accesses the first member when constructing either result. Sparse-path dimensions are taken from that first member. Full-path compatibility is determined by the additions themselves.

## Output

A single matrix containing the sum over all entries. The all-sparse path constructs sparse storage; the other path uses a zeros accumulator initialised from the first entry and then follows MATLAB's addition behaviour.

## References

- Source: [`kernel/overloads/@cell/totsum.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@cell/totsum.m)
- Wiki: [cell/totsum.m](https://spindynamics.org/wiki/index.php?title=cell/totsum.m)
