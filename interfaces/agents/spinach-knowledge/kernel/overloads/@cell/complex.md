# kernel/overloads/@cell/complex.m

- Signature: `A=complex(A)`

For each linear cell index, the overload replaces the content with the result of `complex(A{n})`. The cell array's indexing and shape are retained; conversion behaviour for each value is delegated to MATLAB's `complex` dispatch. The documented use is a cell array of numeric objects. This wrapper has no explicit cell-type or element-type validation and adds no broadcasting or two-argument real/imaginary construction.

## References

- MATLAB source: [`kernel/overloads/@cell/complex.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@cell/complex.m)
- Spinach Wiki: [`cell/complex.m`](https://spindynamics.org/wiki/index.php?title=cell/complex.m)
