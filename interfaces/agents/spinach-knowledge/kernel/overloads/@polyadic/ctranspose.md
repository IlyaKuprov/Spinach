# kernel/overloads/@polyadic/ctranspose.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/ctranspose.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=polyadic/ctranspose.m)

`ctranspose(p)` applies MATLAB's conjugate transpose to each core factor in each term, retaining the term and factor cell ordering. It also forms the new prefix from the old suffix in reverse order, conjugate-transposes those entries, and forms the new suffix analogously from the old prefix. This is the adjoint rule for the composed representation: reversing the left/right operator product while adjointing each factor.

The result remains polyadic; the method does not materialise Kronecker products or contract factors. Each factor's row and column dimensions exchange, so the represented operator's overall row and column dimensions exchange when the stored products are compatible. There is no explicit dimension check, type guard, or singleton-expansion/broadcasting logic in this overload.