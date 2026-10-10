# kernel/overloads/@ttclass/mldivide.m

Direct mapped source: [kernel/overloads/@ttclass/mldivide.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/mldivide.m) · [existing Spinach Wiki entry](https://spindynamics.org/wiki/index.php?title=ttclass/mldivide.m)

## Behaviour

This overload computes tensor-train solution `x` from matrix `A` and vector `y`; both inputs must be `ttclass` objects. The method shrinks `A` and `y`, forms and shrinks `A'*A` and `A'*y`, then calls `amensolve(AA,Ay,1e-6)`. This is the normal-equation system `A'*A*x=A'*y`, not a direct unsymmetrised solve. MATLAB's `'` is conjugate transpose, so the products conjugate complex entries of `A` as part of the adjoint.

The documented roles are a tensor-train matrix `A` and vector `y`. The method's explicit guard checks only that both are `ttclass`; it has no separate mode-size, rank, or vector-shape check. Compatibility is exercised by the products and solver. The returned `x` is a tensor train produced by AMEn; its rank is solver-controlled, with fixed requested tolerance `1e-6`. The method does not convert the result to a dense array or explicitly shrink `x` afterward.

No numeric example or DOI is present in the source or either existing page.
