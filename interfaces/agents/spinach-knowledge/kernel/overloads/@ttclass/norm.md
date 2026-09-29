# kernel/overloads/@ttclass/norm.m

- Signature: `ttnorm=norm(ttrain,norm_type)`
- Source: [`kernel/overloads/@ttclass/norm.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/norm.m)
- Wiki: [`ttclass/norm.m`](https://spindynamics.org/wiki/index.php?title=ttclass/norm.m)

## Core action

Only `norm_type='fro'` is computed. The method calls `pack(ttrain)`, then `ttort(ttrain,-1)`, and returns `abs(ttrain.coeff) * norm(ttrain.cores{1,1}(:),2)`. Thus it reduces the packed and orthogonalised train to the Euclidean norm of the flattened first core, with the coefficient contributing by its magnitude; it does not form a dense represented matrix. MATLAB's vector 2-norm uses modulus for complex entries, so this is the Hermitian Frobenius norm rather than a bilinear sum.

Packing may combine buffered terms into a single train; orthogonalisation is applied before the core norm. The result is a nonnegative scalar, and the method does not return a train or its ranks.

## Supported inputs

The explicit cases for `1`, `inf`, and `2` each raise an error stating that norm is unavailable for `ttclass`. Any other value raises an unrecognised-norm-type error. There is no separate type or shape check in this method.
