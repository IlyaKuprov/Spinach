# kernel/overloads/@ttclass/kron.m

[Mapped MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/kron.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ttclass/kron.m)

## Signature

`c=kron(a,b)`

## Core and coefficient action

The method shrinks both operands first. If shrinking does not leave a `ttclass` object, it converts that operand with `truncate(ttort(pack(x),+1))`, retaining scalar-topology trains as a single normalised train. It then requires equal core counts; a mismatch raises `tensor train structures must be consistent.`.

For each core index `k`, it reshapes each operand core using that train's ranks and two mode sizes, takes `kron(core_of_b,core_of_a)`, permutes the resulting axes, and reshapes to one output core. Its shape is `[a_ranks(k,1)*b_ranks(k,1), b_sizes(k,1)*a_sizes(k,1), b_sizes(k,2)*a_sizes(k,2), a_ranks(k+1,1)*b_ranks(k+1,1)]`. Thus each of the two mode sizes and each boundary rank is the product of the corresponding operand values; the rank ordering is the product of the two operands' ranks.

The output coefficient is `b.coeff*a.coeff`. If either coefficient is zero, output tolerance is set to zero; otherwise it is `abs(a.coeff)*b.tolerance+abs(b.coeff)*a.tolerance`. The method does not conjugate either operand. It returns a `ttclass` with constructed cores, not a full materialised matrix. Its element ordering differs from the flat matrix Kronecker product by row and column permutations, and is consistent with [`ttclass/vec`](https://spindynamics.org/wiki/index.php?title=ttclass/vec.m).

## Inputs and output

- `a`, `b` — tensor-train objects with equal numbers of cores.
- `c` — tensor-train Kronecker result.
