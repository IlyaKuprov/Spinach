# kernel/overloads/@polyadic/minus.m

Source: [kernel/overloads/@polyadic/minus.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/minus.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=polyadic/minus.m)

`polyadic` storage keeps a sum of unopened Kronecker products: each outer `cores` cell is one additive term, and its inner cells are that term's matrix factors. For example, `cores={{A,B,C},{D,E}}` represents `A kron B kron C + D kron E`. `prefix` and `suffix` hold matrix actions on the left and right.

## Behaviour

`minus(a,b)` returns `plus(a,(-1)*b)`. The scalar left multiplication `(-1)*b` scales the smallest-by-`numel` factor in each term of `b` and simplifies the polyadic; `plus` then adds that result to `a`. The method does not itself call `full`, but its delegated `plus` call simplifies immediately. `simplify` can merge adjacent eligible non-`opium` factors with `kron` and unwrap a lone factor from a buffer-free single-term result. Thus the addition is represented through factor terms, but simplification can materialise factor products; it is not uniformly lazy.

`minus` has no separate dimension or type check. The delegated `plus` check rejects different row or column sizes when both operands are non-scalar; if one operand is scalar, that pairwise size check is skipped. This is not a general array-broadcasting rule. The mapped methods introduce no complex conjugation.

## Inputs and output

- `a` and `b`: polyadic objects, as specified by the source help text.
- `c`: the result of delegated addition and simplification; a trivial one-term, buffer-free result may be unwrapped to its underlying matrix.
