# kernel/overloads/@ttclass/revert.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/revert.m`
- Signature: `tt=revert(tt)`
- Total lines: 42

## Purpose

Applies a bit-revert permutation to a tensor train operator by reversing the core order and swapping bond indices. Syntax: tt=revert(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- tt -tensor train operator

## Outputs

- tt -tensor train operator with reversed core order

## Implementation structure

- Applies a bit-revert permutation to a tensor train operator by
- reversing the core order and swapping bond indices. Syntax:
- tt=revert(tt)
- tt -tensor train operator
- tt -tensor train operator with reversed core order
- Read sizes and ranks
- Swap bond indices
- Revert the train direction
- Asking for efficiency and adaptability in the same program is
- like asking for a beautiful and modest wife... we'll probably
- have to settle for one or the other.
- Gerald M. Weinberg, "The psychology of computer programming"
