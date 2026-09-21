# kernel/overloads/@ttclass/transpose.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/transpose.m`
- Signature: `ttrain=transpose(ttrain)`
- Total lines: 38

## Purpose

Transposes a tensor without complex conjugation. Syntax: ttrain=transpose(ttrain)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- ttrain -tensor train representation of a matrix

## Outputs

- ttrain -transpose of the input tensor train

## Implementation structure

- Transposes a tensor without complex conjugation. Syntax:
- ttrain=transpose(ttrain)
- ttrain -tensor train representation of a matrix
- ttrain -transpose of the input tensor train
- Read tensor sizes and ranks
- Swap the middle dimensions of all cores
- "Public welfare" is the welfare of those who do not earn
- it; those who do, are entitled to no welfare.
- Ayn Rand, "Atlas Shrugged"
- #NGRUM
