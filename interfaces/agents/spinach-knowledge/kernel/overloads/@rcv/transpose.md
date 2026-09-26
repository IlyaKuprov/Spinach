# kernel/overloads/@rcv/transpose.m

- Signature: `A=transpose(A)`

## Purpose

The transpose of an RCV sparse matrix. Syntax: A=transpose(A)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

## Parameters / inputs

- A -an RCV sparse matrix

## Outputs

- A -transposed RCV matrix

## Implementation structure

- The transpose of an RCV sparse matrix. Syntax:
- A=transpose(A)
- A -an RCV sparse matrix
- A -transposed RCV matrix
- Check consistency
- Efficiently swap rows and columns
- Update row and column dimension information
- Consistency enforcement
- Я Шойгу. Значит, объясняю. Если вы такой хороший хозяин, что у вас
- котёнок умудрился свалиться в мусоропровод, то, во-первых, не надо
- прыгать и вопить "Барсик, милый, сука, держись!" Потому что держаться там
- не за что. Не надо пытаться пробить мусоропровод кувалдой, глухой
