# kernel/operators/boson_mono.m

- Signature: `B=boson_mono(nlevels)`

## Purpose

Bosonic monomial operators of the following structure: B(k,q)=(Cr^k)*(An^q) obeying the following commutation relations with the po- pulation number operator N: [N,B(k,q)]=(k-q)*B(k,q)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

## Syntax

```matlab
B=boson_mono(nlevels)
```

## Parameters / inputs

- nlevels -number of bosonic ladder population
- levels, k and q go from 0 to nlevels-1

## Outputs

- B -a cell array with the following numbering
- map between (k,q) and a single index:
- (0,0)(0,1)(0,2) (1)(3)(6)
- (1,0)(1,1)(1,2) <=> (2)(5)(8)
- (2,0)(2,1)(2,2) (4)(7)(9)

## Implementation structure

- Bosonic monomial operators of the following structure:
- B(k,q)=(Cr^k)*(An^q)
- obeying the following commutation relations with the po-
- pulation number operator N:
- [N,B(k,q)]=(k-q)*B(k,q)
- B=boson_mono(nlevels)
- nlevels -number of bosonic ladder population
- levels, k and q go from 0 to nlevels-1
- B -a cell array with the following numbering
- map between (k,q) and a single index:
- (0,0)(0,1)(0,2) (1)(3)(6)
- (1,0)(1,1)(1,2) <=> (2)(5)(8)
