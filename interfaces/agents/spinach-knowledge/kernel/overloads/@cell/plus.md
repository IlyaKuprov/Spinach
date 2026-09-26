# kernel/overloads/@cell/plus.m

- Signature: `C=plus(A,B)`

## Purpose

Adds cell arrays element-by-element. Syntax: A=plus(A,B)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- A,B -cell arrays of identical topology

## Outputs

- C -the resulting cell array

## Implementation structure

- Adds cell arrays element-by-element. Syntax:
- A=plus(A,B)
- A,B -cell arrays of identical topology
- C -the resulting cell array
- Check consistency
- Decide the topology
- Add cell-by-cell
- Add to each cell
- Complain and bomb out
- Consistency enforcement
- I came into the room, which was half dark, and presently spotted Lord
- Kelvin in the audience and realized that I was in for trouble at the last
