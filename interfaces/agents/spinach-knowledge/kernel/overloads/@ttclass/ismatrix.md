# kernel/overloads/@ttclass/ismatrix.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/ismatrix.m`
- Signature: `answer=ismatrix(tt)`
- Total lines: 34

## Purpose

Returns TRUE for non-empty tensor train objects. Syntax: answer=ismatrix(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- tt -tensor train object

## Outputs

- answer -logical true for non-empty tensor train objects

## Implementation structure

- Returns TRUE for non-empty tensor train objects. Syntax:
- answer=ismatrix(tt)
- tt -tensor train object
- answer -logical true for non-empty tensor train objects
- Non-empty tensor trains should return true()
- The stronger the house, the greater the immigration.
- The Law of Three Little Pigs
- #NGRUM

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `true()`, `false()`.
