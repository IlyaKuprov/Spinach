# kernel/overloads/@ttclass/minus.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/minus.m`
- Signature: `a=minus(a,b)`
- Total lines: 50

## Purpose

Tensor train subtraction operation. Does not perform the actual subtraction but instead concatenates the operands until such time as recompression beco- mes absolutely necessary. Syntax: c=minus(a,b)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- a -a tensor train object
- b -a tensor train object

## Outputs

- c -a tensor train object

## Implementation structure

- Tensor train subtraction operation. Does not perform the actual subtraction
- but instead concatenates the operands until such time as recompression beco-
- mes absolutely necessary. Syntax:
- c=minus(a,b)
- a -a tensor train object
- b -a tensor train object
- c -a tensor train object
- Validate the input
- Write the difference object
- Filter out zero coeff
- Meraki (Greek, n.) -the soul, creativity or love put into something;
- the essence of yourself that is put into your work.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `elseif()`, `all()`, `sizes()`, `unit_like()`.
