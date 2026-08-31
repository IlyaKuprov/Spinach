# kernel/overloads/@ttclass/minus.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/minus.m`
- Signature: `a=minus(a,b)`
- Total lines: 50

## Purpose

Tensor train subtraction operation. Does not perform the actual subtraction but instead concatenates the operands until such time as recompression beco- mes absolutely necessary. Syntax: c=minus(a,b)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Validate the input; implemented by `if (~isa(a,'ttclass'))||(~isa(b,'ttclass'))`.
- Lines 31-32: Write the difference object; implemented by `a.coeff=[a.coeff -b.coeff]`.
- Lines 36-37: Filter out zero coeff; implemented by `pos=find(a.coeff)`.

### Control flow inferred from the code

- Line 25: conditional branch on `(~isa(a,'ttclass'))||(~isa(b,'ttclass'))`.
- Line 38: conditional branch on `~isempty(pos)`.

### Key state/data transformations

- Lines 32: computes `a.coeff` using `a.coeff=[a.coeff -b.coeff]`.
- Lines 33: computes `a.cores` using `a.cores=[a.cores b.cores]`.
- Lines 34: computes `a.tolerance` using `a.tolerance=[a.tolerance b.tolerance]`.
- Lines 37: computes `pos` using `pos=find(a.coeff)`.
- Lines 43: computes `a` using `a=0*unit_like(a)`.

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
