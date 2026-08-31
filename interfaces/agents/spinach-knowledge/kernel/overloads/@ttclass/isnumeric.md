# kernel/overloads/@ttclass/isnumeric.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/isnumeric.m`
- Signature: `answer=isnumeric(tt)`
- Total lines: 33

## Purpose

Returns TRUE for non-empty tensor train objects. Syntax: answer=isnumeric(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Non-empty tensor trains should return true(); implemented by `if isa(tt,'ttclass')&&(~isempty(tt.cores))`.

### Control flow inferred from the code

- Line 21: conditional branch on `isa(tt,'ttclass')&&(~isempty(tt.cores))`.

### Key state/data transformations

- Lines 22: computes `answer` using `answer=true()`.

## Parameters / inputs

- tt -tensor train object

## Outputs

- answer -logical true for non-empty tensor train objects

## Implementation structure

- Returns TRUE for non-empty tensor train objects. Syntax:
- answer=isnumeric(tt)
- tt -tensor train object
- answer -logical true for non-empty tensor train objects
- Non-empty tensor trains should return true()
- People who think honestly and deeply have a hostile
- attitude towards the public.
- Johann Wolfgang von Goethe

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `true()`, `false()`.
