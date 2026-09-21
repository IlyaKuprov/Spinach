# kernel/overloads/@ttclass/subsref.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/subsref.m`
- Signature: `answer=subsref(ttrain,reference)`
- Total lines: 126

## Purpose

Dot and bracket property specifications for the tensor train class.

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`, `ttclass_ind2sub()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Syntax

```matlab
answer=subsref(ttrain,reference)
```

## Parameters / inputs

- ttrain -tensor train object
- reference -Matlab subscript reference structure

## Outputs

- answer -requested tensor train property, scalar
- matrix element, or nested subscript result

## Implementation structure

- Dot and bracket property specifications for the tensor train class.
- answer=subsref(ttrain,reference)
- ttrain -tensor train object
- reference -Matlab subscript reference structure
- answer -requested tensor train property, scalar
- matrix element, or nested subscript result
- Methods and properties
- Return the output requested
- Matrix element extraction
- Start with zero
- Convert indices
- Multiply up the tensor train

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `reference()`, `isscalar()`, `sizes()`, `islogical()`, `double()`, `grumble()`, `ttclass_ind2sub()`, `siz()`, `ind()`, `jnd()`, `ivec()`.
