# kernel/overloads/@ttclass/subsref.m

- Signature: `answer=subsref(ttrain,reference)`

## Purpose

Dot and bracket property specifications for the tensor train class.

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

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
