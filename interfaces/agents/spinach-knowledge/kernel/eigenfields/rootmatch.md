# kernel/eigenfields/rootmatch.m

- Signature: `[idx1,idx2,idx3]=rootmatch(field1,field2,field3,...`

## Purpose

Match roots across three magnetic-field root lists while preserving their field order. The routine returns indices into the original input lists.

## Parameters / inputs

- `field1`, `field2`, `field3` — finite real vectors of roots at the first, second, and third triangle vertices, respectively; empty vectors are accepted.
- `edge12`, `edge23`, `edge31` — finite positive real scalar distances between the corresponding pairs of triangle vertices.

## Outputs

- `idx1`, `idx2`, `idx3` — indices of matched roots in `field1`, `field2`, and `field3`, respectively.

## Numerical / algorithmic content

The routine sorts each root list by field, then uses dynamic programming to maximize the number of order-preserving three-way matches. Among solutions with the same number of matches, it minimizes the total field-continuation cost. Each matched triple `(f1,f2,f3)` contributes `((f1-f2)/edge12)^2 + ((f2-f3)/edge23)^2 + ((f3-f1)/edge31)^2`.

If any root list is empty, all outputs are empty. If the three lists have equal lengths, roots are matched directly by sorted order. Output indices refer to the original, unsorted lists and are returned in ascending field order.

[Source documentation](https://spindynamics.org/wiki/index.php?title=rootmatch.m)