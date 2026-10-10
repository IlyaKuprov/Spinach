# kernel/eigenfields/rootmatch.m

- Direct source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/eigenfields/rootmatch.m
- Spinach Wiki: https://spindynamics.org/wiki/index.php?title=rootmatch.m

## Signature

`[idx1,idx2,idx3]=rootmatch(field1,field2,field3,edge12,edge23,edge31)`

## Purpose

Match roots across three magnetic-field root lists while preserving their order. The returned values are indices into the original input lists, not the sorted copies.

## Inputs

- `field1`, `field2`, `field3`: finite real numeric vectors of roots at the three triangle vertices. Row and column vectors are accepted; empty vectors are allowed.
- `edge12`, `edge23`, `edge31`: finite positive real numeric scalars giving the respective vertex separations. Each edge must use the same distance unit as the field coordinates it relates.

## Outputs

- `idx1`, `idx2`, `idx3`: row vectors of matched indices in the original input lists, returned in ascending field order. If any root list is empty, all three outputs are empty.

## Algorithm and edge cases

Each list is sorted by field. When all three lists have equal lengths, the routine pairs roots by sorted order and returns immediately. Otherwise, a three-dimensional dynamic programme considers skipping the next root from any one list or matching the next root from all three. It first maximises the number of order-preserving triples, then minimises total continuation cost among solutions with that match count. A triple `(f1,f2,f3)` contributes

`((f1-f2)/edge12)^2 + ((f2-f3)/edge23)^2 + ((f3-f1)/edge31)^2`.

Traceback restores indices to the original lists and the forward field order. The source contains no worked numeric example.
