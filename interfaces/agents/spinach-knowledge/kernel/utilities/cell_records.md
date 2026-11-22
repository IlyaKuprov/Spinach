# kernel/utilities/cell_records.m

At the input boundary, occupied cells are emitted in column-major order with both indices and the original payload. Empty tables, including an unoccupied singleton, produce a zero-by-three cell array. Only one input column needs an occupancy mask.
