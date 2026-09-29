# kernel/utilities/binpack.m

## Purpose

`binpack.m` implements a simple one-dimensional bin packing algorithm. It collects a supplied list of numbers into sublists whose sums are smaller than or equal to a specified bin size. The header comment states that the algorithm is not optimal, but that it does the job.

## Behaviour

The function is called as `bins=binpack(box_sizes,bin_size)`.

1. Input consistency is enforced by an internal `grumble` subfunction, which errors if `box_sizes` is not a numeric, real, finite row vector of positive integers, or if `bin_size` is not a single positive real integer.
2. Boxes are numbered by their indices, `box_index=(1:numel(box_sizes))'`.
3. Boxes larger than `bin_size` (`box_sizes>bin_size`) are each placed into their own bin as singleton index vectors and removed from the packing list.
4. Remaining boxes are packed iteratively: in each pass, `cumsum(box_sizes)<=bin_size` selects the leading run of boxes whose cumulative sum fits within the bin; those indices are appended as a new bin, and the packed boxes are removed from the lists. The loop repeats until no boxes remain.

Because the packing loop consumes boxes in list order via a cumulative-sum cutoff, each bin contains a contiguous prefix of the remaining box list; the header notes the result is not an optimal packing.

## Inputs and outputs

Inputs:

- `box_sizes` — a row vector of box sizes; must be a numeric, real, finite row vector of positive integers.
- `bin_size` — an integer specifying the bin size; must be a single positive real integer.

Output:

- `bins` — a cell array of index vectors specifying the boxes allocated into each bin. Oversized boxes each occupy their own single-element bin.

## References

- Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/binpack.m>
- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=binpack.m>
