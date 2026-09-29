# examples/fundamentals/convention_tests/tensors.m

- MATLAB implementation: [examples/fundamentals/convention_tests/tensors.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/convention_tests/tensors.m)

- Signature: `tensors()`

## Question tested

For a spin-15/2 system (multiplicity 16), do Stevens coefficient vectors and their `stev2sph`-converted IST coefficient vectors build the same operator matrix when expanded in their respective operator sets?

## Construction and comparison

The test covers ranks `k=1,…,6`. For each rank, it draws `2k+1` real coefficients as `rand(2*k+1,1)-0.5` and uses them with `stevens(mult,k,q-k-1)`, so the index `q=1,…,2k+1` runs over component labels `-k,…,k`. It sums those terms over all ranks into `lin_comb_a`. It then transforms each rank's coefficient vector with `stev2sph(k,r{k})`, obtains the rank's operators from `irr_sph_ten(mult,k)`, and sums the converted coefficients against those operators into `lin_comb_b`.

The criterion is `norm(lin_comb_b-lin_comb_a,1)<1e-6`. If it fails, the source displays both full matrices and raises an error; on the other branch it displays a success message. This is a matrix-level check of the coefficient conversion together with reconstruction across the stated ranks; it is not a claim about any run's outcome.
