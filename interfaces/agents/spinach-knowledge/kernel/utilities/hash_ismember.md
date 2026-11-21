# kernel/utilities/hash_ismember.m

`[present,indices]=hash_ismember(A,hash_b)` matches rows of a real, finite double query matrix against prehashed reference rows. `A` may be sparse or full. `hash_b` must be an N-by-32 character array, with row k equal to `md5_hash(full(B(k,:)))` for a real, finite double reference matrix B having the same column count as A.

Only the query rows are hashed, after conversion to full storage. Hash membership returns logical and double column vectors, respectively, with the first matching reference index or zero for an absent row. Empty query and reference arrays are supported, as are zero-column rows. MD5 collisions are not checked. Reference column counts cannot be recovered or validated from the hashes; the caller supplies hashes for the appropriate descriptor columns.

`basis` stores the reusable reference arrays in `spin_system.bas.row_hash`, one cell per substance. This row cache is not the aggregate `bas.basis_hash` used by other caches.
