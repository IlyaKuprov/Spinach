# test_clexpn_descriptors — canonical operator and state conformance
Keywords: test_clexpn_descriptors cce_operator cce_state pauli physical IDs populations ordering Hilbert Liouville

Registered as kernel/clexpn_descriptors. Compares all six supported spin labels against stock Hilbert operators, checks arbitrary-spin normalisation, isotope-group sums and products in permuted observed-first order, and compares a product state with stock state. Explicit nonuniform populations preserve trace and tensor-product order. Liouville commutators and column-major states are checked independently. Empty-bath observed operators remain valid; missing populations, discarded physical IDs and unsupported labels are rejected.
