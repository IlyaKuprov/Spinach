# examples/cce/gcce_two_centres.m

Signature: [answer,diagnostics]=gcce_two_centres()

Two observed electrons coupled to two bath nuclei undergo free_evolution through clexpn with generalized CCE order 2 and a zeeman-hilb retained basis. Standard sys/inter/bas/parameters blocks specify the initial equal superposition of the two-electron all-up and all-down states using transition.spins=[1 2] and transition.ms=[0.5 -0.5;0.5 -0.5], not manually assembled vectors.

answer has one row per time and column-vectorized 4-by-4 observed density matrices (31-by-16). The example plots their four populations. Independent full-bath Kronecker/expm propagation, partial traces, negativity, and every associated physical array are retained in test_cce_examples. The illustrative finite bath is not a material fit or a thermodynamic decay-time calculation.
