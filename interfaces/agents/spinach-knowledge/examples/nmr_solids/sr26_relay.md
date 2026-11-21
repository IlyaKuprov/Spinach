# examples/nmr_solids/sr26_relay.m

Compares direct and relay-only dipolar networks under SR26 and S3. The example checks the symmetry schedule before compiling its segments and retains total coherence orders ±2 between identical arms using `relay_echo`. The schedules come from `homonuclear_recoupling/SR26_CC.in` and `homonuclear_recoupling/S3SR26_SiSi_varydist.in` at dnp-grenoble revision `fa4d4d9ec2f49aed9a01bac7343f9e166121ad34`; the modified silicon network is not a reproduction of either powder benchmark.
