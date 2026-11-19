# cce_nv_diamond

Syntax: [answer,diagnostics]=cce_nv_diamond(). The numeric answer is the projected Hahn-echo amplitude; diagnostics is the separate clexpn record.

cce_fixture('nv_four_c13','6f7282b6') supplies the original full physical tensors, four fixed 13C coordinates, source metadata, and literal 41-point grid through 0.5 ms. The example visibly applies observed-Sz projection, sets the mS=-1/+1 transition, mixed-product bath, order two, graph, and numerical policy, and calls clexpn(spin_system,@hahn_echo,parameters,'labframe'). This is a fixed numerical fixture, not natural-abundance diamond.

The independent projected PyCCE complex reference remains in tests/kernel/test_cce_examples.m and now matches the fixed physical receiver without conjugation. Source tensors, total times, and physical-site order are unchanged.
