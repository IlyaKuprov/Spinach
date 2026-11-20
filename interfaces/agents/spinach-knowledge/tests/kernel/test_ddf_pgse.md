# Finite-lobe diffusion control

A balanced rectangular gradient pair, with feedback disabled, measures diffusion independently of dipolar echo formation. The production sequence runs both the encoded acquisition and its zero-gradient reference. For lobe duration δ and onset separation Δ, the complex ratio is real and equals exp[−D(γGδ)²(Δ−δ/3)] in unrestricted liquid. An instantaneous-kick implementation would miss the δ/3 correction. Detecting only the middle third of the reflecting slab suppresses wall effects; rotating both the slab and receiver tests the same physical attenuation along X and Z.

Reference: [Stejskal and Tanner (1965)](https://doi.org/10.1063/1.1695690).
