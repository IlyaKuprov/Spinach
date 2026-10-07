# kernel/utilities/zte_warr.m

`duration=zte_warr(spin_system,L,rho,projector)` reports and returns a cheap
leakage-based time estimate in seconds. `L` is the fixed square Liouvillian,
`rho` a matching column vector, and `projector` the coordinate embedding
returned by `zte` (including its scalar `1` unchanged-space convention).

For retained coordinates S and discarded coordinates D, let
`d=norm(rho(D))`, `r=norm(rho(S))`, and `b=norm(L(D,S),'fro')`. The
estimated time is `(spin_system.tols.zte_warr-d)/(b*r)`. The Frobenius norm
bounds the block's spectral norm without eigensolvers. Initial discard at
or above the tolerance gives zero; zero estimated leakage otherwise gives
`Inf`. No reduction and zero input states give `Inf`.

This is a linear leakage estimate, not a general guarantee: subsequent
amplification, whether the generator is normal or non-normal, is neglected. An `Inf` estimate
with nonzero initial discard does not certify its later evolution. Only the
supplied vector and fixed generator are considered; other preparations,
frequency-domain `slowpass` spectra, propagation error, and roundoff are not
certified. The positive finite tolerance defaults to `1e-6` and is independent
of the pruning tolerance `zte_tol`.
