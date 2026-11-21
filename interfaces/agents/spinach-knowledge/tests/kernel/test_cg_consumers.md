# CG consumer conventions

`tests/kernel/test_cg_consumers.m` tests signed scalar phases before they can
cancel in a product of two coefficients. It checks the reversed spin-half
singlet and the odd exchange phase of unequal coupled angular momenta.

Fresh SLE product operators for interaction ranks one, two, and three are
compared with products of independently evaluated Wigner rotation matrices
at two non-axial orientations. Only source ranks whose complete product fits
inside the retained basis are used, so this is an exact product identity,
not a truncation comparison. The rotation generators also satisfy their
angular-momentum commutator.

PCS checks use an off-diagonal susceptibility tensor: the rank-zero density
must reproduce the point model, and a rank-fourteen expansion of three
non-axial point sources must reproduce their weighted point-model shifts.
