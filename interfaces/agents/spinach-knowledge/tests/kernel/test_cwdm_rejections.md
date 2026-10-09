# tests/kernel/test_cwdm_rejections.m

Checks T3: a scalar coupling between substances raises
`Spinach:create:crossSubstanceCoupling`; product operators and states spanning
substances raise `Spinach:which_subst:crossSubstance`, including an explicit
identity factor. Both identifiers and messages are asserted. An out-of-range
local symmetry index is rejected with the local spin-label error message.
