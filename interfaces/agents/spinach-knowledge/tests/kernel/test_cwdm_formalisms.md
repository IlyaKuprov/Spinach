# CWDM Zeeman formalism capabilities

`test_cwdm_formalisms` checks local dimensions and cumulative offsets for unequal Zeeman blocks and a spin-free substance. It requires empty descriptor cells for all three Zeeman formalisms. Wavefunction thermal equilibrium, concentration-weighted units, both thermalisation methods, undefined ket concentrations, first-order chemistry, and mass-action chemistry must raise their exact documented identifiers and messages. The chemistry checks exercise both basis and kinetics entry points.

Unweighted `coil_state` storage is compared exactly with the direct sum of local pure kets. The weighted `state` wrapper must reject a segmented wavefunction request with its complete concentration-weighting message.

Missing basis metadata and a missing formalism field must retain the explicit unit-state input-validation message rather than a MATLAB missing-field exception.
