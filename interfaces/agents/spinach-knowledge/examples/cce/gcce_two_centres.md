# Proton-pair contribution to the double-quantum triplet echo

`[loss,axes_hz]=gcce_two_centres()` computes Fig. 3(b) of [Chen et al.](https://doi.org/10.1063/5.0303171): the loss of double-quantum Hahn-echo coherence from one antiparallel proton pair at total evolution time 40 µs. It is not the random-bath average of Fig. 4. The result is a surface with rows indexed by `axes_hz.flip` and columns by `axes_hz.detuning`, both in Hz.

The two electron triplets with magnetic quantum numbers ±1 select proton Hamiltonians with opposite hyperfine detunings. The proton chemical-shift difference supplies the common 5 kHz energy cost. The transverse dipolar term mixes the two antiparallel proton states. The echo disappears from the loss surface when either this mixing or the branch-dependent energy difference vanishes. Order two retains both protons together, so the cluster expansion is exact for the Hamiltonian of this panel. The signs of the two sweep axes leave the magnitude invariant.
