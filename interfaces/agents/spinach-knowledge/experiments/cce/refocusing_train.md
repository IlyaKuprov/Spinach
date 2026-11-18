# refocusing_train

Syntax: [answer,aux]=refocusing_train(spin_system,parameters,H,R,K).

Ideal generalised Hilbert-space refocusing with explicit preparation/detection
following free_evolution. Supply either a refocusing cell row of descriptors
(spins, integrated Hermitian operator or axis/angle, and time or fraction), or
train with name, count, axis, and spins. XY8 has literal xyxyyxyx phases and
uniform half-interval placements; its count is a positive multiple of eight.
Uhrig fractions are sin(pi*k/(2*count+2))^2. These are not renamed CPMG sequences.

The drift acts between pulses, zero elapsed time is the initial state, and an
absolute pulse beyond a sample is omitted. Arbitrary retained-spin operator
supports preserve their declared tensor order. Cartesian exterior rotations
update static product-state mean fields; correlated exterior pulses are refused.
Auxiliary outputs retain density, arithmetic bounds, and compiled segments.
Projected and master-equation trains are not registered by this experiment;
use the already-supported hahn_echo for a single admissible master echo.

Caller-supplied parameters.sequence is rejected, including empty arrays. Use finite_drive, refocusing_train, or the appropriate named echo experiment; internal events are constructed from native physical declarations, never read from that removed field.

Elapsed times and the named-train or ideal-refocusing descriptor fields are checked locally before acquisition.
