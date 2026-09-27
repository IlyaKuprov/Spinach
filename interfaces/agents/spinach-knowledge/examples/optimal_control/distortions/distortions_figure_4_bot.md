# examples/optimal_control/distortions/distortions_figure_4_bot.m

- Signature: `distortions_figure_4_bot()`

## Purpose

Figure 4 (bottom) from the paper by Rasulov and Kuprov:

## Physical / mathematical content

- The script models 100 non-interacting 13C spins with chemical shifts spanning −100 to +100 ppm at a magnetic field of 28.18.
- It optimises two-channel, 125-interval RF controls to map initial states Sx, Sy, and Sz to −Sz, Sy, and Sx. The last five intervals are frozen as dead time.
- The control setup specifies `lbfgs`, and the script runs `fmaxnewton(spin_system,@grape_xy,guess)`, labelled LBFGS-GRAPE in the source. The optimisation includes power levels, NS and SNS penalties, and an ensemble of amplifier-saturation distortions.
- A subsequent simulation evaluates infidelity across grids of RF power and amplifier-saturation factors, then plots its logarithm.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Figure 4 (bottom) from the paper by Rasulov and Kuprov:
- Set the magnetic field
- Put 100 non-interacting spins at equal intervals
- within the [-100,+100] ppm chemical shift range
- Select a basis set -IK-2 keeps complete basis on each
- spin in this case, but ignores multi-spin orders
- Run Spinach housekeeping
- Set up spin states
- Get the control operators
- Get the drift Hamiltonian
- Define control parameters
- Last 5 slices are dead time
