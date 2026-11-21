# Constant-time MAXCRAZED delay sweep

Chen, Branca, and Warren ([Fig. 5](https://doi.org/10.1063/1.4721637)) varied the spacing of the balanced gradient lobes at fixed time before acquisition. The model allocates 60 ms between the end of the mixing pulse and the first acquired point. Three 1.5 ms lobes leave 55.5 ms to divide between the gap separating the balanced lobes and recovery following the returning lobe. This convention differs from assigning a 60 ms *free* inter-lobe gap to the separate Table I benchmark.

For each gap the example returns both complex acquisition dimensions, four peak heights for each encoding relative to X-CRAZED, and the summed mixed/self peak ratio. No parameter is inferred from Fig. 5. The X scan is a timing control: shifting a zero-amplitude lobe should not materially change its peaks. The nominal diffusivities, correlation pitch, solvent magnetisations, and processing assumptions are the same as in the Table I comparison.

At gaps of 0, 15, 30, 45, and 55.5 ms the nominal XZ mixed/self ratio is 0.66, 0.89, 1.34, 2.76, and 10.40. Both X and XY stay near 0.64. Thus the model predicts increasing XZ contrast without a comparable XY change, but overstates the late-time contrast; this is not a quantitative reproduction.
