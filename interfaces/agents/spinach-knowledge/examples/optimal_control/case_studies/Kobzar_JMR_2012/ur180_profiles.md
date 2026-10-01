# examples/optimal_control/case_studies/Kobzar_JMR_2012/ur180_profiles.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/case_studies/Kobzar_JMR_2012/ur180_profiles.m) · [Kobzar et al., J. Magn. Reson. 225, 142–160 (2012)](https://doi.org/10.1016/j.jmr.2012.09.013)

Simulates the author's deposited `UR180_730u` y-axis universal π-rotation waveform from the [BURBOP archive](https://www.ioc.kit.edu/luy/311.php), not the different 2 ms pulse of Skinner et al., JMR 216 (2012). Its original header specifies 730 µs total duration, 0.5 µs steps, 20 kHz total offset bandwidth, five RF scalings spanning ±40%, and 10 kHz nominal maximum RF. The repository data file retains its 1,460 numerical x/y/duration rows.

`[profiles,fig]=ur180_profiles(linspace(-10e3,10e3,100),linspace(0.6,1.4,5))` independently propagates x, y and z spin states through that waveform on a 100-by-5 offset–RF grid. It returns the three-state rotation score and the single-spin SU(2) propagator-overlap score `sqrt((1+3*fidelity)/4)`; the plotted map and profiles use the latter. Both vectors are required inputs; this call uses the source grid. This is a simulation of deposited numerical data, not a reoptimisation or reproduction of the other 2012 paper's Fig. 4.
