# examples/quantum_tech/transmon_duffing_ladder.m

Source: [examples/quantum_tech/transmon_duffing_ladder.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/transmon_duffing_ladder.m)

- Signature: `transmon_duffing_ladder()`

## Model and calculation

The isotope label `T5` represents a single five-level transmon-mode truncation; it is not an electron spin, defect isotope, or EPR system. The mode frequency is `5.0e9` (5 GHz). In the lab-frame Hamiltonian the code adds a Duffing term, using the `CCAA` operator, and sweeps the anharmonicity over 80 values from −400 to −50 MHz (the code stores this range as angular-frequency values, `2*pi*linspace(-400e6,-50e6,80)`). No time-dependent drive or dissipative interaction is configured: each point is a static Hamiltonian calculation.

For each value, the five eigenenergies are sorted and adjacent differences give the 0–1, 1–2, 2–3, and 3–4 transition frequencies. The plot shows those four transitions in GHz against the positive quantity −α/(2π) in MHz, so increasing horizontal coordinate means increasing magnitude of the negative anharmonicity. It is an energy-ladder comparison across the specified parameter sweep, not a measured spectrum or device-performance result.

The source cites Koch et al., *Physical Review A* **76**, 042319 (2007) ([DOI](https://doi.org/10.1103/PhysRevA.76.042319)).
