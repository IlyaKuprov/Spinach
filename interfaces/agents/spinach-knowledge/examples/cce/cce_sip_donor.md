# cce_sip_donor

`[answer,diagnostics]=cce_sip_donor()` calculates the orientation-dependent
Hahn-echo coherence of a phosphorus donor in natural silicon, following
Fig. 5(a) of [Onizhuk and Galli](https://doi.org/10.1002/adts.202100254).
The three columns correspond to [001], [111], and [110] on the same isotope
realisation and 0–1 ms total echo-time grid. The diagnostics record cluster
counts and validity for each direction.

The six-valley Kohn–Luttinger contact term is evaluated in the crystal frame.
The dipolar contribution to the electron–nuclear coupling is added outside
the donor core. The high-field model projects both Zeeman interactions: the retained
hyperfine term is longitudinal, and the secular nuclear dipole retains
the homonuclear flip-flops. Both static Zeeman terms then commute with the
conditional Hamiltonians and are omitted in the rotating frame. The 0.1 T
field scale and the 80 and 10 Å radii come from the authors’ archived
silicon tutorial rather than from the figure caption.

The mixed-product bath permits a native CCE2 calculation in separate edge
batches without changing any cluster Hamiltonian. Each edge appears in one
batch; a monomer shared between batches is divided out of all but the first.
This factorisation uses the unit observed-only ideal Hahn echo and zero
exterior spin means and is not a prescription for sampled or interlaced
bath states.

The published isotope realisation and raw numerical traces are not supplied.
The model is an independent reconstruction, not a replay of those inputs,
and the plotted comparison concerns the paper’s calculated coherence, not
its experimental points.
