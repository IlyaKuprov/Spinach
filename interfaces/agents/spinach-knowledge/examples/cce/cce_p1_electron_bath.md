# cce_p1_electron_bath

`[answer,diagnostics]=cce_p1_electron_bath()` calculates the relative
P1-pair energies and active electron flip-flop channels of Park et al.,
[Fig. 5a,c,d](https://doi.org/10.1038/s41534-022-00605-4), at 500 G along
[111]. A P1 centre comprises an electron and a co-located ¹⁴N nucleus.
The pair diagrams use the two opposite electron-spin configurations and
the nine nuclear product states labelled in Fig. 5b. The flip-flop
operator must connect the two states, and their diagonal high-field
Zeeman/hyperfine energies must match. The resulting channel counts are
9 with hyperfine off, 3 with both Jahn–Teller axes aligned, and 1 with
one axis tilted. The plotted energies are relative; the schematic does
not resolve nuclear quadrupolar splittings.

The separate NV Hahn echo is a projected partition CCE-2 calculation on
one explicit diamond-site pair, with both P1 centres retained as
electron–nucleus groups. Unlike the level diagram, the echo retains
full on-site anisotropic hyperfine (when enabled) and ¹⁴N quadrupole
tensors from Supplementary Note 2. It uses secular inter-centre dipolar
couplings, the NV `m_s=0,+1` transition, and a maximally mixed bath.
This single pair illustrates how the three models affect the echo; it
does not reproduce the bath ensemble or its decay times in Fig. 4.

The columns of `answer` give complex coherence in the model order above.
`diagnostics` contains the total echo times in seconds, validity masks,
cluster counts, nuclear projections in Fig. 5b order, the nine energy
levels in MHz for each opposite electron configuration, and the
corresponding 9-by-9 resonance matrices.
