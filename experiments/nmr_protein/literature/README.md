# Protein NMR literature

This folder contains literature PDFs used to audit the phase-sensitive
`experiments/nmr_protein` pulse sequence files. Where a clean public PDF could
not be downloaded, `DOI_FALLBACKS.md` gives the DOI for manual retrieval.

## Downloaded PDFs

- `kay_ikura_tschudin_bax_jmr_1990_triple_resonance.pdf`
  - Kay, L.E.; Ikura, M.; Tschudin, R.; Bax, A. Three-dimensional
    triple-resonance NMR spectroscopy of isotopically enriched proteins.
    Journal of Magnetic Resonance 89, 496-514 (1990).
  - DOI: `10.1016/0022-2364(90)90333-5`
  - Used for `hnca.m`, `hnco.m`, and the general triple-resonance framework.

- `bax_ikura_jbnmr_1991_hncoca.pdf`
  - Bax, A.; Ikura, M. An efficient 3D NMR technique for correlating the
    proton and 15N backbone amide resonances with the alpha-carbon of the
    preceding residue in uniformly 15N/13C enriched proteins.
    Journal of Biomolecular NMR 1, 99-104 (1991).
  - DOI: `10.1007/BF01874573`
  - Used for `hncoca.m`.

- `bax_clore_driscoll_gronenborn_ikura_kay_jmr_1990_hcch_cosy.pdf`
  - Bax, A.; Clore, G.M.; Driscoll, P.C.; Gronenborn, A.M.; Ikura, M.;
    Kay, L.E. Practical aspects of proton-carbon-carbon-proton
    three-dimensional correlation spectroscopy of 13C-labeled proteins.
    Journal of Magnetic Resonance 87, 620-627 (1990).
  - Used for `hcch_cosy.m`.

- `bax_clore_gronenborn_jmr_1990_hcch_tocsy.pdf`
  - Bax, A.; Clore, G.M.; Gronenborn, A.M. 1H-1H correlation via isotropic
    mixing of 13C magnetization, a new three-dimensional approach for
    assigning 1H and 13C spectra of 13C-enriched proteins.
    Journal of Magnetic Resonance 88, 425-431 (1990).
  - DOI: `10.1016/0022-2364(90)90202-K`
  - Used for `hcch_tocsy.m`.

- `edwards_savostyanov_welderufael_lee_kuprov_jmr_2014_spinach_stitch.pdf`
  - Edwards, L.J.; Savostyanov, D.V.; Welderufael, Z.T.; Lee, D.;
    Kuprov, I. Quantum mechanical NMR simulation algorithm for
    protein-size spin systems. Journal of Magnetic Resonance 243, 107-113
    (2014).
  - DOI: `10.1016/j.jmr.2014.04.002`
  - Used for the Spinach bidirectional propagation and stitching method.

## Sequence coverage

- `hnca.m`: downloaded Kay/Ikura/Tschudin/Bax 1990 PDF.
- `hnco.m`: downloaded Kay/Ikura/Tschudin/Bax 1990 PDF.
- `hncoca.m`: downloaded Bax/Ikura 1991 PDF.
- `hncaco.m`: DOI fallback for Clubb/Thanabal/Wagner 1992.
- `hcanh.m`: DOI fallbacks for Montelione/Wagner 1990 and
  Feng/Rios/Montelione 1996.
- `hcch_cosy.m`: downloaded Bax/Clore/Driscoll/Gronenborn/Ikura/Kay 1990
  PDF, with DOI fallback for the short JACS HCCH-COSY concept paper.
- `hcch_tocsy.m`: downloaded Bax/Clore/Gronenborn 1990 PDF, with DOI
  fallbacks for gradient-enhanced and 4D HCCH-TOCSY variants.

