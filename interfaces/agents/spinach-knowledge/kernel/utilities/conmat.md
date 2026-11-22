# kernel/utilities/conmat.m

The call `pairs=conmat(xyz,images,r0)` returns ordered records
`{from,to,images(to,:)-xyz(from,:)}` strictly inside the cutoff. Both
coordinate arrays contain the same atoms in the same order, in one reference
cell and one translated image cell. Use the same array twice for an isolated
molecule. Self images are excluded; zero cutoff returns no records.

Only occupied cutoff-sized cells and their neighbours are examined. Records
are sorted by destination and then source, without a square coordinate mask.

## References

- Source: [kernel/utilities/conmat.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/conmat.m)
- Spinach Wiki: [conmat.m](https://spindynamics.org/wiki/index.php?title=conmat.m)
