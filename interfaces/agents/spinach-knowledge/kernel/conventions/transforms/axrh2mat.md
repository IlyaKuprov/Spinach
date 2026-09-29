# kernel/conventions/transforms/axrh2mat.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/axrh2mat.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=axrh2mat.m)

## Contract

axrh2mat constructs a real symmetric 3-by-3 interaction matrix from isotropic value, axiality, rhombicity, and three Euler angles. The interaction's unit is inherited from iso, ax, and rh; the source does not name a particular unit. The angles alp, bet, and gam are in radians and are passed to euler2dcm.

## Principal values and frame

Using Mehring ordering xx<=yy<=zz, the source defines iso=(xx+yy+zz)/3, ax=2*zz-(xx+yy), and rh=yy-xx. It reconstructs the principal values as xx=iso-(ax+3*rh)/6, yy=iso-(ax-3*rh)/6, and zz=iso+ax/3. It then forms R*diag([xx yy zz])*R' with R=euler2dcm(alp,bet,gam) and symmetrises the result as (M+M')/2.

Inputs are real numeric scalars; the source requires rh>=0 and ax>=rh. There are no spin-system fields, grid axes, or spatial dimensions in this transform. The source notes that the inverse transformation is ill-defined.

## Source-supported use

The documented call is M=axrh2mat(iso,ax,rh,alp,bet,gam).
