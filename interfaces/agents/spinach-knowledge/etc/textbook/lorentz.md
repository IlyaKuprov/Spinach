# etc/textbook/lorentz.m

- MATLAB implementation: [etc/textbook/lorentz.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/textbook/lorentz.m)

- Signature: `[J,K,Kil]=lorentz(L)`

## Purpose

Constructs the direct-sum (L,0) plus (0,L) matrix representation of the Lorentz group with inversion. Use it when explicit rotation and boost generators are needed; request only two outputs if the Killing form is not needed.

## Input

- L — required real numeric scalar representation rank; it must be an integer or half-integer with L >= 1/2. There are no defaults. The source sets D = 2L+1 and obtains the spin matrices from pauli(D).

## Construction and outputs

For each spin matrix s.x, s.y, s.z, J has s in both diagonal blocks, while K has +i s in the first block and -i s in the second. Each component is returned as a full 2D-by-2D matrix.

- J — structure with rotation generators J.x, J.y, J.z.
- K — structure with boost generators K.x, K.y, K.z.
- Kil — optional 6-by-6 Killing form. The code computes it only when a third output is requested (nargout > 2), by forming adjoint-representation matrices for the six generators and taking pairwise trace products. This is the expensive part; [J,K] avoids it.

## Source

[Spinach Wiki: lorentz.m](https://spindynamics.org/wiki/index.php?title=lorentz.m).
