%KEHL_DIAG_OPS Generates diagonal-frame Kehl relaxation operators.
%
%   [SX_D,SY_D,SZ_D,IX_D,IY_D,IZ_D]=KEHL_DIAG_OPS(SX,SY,SZ,IX,IY,IZ,
%   N_ENDOR) applies the Kehl diagonal-frame transformation to the
%   explicit electron and nuclear spin operators supplied by the caller.
%
%   Inputs:
%
%      SX,SY,SZ - electron spin operators.
%      IX,IY,IZ - cell arrays of nuclear spin operators.
%      N_ENDOR  - number of ENDOR nuclei.
%
%   Outputs:
%
%      SX_D,SY_D,SZ_D - transformed electron spin operators.
%      IX_D,IY_D,IZ_D - transformed nuclear spin operators.
%
%   February 2024 A. Kehl (akehl@gwdg.de)
%   Spinach architecture migration May 2026 Talos

function [Sx_D,Sy_D,Sz_D,Ix_D,Iy_D,Iz_D]=kehl_diag_ops(Sx,Sy,Sz,Ix,Iy,Iz,n_endor)

% Check consistency
grumble(Sx,Sy,Sz,Ix,Iy,Iz,n_endor);

% Preserve the identity transformation matrix
D=eye(size(Sz));

% Transform electron operators explicitly
Sx_D=D\Sx*D;
Sy_D=D\Sy*D;
Sz_D=D\Sz*D;

% Preallocate transformed nuclear operators
Ix_D=cell(1,n_endor);
Iy_D=cell(1,n_endor);
Iz_D=cell(1,n_endor);

% Transform nuclear operators explicitly
for n=1:numel(Ix)
    Ix_D{n}=D/Ix{n}*D;
    Iy_D{n}=D/Iy{n}*D;
    Iz_D{n}=D/Iz{n}*D;
end

end

% Consistency enforcement
function grumble(Sx,Sy,Sz,Ix,Iy,Iz,n_endor)
if (~isnumeric(Sx))||(~isnumeric(Sy))||(~isnumeric(Sz))
    error('electron spin operators must be numeric.');
end
if (~iscell(Ix))||(~iscell(Iy))||(~iscell(Iz))
    error('nuclear spin operators must be cell arrays.');
end
if (~isnumeric(n_endor))||(~isscalar(n_endor))||...
   (n_endor<0)||mod(n_endor,1)~=0
    error('n_endor must be a non-negative integer.');
end
end
