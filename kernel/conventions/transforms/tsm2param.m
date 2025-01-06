% Attempts to convert a traceless symmetric 3x3 interaction matrix into 
% axiality, rhombicity and three Euler angles. The transformation is un-
% stable and should be avoided if at all possible: it is always best to
% just publish the 3x3 matrix as recommended by IUPAC. Syntax:
%
%                      [ax,rh,angles]=tsm2param(M)
%
% Parameters:
%
%    M    -  3x3 matrix or its five independent elements in the 
%            order of [Mxx, Mxy, Mxz, Myy, Myz]
%
% Outputs:
%
%    ax     -  axiality, Mehring order of eigenvalues
%
%    rh     -  rhombicity, Mehring order of eigenvalues
%
%    angles -  Euler angles (one of the eight equivalent 
%              sets), radians
%
% Note: Mehring convention has Z as the largest eigenvalue, and X as 
%       the smallest eigenvalue, this includes signs.
%
% e.suturina@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=tsm2param.m>

function [ax,rh,angles]=tsm2param(M)

% Check consistency
grumble(M);

% Assemble the matrix if needed
if numel(M)==5
    M=[M(1) M(2)  M(3);
       M(2) M(4)  M(5);
       M(3) M(5) -M(1)-M(4)];
end

% Diagonalise the matrix
[V,D]=eig(M); D=diag(D);

% Mehring order eigenvalues
[~,IZ]=max(D); [~,IX]=min(D);
IY=setdiff([1 2 3],[IZ IX]);

% Compute the invariants
ax=2*D(IZ)-(D(IX)+D(IY));
rh=D(IY)-D(IX);

% Compute Euler angles
V=V(:,[IX IY IZ]); 
angles=dcm2euler(V*det(V));

end

% Consistency enforcement
function grumble(M)
if (~isnumeric(M))||(~isreal(M))
    error('M must be a real numeric array.');
end
if (numel(M)~=9)&&(numel(M)~=5)
    error('M must have either nine or five elements.');
end
if (numel(M)==9)&&(~issymmetric(M))
    error('M must be symmetric.');
end
if (numel(M)==9)&&(abs(trace(M))>10*eps)
    error('M must be traceless.');
end
end

% Do that which consists in taking no action, and 
% order will prevail.
%
% The Tao Te Ching

