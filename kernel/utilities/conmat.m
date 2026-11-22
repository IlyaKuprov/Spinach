% Molecular connectivity matrix calculator with N*log(N)
% asymptotic complexity scaling with respect to the num-
% ber or atoms. Syntax:
%
%                conmatrix=conmat(xyz,r0)
%
% Parameters:
%
%       xyz  - an array with N rows and three columns,
%              giving the Cartesian coordinates of
%              each particle
%
%        r0  - the distance below which the particles
%              are to be considered "connected"
%
% Output:
%
%  conmatrix - a sparse logical matrix containing 1
%              at the positions corresponding to the
%              connected particles
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=conmat.m>

function conmatrix=conmat(xyz,r0)

% Check consistency
grumble(xyz,r0);

% Find neighbours and apply the open Cartesian cutoff planes
[from,to]=nearby_points(xyz,xyz,r0);
keep=(from~=to)&all(abs(xyz(from,:)-xyz(to,:))<r0,2);
conmatrix=sparse(from(keep),to(keep),true,size(xyz,1),size(xyz,1));

end

% Consistency enforcement
function grumble(xyz,r0)
if (~isnumeric(xyz))||(size(xyz,2)~=3)||(~isreal(xyz))
    error('xyz parameter should be a real matrix with three columns.');
end
if (~isnumeric(r0))||(~isreal(r0))||(r0<=0)
    error('r0 must be a positive real number.');
end
end

% It's always good to be underestimated.
%
% Donald Trump


