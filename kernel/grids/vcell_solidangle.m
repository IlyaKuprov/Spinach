% Solid angle of a spherical Voronoi cell. Syntax:
%
%            s=vcell_solidangle(P,K,xyz)
%
% Parameters: 
%
%      P    - (3 x m) array with coordinates of the
%              vertices of the Voronoi cell
%
%      K    - (n x 1) cell, each K{j} contains the 
%              indices of the Voronoi cell
%
%      xyz  - optional (3 x n) knot points to guide
%             vcell_solidangle to compute the solid 
%             angle of the "right" cell containing 
%             the node (and not the complement cell)
%
% Outputs:
%
%      S    - the solid angle of the Voronoi cell
%
% brunoluong@yahoo.com
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=vcell_solidangle.m>

function S=vcell_solidangle(P,K,xyz)

% Check consistency
grumble(P);

% Adapt to the input
if nargin<3
    S=cellfun(@(k)one_vcell_solidangle(P(:,k)),K);
else
    S=arrayfun(@(j)one_vcell_solidangle(P(:,K{j}),xyz(:,j)),(1:size(xyz,2)).');
end

end

% Consistency enforcement
function grumble(P)
if any(abs(sum(P.^2,1)-1)>1e-6)
    error('P must contain unit vectors.');
end
end

% Maurice Bowra, the flamboyant warden of Wadham College from 1938 to 1970,
% once argued against the legalisation of homosexuality on the grounds that
% it would take all the fun out of it. Without the risk of being picked up
% by the police, cruising up and down the Cowley Road at one in the morning
% would become rather tedious. He referred to the secret club of powerful
% homosexuals in the British establishment as the "homintern" and prided
% himself on being a high-ranking officer. He liked the fact that there was
% something exotic and clandestine about his sexuality and dreaded the risk
% of embourgeoisement if the law was changed.
%
% Toby Young

