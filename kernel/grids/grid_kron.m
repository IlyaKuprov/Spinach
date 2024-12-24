% Spherical grid direct product. Tiles one grid using the rotations of
% the other. Grids should be supplied using Euler angles in three col-
% umns [alphas betas gammas] in radians. Syntax:
%
%    [angles,weights]=grid_kron(angles1,weights1,angles2,weights2)
%
% Parameters:
%
%    angles1 - Euler angles (ZYZ active) of the first grid,
%              as [alpha beta gamma], radians
%
%   weights1 - weights of the first grid
%
%    angles2 - Euler angles (ZYZ active) of the second grid,
%              as [alpha beta gamma], radians
%
%   weights1 - weights of the second grid
%
% Outputs:
%
%     angles - Euler angles (ZYZ active) of the product grid,
%              as [alpha beta gamma], rad
%
%    weights - weights of the product grid
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=grid_kron.m>

function [angles,weights]=grid_kron(angles1,weights1,angles2,weights2)

% Check consistency
grumble(angles1,weights1,angles2,weights2);

% Convert both grids to quaternions
quats1=angle2quat(angles1(:,1),angles1(:,2),angles1(:,3),'ZYZ');
quats2=angle2quat(angles2(:,1),angles2(:,2),angles2(:,3),'ZYZ');

% Build a table of quaternion products
quats=[kron(quats1,ones(size(quats2,1),1)) kron(ones(size(quats1,1)),quats2)];

% Multiply up quaternions
quats=[quats(:,1).*quats(:,5)-quats(:,2).*quats(:,6)-quats(:,3).*quats(:,7)-quats(:,4).*quats(:,8),...
       quats(:,1).*quats(:,6)+quats(:,2).*quats(:,5)+quats(:,3).*quats(:,8)-quats(:,4).*quats(:,7),...
       quats(:,1).*quats(:,7)-quats(:,2).*quats(:,8)+quats(:,3).*quats(:,5)+quats(:,4).*quats(:,6),...
       quats(:,1).*quats(:,8)+quats(:,2).*quats(:,7)-quats(:,3).*quats(:,6)+quats(:,4).*quats(:,5)];
   
% Convert quaternions into angles
[alphas,betas,gammas]=quat2angle(quats,'ZYZ'); angles=real([alphas betas gammas]);

% Tile the weights
weights=kron(weights1,weights2);

end

% Consistency enforcement
function grumble(angles1,weights1,angles2,weights2)
if (~isnumeric(angles1))||(~isreal(angles1))||(size(angles1,2)~=3)
    error('angles1 must be a real matrix with three columns.');
end
if (~isnumeric(angles2))||(~isreal(angles2))||(size(angles2,2)~=3)
    error('angles2 must be a real matrix with three columns.');
end
if (~isnumeric(weights1))||(~isreal(weights1))||...
   any(~isfinite(weights1))||(size(weights1,2)~=1)
    error('weights1 must be a column vector of real numbers.');
end
if (~isnumeric(weights2))||(~isreal(weights2))||...
   any(~isfinite(weights2))||(size(weights2,2)~=1)
    error('weights1 must be a column vector of real numbers.');
end
if size(angles1,1)~=size(weights1,1)
    error('the number of rows in angles1 and weights1 must be the same.');
end
if size(angles2,1)~=size(weights2,1)
    error('the number of rows in angles2 and weights2 must be the same.');
end
end

% Dostoevsky's lack of taste, his monotonous dealings with persons
% suffering with pre-Freudian complexes, the way he has of wallow-
% ing in the tragic misadventures of human dignity -- all this is
% difficult to admire.
%
% Vladimir Nabokov

