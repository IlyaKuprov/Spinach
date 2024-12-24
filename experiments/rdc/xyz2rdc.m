% Converts Cartesian coordinates of a pair of nuclei and an 
% order matrix into residual dipolar coupling. Syntax:
%
%         rdc=xyz2rdc(spin_a,spin_b,xyz_a,xyz_b,chi)
%
% Parameters:
%
%    spin_a, spin_b - character strings indicating spin
%                     type, for example '13C'
%
%    xyz_a, xyz_b   - three-element vectors specifying
%                     Cartesian coordinates of the two
%                     spins in Angstroms
%
%    order_spec     - {S,'saupe'} uses Saupe order mat-
%                     rix, S is a traceless symmetric
%                     3x3 matrix, dimensionless
%
% Outputs:
%
%    rdc            - residual dipolar coupling in the
%                     heteronuclear case, Hz
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=xyz2rdc.m>

function rdc=xyz2rdc(spin_a,spin_b,xyz_a,xyz_b,order_spec)

% Check consistency
grumble(spin_a,spin_b,xyz_a,xyz_b,order_spec);

% Choose the path
switch order_spec{2}

    % Saupe order matrix
    case 'saupe'

        % Get the dipole coupling tensor in rad/s
        [~,~,~,~,D]=xyz2dd(xyz_a,xyz_b,spin_a,spin_b);

        % Get weak heteronuclear RDC in Hz
        rdc=(2/3)*trace(order_spec{1}*D)/(2*pi);
    
    % For future options
    otherwise

        % Complain and bomb out
        error('unknown order specification type');

end

end

% Consistency enforcement
function grumble(spin_a,spin_b,xyz_a,xyz_b,order_spec)
if (~iscell(order_spec))||(~ismember(numel(order_spec),[2 4]))
    error('order_spec must be a cell array with 2 or 4 elements.');
end
if (~isnumeric(order_spec{1}))||(~ismatrix(order_spec{1}))||...
   (any(size(order_spec{1})~=[3 3]))||(~isreal(order_spec{1}))
    error('order_spec{1} must be a real 3x3 matrix.');
end
if trace(order_spec{1})>10*norm(order_spec{1},2)*eps()
    error('order_spec{1} must be traceless.');
end
if norm(order_spec{1}-order_spec{1}.',1)>10*norm(order_spec{1},2)*eps()
    error('order_spec{1} must be symmetric.');
end
if (~ischar(spin_a))||(~ischar(spin_b))
    error('spin_a and spin_b must be character strings.');
end
if strcmp(spin_a,spin_b)
    error('spins A and B must be different.');
end
if (~isreal(xyz_a))||(~isnumeric(xyz_a))||(numel(xyz_a)~=3)||...
   (~isreal(xyz_b))||(~isnumeric(xyz_b))||(numel(xyz_b)~=3)
    error('xyz_a and xyz_b must be real vectors with three elements each.');
end
end

% What cannot be said will be wept.
%
% Sappho

