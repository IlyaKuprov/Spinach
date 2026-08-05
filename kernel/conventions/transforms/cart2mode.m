% Converts Cartesian derivatives of spin Hamiltonian parameters,
% as produced by electronic structure theory packages, into the
% derivatives with respect to dimensionless mode coordinates that
% the bosonic mode specification interface of create.m expects in
% inter.modes.coupling_mod and inter.modes.zeeman_mod fields. The
% dimensionless coordinate of each mode is (a+a')/sqrt(2), and the
% Cartesian displacement that it produces on a degree of freedom
% with mass m is scaled by sqrt(hbar/(m*omega)), where omega is
% the angular frequency of the mode. Syntax:
%
%     mode_derivs=cart2mode(cart_derivs,eigvecs,masses,frqs)
%
% Parameters:
%
%   cart_derivs - first derivatives of an interaction with
%                 respect to Cartesian displacements, in Hz
%                 per Angstrom, an array of dimension
%                 [d1 d2 3N] where N is the number of atoms
%                 and the Cartesian degrees of freedom are
%                 ordered [x1 y1 z1 x2 y2 z2 ...]; or second
%                 derivatives in Hz per Angstrom squared, an
%                 array of dimension [d1 d2 3N 3N]
%
%   eigvecs     - orthonormal mass-weighted normal mode
%                 eigenvector, a [3N 1] column vector for
%                 the first order case; two such vectors
%                 as a [3N 2] array for the second order
%                 case
%
%   masses      - atomic masses in unified atomic mass
%                 units, an [N 1] column vector
%
%   frqs        - mode frequency in Hz, a positive scalar
%                 for the first order case; a [1 2] vector
%                 for the second order case
%
% Outputs:
%
%   mode_derivs - derivatives with respect to dimensionless
%                 mode coordinates, in Hz, a [d1 d2] array
%                 to be placed into the corresponding cell
%                 of inter.modes.coupling_mod (d1=3, d2=3)
%                 or inter.modes.zeeman_mod (d1=1, d2=3)
%
% Note: raw Taylor derivatives are returned; the 1/2 factors of
%       the Taylor expansion are applied by Spinach internally.
%       Derivative data in wavenumbers or meV should be conver-
%       ted into Hz with icm2hz.m or mev2hz.m beforehand. Zero
%       and negative frequency modes are rejected because their
%       zero-point scaling is undefined.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=cart2mode.m>

function mode_derivs=cart2mode(cart_derivs,eigvecs,masses,frqs)

% Check consistency
grumble(cart_derivs,eigvecs,masses,frqs);

% Expand atomic masses over Cartesian degrees of freedom
dof_masses=kron(masses,[1;1;1])*1.66053906892e-27;

% Zero-point displacement scale vectors in Angstrom
scales=1e10*eigvecs.*sqrt((6.62607015e-34/(2*pi))./(dof_masses*(2*pi*frqs)));

% Contract the derivatives with the scale vectors
if ndims(cart_derivs)==3
    mode_derivs=sum(cart_derivs.*reshape(scales,1,1,[]),3);
else
    mode_derivs=sum(sum(cart_derivs.*reshape(scales(:,1),1,1,[]).*...
                                     reshape(scales(:,2),1,1,1,[]),4),3);
end

end

% Consistency enforcement
function grumble(cart_derivs,eigvecs,masses,frqs)
if (~isnumeric(cart_derivs))||(~isreal(cart_derivs))||...
   any(~isfinite(cart_derivs),'all')
    error('cart_derivs must be an array of finite real numbers.');
end
if (ndims(cart_derivs)~=3)&&(ndims(cart_derivs)~=4)
    error('cart_derivs must have three (first order) or four (second order) dimensions.');
end
if (~isnumeric(masses))||(~isreal(masses))||any(~isfinite(masses),'all')||...
   any(masses<=0,'all')||(~iscolumn(masses))
    error('masses must be a column vector of positive real numbers.');
end
if (~isnumeric(eigvecs))||(~isreal(eigvecs))||any(~isfinite(eigvecs),'all')
    error('eigvecs must be an array of finite real numbers.');
end
if size(eigvecs,1)~=3*numel(masses)
    error('the number of rows in eigvecs must be three times the number of masses.');
end
if (~isnumeric(frqs))||(~isreal(frqs))||any(~isfinite(frqs),'all')||any(frqs<=0,'all')
    error('frqs must be an array of positive real numbers.');
end
if ndims(cart_derivs)==3
    if size(cart_derivs,3)~=size(eigvecs,1)
        error('the third dimension of cart_derivs must match the number of rows in eigvecs.');
    end
    if size(eigvecs,2)~=1
        error('first order conversion requires eigvecs with one column.');
    end
    if ~isscalar(frqs)
        error('first order conversion requires a scalar frqs.');
    end
else
    if (size(cart_derivs,3)~=size(eigvecs,1))||(size(cart_derivs,4)~=size(eigvecs,1))
        error('the last two dimensions of cart_derivs must match the number of rows in eigvecs.');
    end
    if size(eigvecs,2)~=2
        error('second order conversion requires eigvecs with two columns.');
    end
    if ~isequal(size(frqs),[1 2])
        error('second order conversion requires frqs to be a 1x2 vector.');
    end
end
if any(abs(sqrt(sum(eigvecs.^2,1))-1)>1e-3)
    error('eigvecs columns must be normalised mass-weighted eigenvectors with unit 2-norm.');
end
end

% The career of a young theoretical physicist consists
% of treating the harmonic oscillator in ever-increasing
% levels of abstraction.
%
% Sidney Coleman


