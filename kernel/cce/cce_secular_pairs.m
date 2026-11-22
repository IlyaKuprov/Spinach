% Secular CCE-2 Hahn echo for an unpolarised bath of identical spin-1/2
% nuclei with longitudinal dipolar hyperfine couplings. Syntax:
%
%  echo=cce_secular_pairs(xyz,field_axis,hyperfine,dipolar,times,abs_tol)
%
% Parameters:
%
%     xyz        - N-by-3 array of centre-relative nuclear sites, in metres
%
%     field_axis - 3-by-1 magnetic field unit vector
%
%     hyperfine  - point-dipole hyperfine prefactor, in rad s^-1 m^3
%
%     dipolar    - nuclear point-dipole prefactor, in rad s^-1 m^3
%
%     times      - column of total Hahn echo times, in seconds
%
%     abs_tol    - upper bound on echo loss omitted by pair screening
%
% Outputs:
%
%     echo       - column of dimensionless normalised echo amplitudes
%
% ilya.kuprov@weizmann.ac.il

function echo=cce_secular_pairs(xyz,field_axis,hyperfine,dipolar,times,abs_tol)

% Check consistency
grumble(xyz,field_axis,hyperfine,dipolar,times,abs_tol);

% Longitudinal hyperfine; central and common nuclear Zeeman phases refocus
distance=vecnorm(xyz,2,2);
couplings=hyperfine*(1-3*((xyz*field_axis)./distance).^2)./distance.^3;

% Evaluate the same secular model used by the compact Hahn context
n_spins=size(xyz,1);
echo=cce_pair_product(xyz,field_axis,couplings,dipolar,sparse(n_spins,n_spins),...
                      times,abs_tol,Inf);

end

% Consistency enforcement
function grumble(xyz,field_axis,hyperfine,dipolar,times,abs_tol)
if ~isnumeric(xyz)||~isreal(xyz)||~ismatrix(xyz)||size(xyz,2)~=3||any(~isfinite(xyz(:)))
    error('xyz must be a finite real N-by-3 array.');
end
if any(all(xyz==0,2))||size(unique(xyz,'rows'),1)~=size(xyz,1)
    error('xyz must contain distinct sites away from the central spin.');
end
if ~isnumeric(field_axis)||~isreal(field_axis)||~isequal(size(field_axis),[3 1])||...
   ~(abs(norm(field_axis)-1)<=4*eps)
    error('field_axis must be a real three-element unit column vector.');
end
if ~isnumeric(hyperfine)||~isreal(hyperfine)||~isscalar(hyperfine)||~isfinite(hyperfine)
    error('hyperfine must be a finite real scalar.');
end
if ~isnumeric(dipolar)||~isreal(dipolar)||~isscalar(dipolar)||~isfinite(dipolar)
    error('dipolar must be a finite real scalar.');
end
if ~isnumeric(times)||~isreal(times)||~iscolumn(times)||isempty(times)||...
   any(~isfinite(times))||any(times<0)
    error('times must be a non-empty finite non-negative column vector.');
end
if ~isnumeric(abs_tol)||~isreal(abs_tol)||~isscalar(abs_tol)||~isfinite(abs_tol)||abs_tol<0
    error('abs_tol must be a finite non-negative real scalar.');
end
end


