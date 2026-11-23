% Secular distant dipolar field of a periodic one-axis magnetisation
% grating in an axially symmetric sample. Syntax:
%
%              field=ddf_field(magn,direction,demag_factor)
%
% Parameters:
%
%           magn - real N-by-3 Cartesian magnetisation array in A/m,
%                  sampled at equally spaced positions over whole periods
%
%      direction - unit three-vector along the spatial modulation axis
%                  in the laboratory frame (the main field is along Z)
%
%   demag_factor - axial SI demagnetising factor of the sample, in [0,1]
%
% Outputs:
%
%          field - N-by-3 Cartesian rotating-frame field array, tesla
%
% The spatial mean sees the sample shape; all nonzero Fourier modes see
% the modulation direction. This is not a finite-object field solver.
% See https://doi.org/10.1063/1.4721637 for the underlying physics.

function field=ddf_field(magn,direction,demag_factor)

% Check consistency
grumble(magn,direction,demag_factor);

% Separate the sample-shape and modulation contributions
uniform=mean(magn,1);
bulk=(3*demag_factor-1)/2;
grating=(3*direction(3)^2-1)/2;

% Apply the secular dipolar tensor using SI vacuum permeability
field=(1.25663706127e-6/3)*(bulk*uniform+grating*(magn-uniform)).*[1 1 -2];

end

% Consistency enforcement
function grumble(magn,direction,demag_factor)
if (~isfloat(magn))||(~isreal(magn))||(~ismatrix(magn))||...
   isempty(magn)||(size(magn,2)~=3)||any(~isfinite(magn(:)))
    error('magn must be a finite real floating-point N-by-3 array.');
end
if (~isfloat(direction))||(~isreal(direction))||(~isvector(direction))||...
   (numel(direction)~=3)||any(~isfinite(direction))||...
   (abs(norm(direction)-1)>3*eps(class(direction)))
    error('direction must be a finite real unit three-vector.');
end
if (~isfloat(demag_factor))||(~isreal(demag_factor))||(~isscalar(demag_factor))||...
   (~isfinite(demag_factor))||(demag_factor<0)||(demag_factor>1)
    error('demag_factor must be a real scalar between zero and one.');
end
end


