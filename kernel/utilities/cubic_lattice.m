% Creates a periodic volume-centered cubic lattice with user- 
% supplied parameters. Syntax:
%
%    [sys,inter]=cubic_lattice(isotope,spacing,n_periods)
%
% Parameters:
%
%   isotope   - character string specifying the isotope,
%               for example, '13C'
%
%   spacing   - lattice spacing in Angstroms
%
%   n_periods - number of lattice periods in each of the
%               three spatial dimensions
%
% Outputs:
%
%   sys, inter - Spinach input data structures with the
%                following fields set:
%
%                sys.isotopes, inter.coordinates, inter.pbc
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=cubic_lattice.m>

function [sys,inter]=cubic_lattice(isotope,spacing,n_periods)

% Check consistency
grumble(isotope,spacing,n_periods);

% Generate the isotope list
sys.isotopes=cell(1,n_periods^3);
for n=1:n_periods^3
    sys.isotopes{n}=isotope;
end

% Generate coordinates
inter.coordinates=cell(n_periods^3,1);
for n=1:n_periods
    for k=1:n_periods
        for m=1:n_periods
            inter.coordinates{sub2ind([n_periods n_periods n_periods],m,k,n)}=spacing*[(n-1) (k-1) (m-1)];
        end
    end
end

% Generate lattice translation vectors
inter.pbc={spacing*n_periods*[1 0 0],...
           spacing*n_periods*[0 1 0],...
           spacing*n_periods*[0 0 1]};
       
end

% Consistency enforcement
function grumble(isotope,spacing,n_periods)
if ~ischar(isotope)
    error('isotope must be a character string.');
end
if (~isnumeric(spacing))||(~isreal(spacing))||...
   (numel(spacing)~=1)||(~isfinite(spacing))||(spacing<=0)
    error('spacing must be a positive real number.');
end
if (~isnumeric(n_periods))||(~isreal(n_periods))||...
   (numel(n_periods)~=1)||(~isfinite(n_periods))||...
   (n_periods<1)||(mod(n_periods,1)~=0)
    error('n_periods must be a positive real integer.');
end
end

% There are only two ways of telling the complete 
% truth - anonymously and posthumously.
%
% Thomas Sowell

