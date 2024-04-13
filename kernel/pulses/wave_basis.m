% Common basis sets for the expansion of pulse waveforms. Returns the wave-
% form basis functions as columns of a matrix. Syntax:
%
%           basis_waves=wave_basis(basis_type,n_functions,n_steps)
%
% Parameters:
%
%       basis_type     - may be set to 'sine_waves', 'cosine_waves',
%                        and 'legendre'. The sine and the cosine op-
%                        tions return the corresponding functions in
%                        the [-pi,pi] interval, legendre option re-
%                        turns legendre polynomials in the [-1,1] in-
%                        terval.
%
%       n_func         - the number of functions to return (integer
%                        frequencies starting from zero on the case
%                        of cosines, integer frequencies starting 
%                        from 1 inthe case of sines, legendre poly-
%                        nomial ranks in the case of legendre func-
%                        tion basis set.
%
%       n_points       - number of discretization points.
%
% Outputs:
%
%       basis_waves    - a matrix with the basis waves in columns
%
% Note: because the resulting waveforms are discretised, they are not pre-
%       cisely orthogonal under the standard scalar multiplication. An ex-
%       tra orthogonalisation step is therefore applied to make them ortho-
%       gonal as vectors. As a result, some functions may be upside-down.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=wave_basis.m>

function basis_waves=wave_basis(basis_type,n_func,n_points)

% Check consistency
grumble(basis_type,n_func,n_points);

% Preallocate the array
basis_waves=zeros(n_func,n_points);

% Generate the functions
switch basis_type
    
    case 'sine_waves'
        
        % Fill the array
        parfor n=1:n_func
            basis_waves(n,:)=sin(n*linspace(-pi,pi,n_points));
        end

    case 'cosine_waves'

        % Fill the array
        parfor n=1:n_func
            basis_waves(n,:)=cos((n-1)*linspace(-pi,pi,n_points));
        end
        
    case 'legendre'
        
        % Fill the array
        parfor n=1:n_func
            basis_waves(n,:)=legendreP(n-1,linspace(-1,1,n_points));
            basis_waves(n,:)=basis_waves(n,:)/norm(basis_waves(n,:),2);
        end
        
    otherwise
        
        % Complain and bomb out
        error('unrecognized basis function type.');
        
end

% Orthogonalize the array
basis_waves=orth(basis_waves');

% Check outgoing dimension
if size(basis_waves,2)~=n_func
    error('basis is linearly dependent, reduce n_func.');
end
        
end

% Consistency enforcement
function grumble(basis_type,n_func,n_points)
if ~ischar(basis_type)
    error('basis_type parameter must be a character string.');
end
if (numel(n_func)~=1)||(~isnumeric(n_func))||(~isreal(n_func))||...
   (n_func<1)||(mod(n_func,1)~=0)
    error('n_func parameter must be a positive real integer greater than 1.');
end
if (numel(n_points)~=1)||(~isnumeric(n_points))||(~isreal(n_points))||...
   (n_points<1)||(mod(n_points,1)~=0)
    error('n_points parameter must be a positive real integer greater than 1.');
end
end

% Self-esteem is the reputation we acquire with ourselves. 
%
% Nathaniel Branden

