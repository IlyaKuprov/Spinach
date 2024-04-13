% PMLG5 phase sequence as described in the paper by Vinogradova,
% Madhu and Vega (https://doi.org/10.1016/S0009-2614(99)01174-4).
% Syntax:
%                            phi=spinal(n)
%
% Parameters:
%
%    n   - a positive integer number
%
% Outputs:
%
%    phi - the phase of the n-th pulse in 
%          PMLG sequence, radians
%
% i.kuprov@soton.ac.uk
% m.carravetta@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=pmlg5.m>

function phi=pmlg5(n)

% Check consistency
grumble(n);

% PMLG5 phase sequence
phi_sequence=[339.22 297.65 256.08 214.51 172.94 ...
              352.94  34.51  76.08 117.65 159.22 ...
              159.22 117.65  76.08  34.51 352.94 ...
              172.94 214.51 256.08 297.65 339.22];

% Loop correctly over
phi=(pi/180)*phi_sequence(mod(n-1,20)+1);

end

% Consistency enforcement
function grumble(n)
if (~isnumeric(n))||(~isreal(n))||...
   (~isscalar(n))||(n<1)||(mod(n,1)~=0)
    error('n must be a postive integer scalar.');
end
end

% One man's crappy software is another 
% man's full time job.
%
% Jessica Gaston

