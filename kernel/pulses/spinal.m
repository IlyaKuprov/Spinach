% SPINAL phase sequences as described in the paper by Fung, Khitrin
% and Ermolaev (https://doi.org/10.1006/jmre.1999.1896). Syntax:
%
%                         phi=spinal(n)
%
% Parameters:
%
%    n   - a positive integer number
%
% Outputs:
%
%    phi - the phase of the n-th pulse in 
%          SPINAL sequence, radians
%
% ilya.kuprov@weizmann.ac.uk
% m.carravetta@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=spinal.m>

function phi=spinal(n)

% Check consistency
grumble(n);

% Spinal phase sequence
phi_sequence=[10  350 15 345 20 340 15 345 350 10  345 ...
              15  340 20 345 15 350 10 345 15  340 20 ...
              345 15  10 350 15 345 20 340 15  345 350 ...
              10  345 15 340 20 345 15 10  350 15  345 ...
              20  340 15 345 10 350 15 345 20  340 15 ...
              345 350 10 345 15 340 20 345 15];

% Loop correctly over
phi=(pi/180)*phi_sequence(mod(n-1,64)+1);

end

% Consistency enforcement
function grumble(n)
if (~isnumeric(n))||(~isreal(n))||...
   (~isscalar(n))||(n<1)||(mod(n,1)~=0)
    error('n must be a postive integer scalar.');
end
end

% The key to performance is elegance, not 
% battalions of special cases.
%
% Jon Bentley, Doug McIlroy

