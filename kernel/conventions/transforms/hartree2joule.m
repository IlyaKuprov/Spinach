% Converts Hartree energy units into J/mol. A Hartree is twice the 
% ground state ionisation energy of the hydrogen atom. Syntax:
%
%           energy=hartree2joule(energy)
%
% Parameters:
%
%     energy  - a numerical array of energies in
%               Hartree units
%
% Outputs:
%
%     energy  - a numerical array of energies in 
%               Joules
%
% ilya.kuprov@weizmann.ac.uk
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=hartree2joule.m>

function energy=hartree2joule(energy)

% Check consistency
grumble(energy);

% Perform the conversion
energy=2625499.62*energy;

end

% Consistency enforcement
function grumble(energy)
if (~isnumeric(energy))||(~isreal(energy))
    error('the argument must be an array of real numbers.');
end
end

% "We only need to be lucky once. You need to be
%  lucky every time."
%
% The IRA to Margaret Thatcher, after 
% a failed assassination attempt.

