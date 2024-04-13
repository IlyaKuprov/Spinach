% Returns a two-spin singlet state. Syntax:
%
%          rho=singlet(spin_system,spin_a,spin_b)
%
% Arguments:
%
%   spin_a  - the number of the first spin in the 
%             singlet state
%
%   spin_b  - the number of the second spin in the 
%             singlet state
%
% Outputs:
%
%   rho     - a density matrix (Hilbert space) or
%             a state vector (Liouville space) 
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=singlet.m>

function rho=singlet(spin_system,spin_a,spin_b)

% Check consistency
grumble(spin_system,spin_a,spin_b);

% Build the component operators
unit=state(spin_system,{'E' ,'E' },{spin_a,spin_b},'exact');
LzSz=state(spin_system,{'Lz','Lz'},{spin_a,spin_b},'exact');
LmSp=state(spin_system,{'L-','L+'},{spin_a,spin_b},'exact');
LpSm=state(spin_system,{'L+','L-'},{spin_a,spin_b},'exact');

% Build the singlet state
rho=unit/4-(LzSz+0.5*(LpSm+LmSp)); 

end

% Consistency enforcement
function grumble(spin_system,spin_a,spin_b)
if (~isnumeric(spin_a))||(~isnumeric(spin_b))||(~isscalar(spin_a))||...
   (~isscalar(spin_b))||(mod(spin_a,1)~=0)||(mod(spin_b,1)~=0)||...
   (spin_a<1)||(spin_b<1)||(spin_a==spin_b)
    error('spin_a and spin_b must be different positive real integers.');
end
if (spin_a>spin_system.comp.nspins)||(spin_b>spin_system.comp.nspins)
    error('spin_a and spin_b must be smaller or equal to the number of spins in the system.');
end
end

% Physics is, hopefully, simple. Physicists are not.
% 
% Edward Teller

