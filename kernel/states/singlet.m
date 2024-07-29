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
%   S       - a density matrix (Hilbert space) or
%             a state vector (Liouville space) 
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=singlet.m>

function S=singlet(spin_system,spin_a,spin_b)

% Check consistency
grumble(spin_system,spin_a,spin_b);

% Build the component operators
EE=state(spin_system,{'E' ,'E' },{spin_a,spin_b});
XX=state(spin_system,{'Lx','Lx'},{spin_a,spin_b});
YY=state(spin_system,{'Ly','Ly'},{spin_a,spin_b});
ZZ=state(spin_system,{'Lz','Lz'},{spin_a,spin_b});

% Build the singlet state
S=EE/4-(XX+YY+ZZ); 

end

% Consistency enforcement
function grumble(spin_system,spin_a,spin_b)
if (~isnumeric(spin_a))||(~isnumeric(spin_b))||...
   (~isscalar(spin_a))||(~isscalar(spin_b))||...
   (mod(spin_a,1)~=0)||(mod(spin_b,1)~=0)||...
   (spin_a<1)||(spin_b<1)||(spin_a==spin_b)
    error('spin indices must be different positive integers.');
end
if (spin_a>spin_system.comp.nspins)||...
   (spin_b>spin_system.comp.nspins)
    error('spin index exceeds the number of spins.');
end
end

% Physics is, hopefully, simple. Physicists are not.
% 
% Edward Teller

