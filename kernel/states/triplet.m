% Returns the components of the two-spin triplet state. Syntax:
%
%         [Tp,T0,Tm]=triplet(spin_system,spin_a,spin_b)
%
% Arguments:
%
%   spin_a    - the number of the first spin in the 
%               triplet state
%
%   spin_b    - the number of the second spin in the 
%               triplet state
%
% Outputs:
%
%   TU,T0,TD  - density matrices (Hilbert space) or
%               state vectors (Liouville space) of 
%               TU, T0, and TD projections
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=triplet.m>

function [TU,T0,TD]=triplet(spin_system,spin_a,spin_b)

% Check consistency
grumble(spin_system,spin_a,spin_b);

% Build the component operators
EE=state(spin_system,{'E' ,'E' },{spin_a,spin_b});
ZE=state(spin_system,{'Lz','E'},{spin_a,spin_b});
EZ=state(spin_system,{'E','Lz'},{spin_a,spin_b});
XX=state(spin_system,{'Lx','Lx'},{spin_a,spin_b});
YY=state(spin_system,{'Ly','Ly'},{spin_a,spin_b});
ZZ=state(spin_system,{'Lz','Lz'},{spin_a,spin_b});

% Build the triplet states
TU=EE/4+(ZE+EZ)/2+ZZ;  % up
T0=EE/4+XX+YY-ZZ;      % middle
TD=EE/4-(ZE+EZ)/2+ZZ;  % down

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

% Everybody in the Galaxy tries to take over the Galaxy, the 
% trick is to be left alone by whoever succeeds.
%
% Rick Sanchez

