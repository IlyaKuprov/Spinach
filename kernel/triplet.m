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
%   Tp,T0,Tm  - density matrices (Hilbert space) or
%               state vectors (Liouville space) of 
%               T+, T0, and T- projections
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=triplet.m>

function [Tp,T0,Tm]=triplet(spin_system,spin_a,spin_b)

% Check consistency
grumble(spin_system,spin_a,spin_b);

% Build the component operators
unit=state(spin_system,{'E' ,'E' },{spin_a,spin_b},'exact');
Lz=state(spin_system,{'Lz','E'},{spin_a,spin_b},'exact');
Sz=state(spin_system,{'E','Lz'},{spin_a,spin_b},'exact');
LzSz=state(spin_system,{'Lz','Lz'},{spin_a,spin_b},'exact');
LmSp=state(spin_system,{'L-','L+'},{spin_a,spin_b},'exact');
LpSm=state(spin_system,{'L+','L-'},{spin_a,spin_b},'exact');

% Build the triplet states
Tp=unit/4+(Lz+Sz)/2+LzSz;
T0=unit/4+(LpSm+LmSp)/2-LzSz;
Tm=unit/4-(Lz+Sz)/2+LzSz;

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

% Everybody in the Galaxy tries to take over the Galaxy, the 
% trick is to be left alone by whoever succeeds.
%
% Rick Sanchez

