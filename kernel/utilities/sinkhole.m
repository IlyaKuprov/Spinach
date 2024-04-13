% Turns the specified states into sinkholes -- any population reaching
% them will be summed up and stored forever in a frozen state. This is 
% useful for state space restriction diagnostics. Syntax:
%
%                   L=sinkhole(spin_system,L,states)
%
% Parameters:
%
%      L       - Liovillian matrix
%
%      states  - a vector of integers specifying the 
%                numbers of the states to be set up
%                as sinkholes
%
% Output:
%
%      L       - updated Liouvillian matrix
%
% Note: this functionality is only available in sphten-liouv formalism.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=sinkhole.m>

function L=sinkhole(spin_system,L,states)

% Check consistency
grumble(spin_system,L,states);

% Zero columns corresponding to sinkhole states
L(:,states)=0;

end

% Consistency enforcement
function grumble(spin_system,L,states)
if ~strcmp(spin_system.bas.formalims,'sphten-liouv')
    error('this function is only applicable to sphten-liouv formalism.');
end
if (~isnumeric(L))||(size(L,1)~=size(L,2))
    error('L must be a square matrix.');
end
if any(size(L)~=size(spin_system.bas.basis,1))
    error('dimension of the Liouvillian must match the dimension of the basis set.');
end
if (~isnumeric(states))||(~isvector(states))||(~isreal(states))||...
   any(mod(states,1)~=0)||any(states<=0)
    error('states must be a vector of positive integers.');
end
if any(states>size(spin_system.bas.basis,1))
    error('an element of the states vector exceeds the state space dimension.');
end
end

% Acting in a bad play is like spitting into eternity.
%
% Faina Ranevskaya

