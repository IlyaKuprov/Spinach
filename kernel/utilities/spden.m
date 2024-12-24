% Lorentzian spectral density function for rotational 
% diffusion at the user-specified frequency. Syntax:
%
%                  J=spden(L,D,omega)
%
% Parameters:
%
%    L     - spherical rank, use 2 for common NMR
%            mechanisms such as dipolar relaxation
%
%    D     - rotational diffusion coefficient, s^{-1}
%
%    omega - frequency, rad/s
%
% Outputs:
%
%    J     - spectral density function value
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=spden.m>

function J=spden(L,D,omega)

% Check consistency
grumble(L,D,omega);

% Get the correlation time
tau_c=1/(L*(L+1)*D);

% Get the spectral density
J=(tau_c/(2*L+1))/(1+(tau_c*omega)^2);

end

% Consistency enforcement
function grumble(L,D,omega)
if (~isnumeric(L))||(~isscalar(L))||...
   (~isreal(L))||(L<1)||(mod(L,1)~=0)
    error('L must be a positive real integer.');
end
if (~isnumeric(D))||(~isscalar(D))||...
   (~isreal(D))||(D<=0)
    error('D must be a positive real scalar.');
end
if (~isnumeric(omega))||(~isscalar(omega))||...
   (~isreal(omega))
    error('omega must be a real scalar.');
end
end

% Sometimes the "rules" aren't really even rules. Gordon Bruce, 
% the former CIO of the city of Honolulu, explained to me that
% when he entered government from the private sector and tried
% to make changes, he was told, "That's against the law." His
% reply was "OK. Show me the law." "Well, it isn't really a law.
% It's a regulation." "OK. Show me the regulation." "Well, it
% isn't really a regulation. It's a policy that was put in place
% by Mr. Somebody twenty years ago." "Great. We can change that!"
% 
% Tim O'Reilly

