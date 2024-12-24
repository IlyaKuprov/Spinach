% Sorensen bound for the maximum transfer efficiency between
% two states under arbitrary control operators. Equation 186
% from https://doi.org/10.1016/0079-6565(89)80006-8. Syntax:
%
%               b=sorensen(rho_init,rho_targ)
%
% Parameters:
%
%    rho_init  - initia ldensity matrix, Hilbert space
%
%    rho_targ  - target density matrix, Hilbert space
%
% Output:
%
%    b         - Sorensen bound
%
% Note: this is an exact unitary bound; the amount reachable
%       with realistically available instrumental controls
%       may be smaller, see the detailed analysis here:
%       https://doi.org/10.1080/00268979909483117
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=sorensen.m>

function b=sorensen(rho_init,rho_targ)

% Check consistency
grumble(rho_init,rho_targ);

% Diagonalise both matrices
[~,sigma_init]=eig(full(rho_init));
[~,sigma_targ]=eig(full(rho_targ));

% Sort eigenvalues
sigma_init=sort(diag(sigma_init));
sigma_targ=sort(diag(sigma_targ));

% Compute Sorensen bound
b=(sigma_init'*sigma_targ)/trace(rho_targ^2);

end

% Consistency enforcement
function grumble(rho_init,rho_targ)
if (~isnumeric(rho_init))||(~isnumeric(rho_targ))||...
   (~ishermitian(rho_init))||(~ishermitian(rho_targ))||...
   (numel(rho_init)~=numel(rho_targ))
   error('the inputs must be Hermitian matrices of the same size.');
end

end

% People think they do not understand mathematics, but that rather
% depends on how you explain. If you ask an alcoholic which is big-
% ger, 2/3 or 3/5, he would not be able to say. But if you reformu-
% late the question and ask which is better, two bottles of vodka
% for three people or three bottles for five people, he would imme-
% diately reply: "of course it's two bottles for three people".
%
% Israel Gelfand

