% A floating-point precision consistent check for whether
% a particular matrix is traceless. Syntax:
%
%                     A=istraceless(M)
% 
% Parameters:
%
%     M - a matrix of any dimension
%
% Outputs:
%
%     A - true if the matrix is traceless to ap-
%         propriate precision, false otherwise
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=istraceless.m>

function A=istraceless(M)

% Check consistency
grumble(M);

% Working precision
precision=eps(class(M));

% Cheapest norm of M
norm_m=cheap_norm(M);

% Decide if M is traceless
A=(abs(trace(M))<precision*norm_m);

end

% Consistency enforcement
function grumble(M)
if ~isnumeric(M)
    error('M must be numeric.');
end
end

% The College asked me to chair the Size and Shape
% Committee. My wife could not stop laughing.
%
% Oliver Taplin

