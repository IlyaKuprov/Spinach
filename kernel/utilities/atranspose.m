% Anti-diagonal array transpose. Syntax:
%
%                       M=atranspose(M)
%
% Parameters:
%
%    M - a transposable array
%
% Outputs:
%
%    M - a transposable array
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=atranspose.m> 

function M=atranspose(M)

% Check consistency
grumble(M);

% Rotate and transpose
M=transpose(rot90(M,2));

end

% Consistency enforcement
function grumble(M)
if ~isnumeric(M)
    error('M must be a numeric array.');
end
end

% Overheard at the reception following David Deutch's lecture on
% Constructor Theory at the Oxford Physics Department in 2012:
%
% A - "It is not very often that you see so clearly
%      what is wrong with modern physics."
%
% B - "And what would that be?"
% 
% A - "The existence of this man. The possibility
%      of his existence."

