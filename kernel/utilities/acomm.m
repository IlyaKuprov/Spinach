% A simple shorthand for the anticommutator of two 
% matrices. Syntax:
%
%                   C=acomm(A,B)
%
% Parameters:
%
%   A,B - square matrices
%
% Outputs:
%
%   C   - a square matrix
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=acomm.m>

function C=acomm(A,B)

% Do the deed
C=A*B+B*A;

end

% Enough of all this academic chatter, back
% again to devilry!
%
% Mephisto

% #NGRUM