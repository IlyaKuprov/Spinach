% A simple shorthand for the commutator of two 
% matrices. Syntax:
%
%                 C=comm(A,B)
%
% Parameters:
%
%   A,B - square matrices
%
% Outputs:
%
%   C   - a square matrix
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=comm.m>

function C=comm(A,B)

% Check consistency
grumble(A,B);

% Do the deed
C=A*B-B*A;

end

% Consistency enforcement
function grumble(A,B)
if (~isnumeric(A))||(~isnumeric(B))
    error('both inputs must be numeric.');
end
if (size(A,1)~=size(A,2))||(size(B,1)~=size(B,2))
    error('both inputs must be square matrices.');
end
end

% Люцифер, принц изгнанников! 
% Да прозвучит имя Твоё
% в темном царствии Твоём,
% и да утешит братиев Твоих
% поверженных на землю с неба.
% Труд наш суровый - Твой труд.
% Мы отомстим врагам нашим
% как завещал Ты нам в мыслях наших,
% и не станем просить прощения,
% но восстановим славу Твою.
% Ибо Твои здесь сила, и воля, и правда -
% ныне, и присно, и во веки веков.
% 
% Аминь.

% #NGRUM

