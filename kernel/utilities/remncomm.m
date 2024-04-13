% Removes from the Hermitian operator A the part that does not com-
% mute with the Hermitian operator B. Syntax:
%
%                         C=remncomm(A,EvB)
%
% Parameters:
%
%    A     -  a square matrix
%
%    EvB   -  a square matrix containing eigenvectors 
%             of B in columns
%
% Outputs:
%
%    C     -  a square matrix
%
% Note: when the matrix B is diagonal, the part of A that commutes
%       with it is the diagonal part, so just use diag(diag(A))
%       
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=remncomm.m>

function A=remncomm(A,EvB)

% Check consistency
grumble(A,EvB);

% The standard projector expression
A=EvB*(diag(diag(EvB'*A*EvB)))*EvB';

end

% Consistency enforcement
function grumble(A,EvB)
if (~isnumeric(A))||(size(A,1)~=size(A,2))||...
   (~ishermitian(A))
    error('A must be a Hermitian matrix.');
end
if (~isnumeric(EvB))||(size(EvB,1)~=size(EvB,2))
    error('EvB must be a square array of column vectors.');
end
end

% The first scientific measurement of the speed of electricity was 
% conducted in 1764 by French physicist Jean-Antoine Nollet. He ar-
% ranged two hundred monks into a large circle, and connected their
% hands with iron wire. He then discharged a large Leyden Jar bat-
% tery into the wire. Nollet was unable to measure the actual speed
% of electricity because all monks reacted simultaneously. He noted
% that the transmission speed of electricity was very high. Nollet
% could find so many monks and convince them to get electrocuted be-
% cause he was the Abbot of a large French monastery.

