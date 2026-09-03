% Removes from the Hermitian operator A the part that does not com-
% mute with the Hermitian operator B. Syntax:
%
%                       C=remncomm(A,EvB,evals_b)
%
% Parameters:
%
%    A     -  a square matrix
%
%    EvB   -  a square matrix containing eigenvectors 
%             of B in columns
%
%    evals_b - a column vector containing the eigenvalues
%              of B in the same order as the columns of EvB
%
% Outputs:
%
%    C     -  a square matrix
%
% Note: within a degenerate eigenspace of B, every Hermitian operator
%       supported on that eigenspace commutes with B, so the corres-
%       ponding block of A (not just its diagonal) is kept
%       
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=remncomm.m>

function A=remncomm(A,EvB,evals_b)

% Check consistency
grumble(A,EvB,evals_b);

% Move A into the eigenbasis of B
A=EvB'*A*EvB;

% Zero out elements linking non-degenerate eigenvalues of B
degen_mask=abs(evals_b-evals_b.')<=1e-10*max(abs(evals_b));
A=A.*degen_mask;

% Move the commuting part back into the original basis
A=EvB*A*EvB';

end

% Consistency enforcement
function grumble(A,EvB,evals_b)
if (~isnumeric(A))||(size(A,1)~=size(A,2))||...
   (~ishermitian(A))
    error('A must be a Hermitian matrix.');
end
if (~isnumeric(EvB))||(size(EvB,1)~=size(EvB,2))
    error('EvB must be a square array of column vectors.');
end
if (~isnumeric(evals_b))||(~iscolumn(evals_b))||(~isreal(evals_b))||...
   (~all(isfinite(evals_b)))||(numel(evals_b)~=size(EvB,2))
    error('evals_b must be a finite real column vector with as many elements as EvB has columns.');
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

