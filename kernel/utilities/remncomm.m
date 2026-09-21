% Removes from the Hermitian operator A the part that does not com-
% mute with the Hermitian operator B. Syntax:
%
%                       C=remncomm(A,EvecB,EvalB)
%
% Parameters:
%
%    A     -  a square matrix
%
%    EvecB -  a square matrix containing eigenvectors
%             of B in columns
%
%    EvalB -  a column vector containing the eigenvalues
%             of B in the same order as the columns of EvecB
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

function A=remncomm(A,EvecB,EvalB)

% Check consistency
grumble(A,EvecB,EvalB);

% Move A into the eigenbasis of B
A=EvecB'*A*EvecB;

% Zero out elements linking eigenvalues of B that differ by more than eigensolver roundoff
degen_mask=abs(EvalB-EvalB.')<=100*numel(EvalB)*eps(max(EvalB)-min(EvalB));
A=A.*degen_mask;

% Move the commuting part back into the original basis
A=EvecB*A*EvecB';

end

% Consistency enforcement
function grumble(A,EvecB,EvalB)
if (~isnumeric(A))||(size(A,1)~=size(A,2))||...
   (~ishermitian(A))
    error('A must be a Hermitian matrix.');
end
if (~isnumeric(EvecB))||(size(EvecB,1)~=size(EvecB,2))
    error('EvecB must be a square array of column vectors.');
end
if (~isfloat(EvalB))||(~iscolumn(EvalB))||(~isreal(EvalB))||...
   (~all(isfinite(EvalB)))||(numel(EvalB)~=size(EvecB,2))
    error('EvalB must be a finite real floating-point column vector with as many elements as EvecB has columns.');
end
if ~isfinite(max(EvalB)-min(EvalB))
    error('the spread of EvalB must be representable in floating point.');
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

