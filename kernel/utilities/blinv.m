% Blicharski's relaxation theory invariants, as given by Equations
% 20-21 in http://doi.org/10.1515/zna-1972-1012, with an error and
% a typo corrected in Equation 21. Syntax:
%
%                        [Lsq,Dsq]=blinv(A)
%
% where A is the interaction matrix. This function is not sensitive
% to the trace of the matrix. Parameters:
%
%     A   -   a real 3x3 matrix
%
% Outputs:
%
%     Lsq -   first rank invariant
%
%     Dsq -   second rank invariant
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=blinv.m>

function [Lsq,Dsq]=blinv(A)

% Check consistency
grumble(A);

% First rank invariant
Lsq=(A(1,2)-A(2,1))^2+...
    (A(1,3)-A(3,1))^2+...
    (A(2,3)-A(3,2))^2;

% Second rank invariant
Dsq=A(1,1)^2+...
    A(2,2)^2+...
    A(3,3)^2-...
    A(1,1)*A(2,2)-...
    A(1,1)*A(3,3)-...
    A(2,2)*A(3,3)+...
    (3/4)*((A(1,2)+A(2,1))^2+...
           (A(1,3)+A(3,1))^2+...
           (A(2,3)+A(3,2))^2);

end

% Consistency enforcement
function grumble(A)
if (~isnumeric(A))||(~isreal(A))||...
   (~ismatrix(A))||any(size(A)~=[3 3])
    error('A must be a real 3x3 matrix.');
end
end

% His fight came to an abrupt end with the realisation that truth is ir-
% relevant, and all that matters in society is belief. Those that see the
% truth outside accepted beliefs can suffer in silence or cry out in fu-
% tility one last time before they disappear.
%
% Tides of Numenera, a computer game

