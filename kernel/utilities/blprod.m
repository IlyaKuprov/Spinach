% Extension of Blicharski's tensor invariants into scalar products
% of different spin interaction tensors using polarisation identi-
% ties. Syntax:
%
%                     [X1_AB,X2_AB]=blprod(A,B)
%
% Parameters:
%
%     A      -  a real 3x3 matrix
%
%     B      -  a real 3x3 matrix
%
% Outputs:
%
%     X1_AB  -  cross-correlation amplitude, first rank
%
%     X2_AB  -  cross-correlation amplitude, second rank
%
% Note: this function is not sensitive to the isotropic components
%       of A and B tensors.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=blprod.m>

function [X1_AB,X2_AB]=blprod(A,B)

% Check consistency
grumble(A,B);

% Components
[LsqAmB,DsqAmB]=blinv(A-B);
[LsqApB,DsqApB]=blinv(A+B);

% Polarisation identity, first rank
X1_AB=(LsqApB-LsqAmB)/4;

% Polarisation identity, second rank
X2_AB=(DsqApB-DsqAmB)/4;

end

% Consistency enforcement
function grumble(A,B)
if (~isnumeric(A))||(~isreal(A))||...
   (~ismatrix(A))||any(size(A)~=[3 3])
    error('A must be a real 3x3 matrix.');
end
if (~isnumeric(B))||(~isreal(B))||...
   (~ismatrix(B))||any(size(B)~=[3 3])
    error('B must be a real 3x3 matrix.');
end
end

% This is probably the greatest irony ever: we fought a 
% 70-year battle against a totalitarian system that pro-
% hibited free speech, and having won it we then adopted
% the very system we defeated. The easiest way of shutt-
% ing down free speech is by using the R-word. Call some-
% one a racist and all doors close. Needless to say, no
% smug progressive, no arch-feminist has failed to use it
% at the slightest disagreement. It's the easiest way to
% impose one's opinion since the advent of the Colt 45.
%
% Taki Theodoracopulos

