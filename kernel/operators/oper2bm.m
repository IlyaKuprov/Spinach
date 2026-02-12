% Bosonic monomial operator expansion of a user-specified 
% square matrix. Syntax:
%
%                [states,coeffs]=oper2bm(A)
%
% Parameters:
%
%     A      - a square matrix
%
% Outputs:
%
%     states - states, in the Spinach BM basis index-
%              ing convention, that contribute to the
%              operator in question; use lin2kq() fun-
%              ction to convert to K,Q bosonic monomi-
%              al indices
%
%     coeffs - coefficients with which the BMs enter
%              the linear combination
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=oper2bm.m>

function [states,coeffs]=oper2bm(A)

% Check consistency
grumble(A);

% Bosonic monomials
BM=boson_mono(size(A,1));
nlevels=size(A,1);

% Compute the overlap matrix
S=zeros(nlevels^2,nlevels^2);
for n=1:nlevels^2
    for k=1:nlevels^2
        S(n,k)=hdot(BM{n},BM{k});
    end
end

% Get the expansion coefficients
coeffs=zeros(nlevels^2,1);
for n=1:nlevels^2
    coeffs(n)=hdot(BM{n},A);
end
coeffs=S\full(coeffs);

% List the states
states=transpose(0:(numel(BM)-1));

% Drop negligible states
idx=(abs(coeffs)>10*eps('double'));
coeffs=coeffs(idx); states=states(idx);

end

% Consistency enforcement
function grumble(A)
if (~isnumeric(A))||(~ismatrix(A))||...
   (size(A,1)~=size(A,2))
    error('A must be a square matrix.');
end
end

% The greatest trainers can teach 
% important skills, but the best 
% instructor of all may be pain.
%
% Frank Herbert, in the Dune series
