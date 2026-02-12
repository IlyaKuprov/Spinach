% Irreducible spherical tensor operator expansion of a user-
% specified square matrix. Syntax:
%
%                [states,coeffs]=oper2ist(A)
%
% Parameters:
%
%     A      - a square matrix
%
% Outputs:
%
%     states - states, in the Spinach IST basis index-
%              ing convention, that contribute to the
%              operator in question; use lin2lm() fun-
%              ction to convert to L,M spherical tens-
%              or indices
%
%     coeffs - coefficients with which the ISTs enter
%              the linear combination
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=oper2ist.m>

function [states,coeffs]=oper2ist(A)

% Check consistency
grumble(A);

% Spherical tensors
T=irr_sph_ten(size(A,1));

% Get expansion coefficients and states
coeffs=cellfun(@(X)full(hdot(X,A)/hdot(X,X)),T);
states=transpose(0:(numel(T)-1));

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

% The obedient always think of themselves 
% as virtuous rather than cowardly.
%
% George Carlin

