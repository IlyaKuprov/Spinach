% Irreducible spherical tensor expansion of a user-specified
% bosonic operator product. Syntax:
%
%          [states,coeffs]=bos2ist(prod_spec,lvl_num)
%
% Parameters:
%
%     prod_spec - bosonic operator product specification
%                 in which 'C' stands for creation opera-
%                 tor and 'A' for annihilation operator,
%                 for example 'CCAA'
%
%     nlevels   - number of energy levels in the trunca-
%                 ted bosonic mode
%
% Outputs:
%
%     states  - states, in the Spinach IST basis index-
%               ing, that contribute to the operator in
%               question; use lin2lm to convert to L,M
%               spherical tensor indices
%
%     coeffs  - coefficients with which the ISTs enter
%               the linear combination
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=bos2ist.m>

function [states,coeffs]=bos2ist(prod_spec,nlevels)

% Check consistency
grumble(prod_spec,nlevels);

% Weyl operators
W=weyl(nlevels);

% Compute the product
P=speye(nlevels,nlevels);
for n=1:numel(prod_spec)
    switch prod_spec(n)

        case 'C'

            % Creation
            P=P*W.c;

        case 'A'

            % Annihilation
            P=P*W.a;

        case 'N'

            % Number
            P=P*W.n;

        otherwise

            % Complain and bomb out
            error('unknown operator specification.');

    end
end

% Spherical tensor expansion
[states,coeffs]=oper2ist(P);

end

% Consistency enforcement
function grumble(prod_spec,nlevels)
if ~char(prod_spec)
    error('prod_spec must be a character string.');
end
if (~isnumeric(nlevels))||(~isscalar(nlevels))||...
   (~isreal(nlevels))||(nlevels<1)
    error('nlevels must be a positive real integer.');
end
end

% Wenn sich der Most auch ganz absurd gebärdet,
% Es gibt zuletzt doch noch n’ Wein.
%
% Goethe, in Faust


