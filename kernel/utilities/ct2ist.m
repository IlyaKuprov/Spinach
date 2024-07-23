% Irreducible spherical tensor expansion of central transition
% operators of half-integer spins. Syntax:
%
%               [states,coeffs]=ct2ist(mult,type)
%
% Parameters:
%
%     mult   - multipicity of the spin in question, an
%              even positive integer
%
%     type   - operator type: 'z' for polarisation, '+'
%              for raising, '-' for lowering
%
% Outputs:
%
%     states - states, in the Spinach IST basis index-
%              ing, that contribute to the operator in
%              question; use lin2lm to convert to L,M
%              spherical tensor indices
%
%     coeffs - coefficients with which the ISTs enter
%              the linear combination
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ct2ist.m>

function [states,coeffs]=ct2ist(mult,type)

% Check consistency
grumble(mult,type);

% CT operator in the Zeeman basis
CT_Z=centrans(mult,type);

% Spherical tensors in Zeeman basis
IST_Z=irr_sph_ten(mult);

% Compute all expansion coefficients and list all states
coeffs=cellfun(@(A)trace(A'*CT_Z)/trace(A'*A),IST_Z);
states=transpose(0:(numel(IST_Z)-1));

% Drop negligible states
idx=(abs(coeffs)>10*eps('double'));
coeffs=coeffs(idx); states=states(idx);

end

% Consistency enforcement
function grumble(mult,type)
if (~isnumeric(mult))||(~isscalar(mult))||...
   (~isreal(mult))||(mult<2)||(mod(mult,2)~=0)
    error('mult must be an even positive integer.');
end
if (~ischar(type))||(~ismember(type,{'x','y','z','+','-'}))
    error('type must be ''x'', ''y'', ''z'', ''+'', or ''-''.');
end
end

% A hundred miles is not a detour to a rabid dog.
%
% A Russian saying

