% Irreducible spherical tensor expansion of the density matrix of a
% coherent state of a truncated bosonic mode. Syntax:
%
%        [states,coeffs,lost_weight]=coh2ist(alpha,nlevels)
%
% Parameters:
%
%     alpha   - coherent state amplitude, a complex scalar
%
%     nlevels - number of energy levels in the truncated
%               bosonic mode
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
%     lost_weight - the fraction of the Poisson distribution
%                   norm chopped off by the mode truncation
%
% Note: the Fock space truncation of the mode chops the tail of
%       the Poisson distribution; the state is renormalised af-
%       ter the truncation, and the lost weight is returned.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=coh2ist.m>

function [states,coeffs,lost_weight]=coh2ist(alpha,nlevels)

% Check consistency
grumble(alpha,nlevels);

% Fock state amplitudes up to the truncation level
amps=alpha.^(0:(nlevels-1))./sqrt(factorial(0:(nlevels-1)));

% Poisson distribution weight lost to the truncation
lost_weight=1-exp(-abs(alpha)^2)*sum(abs(amps).^2);

% Normalise the truncated state
amps=amps/norm(amps,2);

% Spherical tensor expansion of the mode density matrix
[states,coeffs]=oper2ist(amps.'*conj(amps));

end

% Consistency enforcement
function grumble(alpha,nlevels)
if (~isnumeric(alpha))||(~isscalar(alpha))||(~isfinite(alpha))
    error('alpha must be a finite complex scalar.');
end
if (~isnumeric(nlevels))||(~isscalar(nlevels))||(~isreal(nlevels))||...
   (mod(nlevels,1)~=0)||(nlevels<1)
    error('nlevels must be a positive real integer.');
end
end

