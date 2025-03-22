% Chebyshev phase corrector for 1D NMR spectra. Views the
% phase profile across the spectral window as a slowly va-
% rying function and approximates it with a linear combi-
% nation of low-order Chebyshev polynomials. Syntax:
%
%        [spec,cheb_coeffs]=autophase(spec,guess)
%
% Parameters:
%
%     spec - 1D NMR spectrum, a complex vector
%
%    guess - initial guess for the Chebyshev polynomi-
%            al coefficients, radians. [phi 0 0] is a 
%            good start, where phi is the zero-order
%            phase correction guess.
%
% Outputs:
%
%     spec - phased NMR spectrum, a column vector
%
%   coeffs - Chebyshev polynomial coefficients of the 
%            phase profile across the spectrum with
%            the window treated as a [-1,1] interval
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=autophase.m>

function [spec,cheb_coeffs]=autophase(spec,guess)

% Check consistency
grumble(spec,guess);

% Scale the spectrum by its standard deviation
scaling_factor=std(spec); spec=spec/scaling_factor;

% Shift the 4-norm from imaginary to real as much as possible
options=optimoptions('fminunc','Display','iter','GradObj','off',...
                     'MaxIterations',Inf,'FunctionTolerance',1e-12,...
                     'OptimalityTolerance',1e-12,'StepTolerance',1e-12,...
                     'MaxFunctionEvaluations',Inf,'DerivativeCheck','off',...
                     'FinDiffType','central');
cheb_coeffs=fminunc(@(x)objective(spec,x),guess,options);

% Apply the correction
spec=apply_phases(spec,cheb_coeffs);

% Undo spectrum scaling
spec=scaling_factor*spec;

end

% Objective function
function obj=objective(spec,phis)

% Apply the phases
spec=apply_phases(spec,phis);

% Move the 4-norm from imag to real
obj=norm(imag(spec),4)-norm(real(spec),4);

end

% Chebyshev phase profile applicator
function spec=apply_phases(spec,phis)

% Get Chebyshev polynomials
cheb=ones(numel(phis),numel(spec));
cheb(2,:)=linspace(-1,1,numel(spec));
for n=3:numel(phis)
    cheb(n,:)=2*cheb(2,:).*cheb(n-1,:)-cheb(n-2,:);
end

% Get phase multipliers
phi_mults=exp(1i*phis*cheb);

% Apply phase multipliers
spec=phi_mults(:).*spec(:);

end

% Consistency enforcement
function grumble(spec,guess)
if (~isnumeric(spec))||(~isvector(spec))
    error('spec must be a vector.');
end
if (~isnumeric(guess))||(~isreal(guess))||(~isvector(guess))
    error('guess must be a real vector.');
end
end

% Human beings are born with different capacities. If 
% they are free, they are not equal. And if they are
% equal, they are not free.
%
% Alexander Solzhenitsyn

