% Chebyshev phase corrector for 1D NMR spectra.
%
% i.kuprov@soton.ac.uk

function [spec,cheb_coeffs]=autophase(spec,guess)

% Scale by standard deviation
spec=spec/std(spec);

% Minimise the squared norm of the imaginary part
options=optimoptions('fminunc','Display','iter','GradObj','off',...
                     'MaxIterations',Inf,'FunctionTolerance',1e-12,...
                     'OptimalityTolerance',1e-12,'StepTolerance',1e-12,...
                     'MaxFunctionEvaluations',Inf,'DerivativeCheck','off',...
                     'FinDiffType','central');
cheb_coeffs=fminunc(@(x)objective(spec,x),guess,options);

% Apply the correction
spec=apply_phases(spec,cheb_coeffs);

end

% Objective function
function obj=objective(spec,phis)

% Apply the phases
spec=apply_phases(spec,phis);

% Move the 4-norm from imag to real
obj=norm(imag(spec),4)-norm(real(spec),4);

end

% Chebyshev phase applicator
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

% Human beings are born with different capacities. If 
% they are free, they are not equal. And if they are
% equal, they are not free.
%
% Alexander Solzhenitsyn

