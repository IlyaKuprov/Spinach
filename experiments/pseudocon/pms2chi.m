% Runs a least squares fitting procedure on top of Equation 10 from
% http://dx.doi.org/10.1039/c4cp03106g to extract the susceptibility
% tensor from DFT hyperfine coupling tensors and experimentally ob-
% served paramagnetic (contact + pseudocontact) shifts. Syntax:
% 
%              chi=pms2chi(hfcs,shifts,isotopes,nel)
%
% Parameters:
%
%        hfcs   - cell array of 3x3 hyperfine coupling tensors
%                 in Gauss
%
%      shifts   - vector of the observed pseudocontact shifts,
%                 excluding the diamagnetic contribution, ppm
%
%    isotopes   - cell array of character strings specifying
%                 isotopes that exhibit each of the chemical
%                 shifts supplied, for example {'1H','13C'}
%
%         nel   - number of unpaired electrons involved
%
% Outputs:
%
%         chi   - the fitted magnetic susceptibility tensor,
%                 cubic Angstroms
%
%         err   - least squares error
%
% Note: Gauss units are used for hyperfine couplings because they do
%       not depend on the electron g-tensor.
%
% ilya.kuprov@weizmann.ac.uk
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=pms2chi.m>

function [chi,err]=pms2chi(hfcs,shifts,isotopes,nel)

% Check consistency
grumble(hfcs,shifts,isotopes,nel);

% Set minimizer parameters
options=optimoptions('fminunc','Algorithm','quasi-newton','Display','iter',...
                     'GradObj','off','HessUpdate','bfgs','MaxIter',100,...
                     'MaxFunEvals',Inf,'UseParallel',true);
                 
% Set the initial guess
guess=[0 0 0 0 0 0];

% Run the optimisation
[chi,err]=fminunc(@(x)lsq_err([x(1) x(2) x(3); 
                               x(2) x(4) x(5); 
                               x(3) x(5) x(6)],...
          hfcs,shifts,isotopes,nel),guess,options);

% Form the answer
chi=[chi(1) chi(2) chi(3);
     chi(2) chi(4) chi(5);
     chi(3) chi(5) chi(6)];

end

% Least squares error function
function err=lsq_err(chi,hfcs,shifts,isotopes,nel)

% Get the error going
err=0;

% Build the list squares error
for n=1:numel(hfcs)
    shift_expt=shifts(n);
    shift_theo=hfc2pms(hfcs{n},chi,isotopes{n},nel);
    err=err+(shift_expt-shift_theo)^2;
end

end

% Consistency enforcement
function grumble(hfcs,shifts,isotopes,nel)
if ~iscell(hfcs)
    error('hfcs must be a cell array of matrices.');
end
for n=1:numel(hfcs)
    if (~isnumeric(hfcs{n}))||(~isreal(hfcs{n}))||...
       (~issymmetric(hfcs{n}))||(any(size(hfcs{n})~=3))
        error('all elements of hfcs cell array must be real symmetric 3x3 matrices.');
    end
end
if numel(hfcs)~=numel(shifts)
    error('hfcs and shifts arrays must have the same number of elements.');
end
if (~isnumeric(shifts))||(~isreal(shifts))
    error('shifts must be a vector of real numbers.');
end
if ~iscell(isotopes)
    error('isotopes must be a cell array of character strings.');
end
if numel(hfcs)~=numel(isotopes)
    error('hfcs and isotopes arrays must have the same number of elements.');
end
for n=1:numel(isotopes)
    if ~ischar(isotopes{n})
        error('all elements of the isotopes cell array must be character strings.');
    end
end
if (~isnumeric(nel))||(~isreal(nel))||...
   (numel(nel)~=1)||(mod(nel,1)~=0)||(nel<1)
    error('nel must be a non-negative real integer.');
end
end

% The problem with socialism is that you eventually 
% run out of other peoples' money.
%
% Margaret Thatcher

