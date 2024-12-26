% Converts full width at half-maximum (FWHM) of an NMR 
% signal into an approximation of the R2 rate. Syntax:
%
%                 r2rate=fwhm2rlx(fwhm)
%
% Parameters:
%
%     fwhm  - full width at half-maximum, Hz
%
% Outputs:
%
%   r2rate  - approximate R2 relaxation rate, Hz
%
% Note: FWHM is not a reliable measure of the transverse
%       relaxation rate. The value obtained from this 
%       function should be treated as an upper bound.
%
% Note: Lorentzian line shape is assumed.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=fwhm2rlx.m>

function r2rate=fwhm2rlx(fwhm)

% Check consistency
grumble(fwhm);

% Run the conversion
r2rate=pi*fwhm;

end

% Consistency enforcement
function grumble(fwhm)
if (~isnumeric(fwhm))||(~isreal(fwhm))||(any(fwhm(:)<=0))
    error('fwhm must be an array of positive real numbers.');
end
end

% Any fool can use a computer. Many do.
%
% Ted Nelson

