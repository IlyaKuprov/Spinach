% Converts hyperfine couplings from Gauss to MHz (linear
% frequency). The Gauss specification may be defined as 
% "the magnetic field at which the electron frequency is
% equal to the frequency provided". Syntax:
%
%              hfc_mhz=gauss2mhz(hfc_gauss,g)
%
% Arrays of any dimensions are supported. Parameters:
%
%   hfc_gauss  - an array of values in Gauss
%
%   g          - electron g-factor; if this parameter
%                is skipped, free electron g-factor is
%                used for conversion
%
% Outputs:
%
%   hfc_mhz    - an array of values in MHz
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=gauss2mhz.m>

function hfc_mhz=gauss2mhz(hfc_gauss,g)

% Set the defaults
if ~exist('g','var')
    disp('gauss2mhz: using free electron g-factor.');
    g=2.0023193043622;
end

% Check consistency
grumble(hfc_gauss,g);

% Do the conversion
muB=9.274009994*10^-24;
hbar=1.054571628e-34;
C=1e-10*g*muB/(hbar*2*pi);
hfc_mhz=C*hfc_gauss;

end

% Consistency enforcement
function grumble(hfc_gauss,g)
if (~isnumeric(hfc_gauss))||(~isreal(hfc_gauss))
    error('hfc_gauss must be an array of real numbers.');
end
if (~isnumeric(g))||(~isreal(g))||(numel(g)~=1)
    error('g must be a real number.');
end
end

% A celibate clergy is an especially good idea, because 
% it tends to suppress any hereditary propensity toward
% fanaticism.
%
% Carl Sagan

