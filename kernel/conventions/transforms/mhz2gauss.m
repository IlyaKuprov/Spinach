% Converts hyperfine couplings from MHz (linear frequency)
% to Gauss. The Gauss specification may be defined as 
% "the magnetic field at which the electron frequency is
% equal to the frequency provided". Syntax:
%
%               hfc_gauss=mhz2gauss(hfc_mhz,g)
%
% Arrays of any dimensions are supported. Parameters:
%
%   hfc_mhz    - an array of values in MHz
%
%   g          - electron g-factor; if this parameter
%                is skipped, free electron g-factor is
%                used for conversion
%
% Outputs:
%
%   hfc_gauss  - an array of values in Gauss
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=mhz2gauss.m>

function hfc_gauss=mhz2gauss(hfc_mhz,g)

% Set the defaults
if ~exist('g','var')
    disp('mhz2gauss: using free electron g-factor.');
    g=2.0023193043622;
end

% Check consistency
grumble(hfc_mhz,g);

% Do the conversion
muB=9.274009994*10^-24;
hbar=1.054571628e-34;
C=1e-10*g*muB/(hbar*2*pi);
hfc_gauss=hfc_mhz/C;
    
end

% Consistency enforcement
function grumble(hfc_mhz,g)
if (~isnumeric(hfc_mhz))||(~isreal(hfc_mhz))
    error('hfc_gauss must be an array of real numbers.');
end
if (~isnumeric(g))||(~isreal(g))||(numel(g)~=1)
    error('g must be a real number.');
end
end

% "Trillian had come to suspect that the main reason [Zaphood] had 
% had such a wild and successful life was that he never really un-
% derstood the significance of anything he did."
%
% Douglas Adams, "The Hitchhiker's Guide to the Galaxy"

