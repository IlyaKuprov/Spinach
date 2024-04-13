% Converts hyperfine couplings from milliTesla to Hz (linear
% frequency). The milliTesla specification may be defined as
% "the magnetic field at which the electron frequency is equ-
% al to the frequency provided". Syntax:
%
%                  hfc_hz=mt2hz(hfc_mt,g)
%
% Arrays of any dimension are supported. Parameters:
%
%   hfc_mt     - an array of values in mT
%
%   g          - electron g-factor; if this parameter
%                is skipped, free electron g-factor is
%                used for conversion
%
% Outputs:
%
%   hfc_hz     - an array of values in Hz
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=mt2hz.m>

function hfc_hz=mt2hz(hfc_mt,g)

% Set the defaults
if ~exist('g','var')
    disp('mt2hz: using free electron g-factor.');
    g=2.0023193043622;
end

% Check consistency
grumble(hfc_mt,g);

% Do the conversion
muB=9.274009994*10^-24;
hbar=1.054571628e-34;
C=1e-3*g*muB/(hbar*2*pi);
hfc_hz=C*hfc_mt;
   
end

% Consistency enforcement
function grumble(hfc_mt,g)
if (~isnumeric(hfc_mt))||(~isreal(hfc_mt))
    error('hfc_mt must be an array of real numbers.');
end
if (~isnumeric(g))||(~isreal(g))||(numel(g)~=1)
    error('g must be a real number.');
end
end

% You know that I write slowly. This is chiefly because I 
% am never satisfied until I have said as much as possible
% in a few words, and writing briefly takes far more time
% than writing at length.
%
% Carl Friedrich Gauss

