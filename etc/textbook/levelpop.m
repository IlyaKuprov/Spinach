% Equilibrium populations of the energy levels of a user-specified spin at
% the user-specified temperature. Energies are reported as fractions of kT
% at the temperature specified. Syntax:
%
%                [E,P,dP]=levelpop(isotope,field,temperature)
%
% Parameters:
%
%       isotope      - character string specifying the isotope.
%                      e.g. '1H', '13C', 'E', etc.
%
%       field        - primary magnet field in Tesla
%
%       temperature  - spin temperature, Kelvin
%
% Outputs:
%
%       E            - vector of level energies, frac-
%                      tions of kT at the temperature
%                      specified
%
%       P            - vector of level populations
%
%       dP           - vector of population differences
%                      for adjacent levels
%
% Notes: the function is sensitive to the sign of the magnetogyric
%        ratio - negative for electrons, positive for protons, etc.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=levelpop.m>

function [E,P,dP]=levelpop(isotope,field,temperature)

% Check consistency
grumble(isotope,field,temperature);

% Fundamental constants
h_bar=6.62607015e-34/(2*pi); % J*s, exact number
k_bol=1.380649e-23;          % J/K, exact number

% Get the spin data
[mg_ratio,multipl]=spin(isotope);

% Get Pauli matrices
S=pauli(multipl);

% Build Zeeman Hamiltonian
H=-mg_ratio*field*S.z; H=full(H);

% Get the energies
E=h_bar*diag(H)/(k_bol*temperature);

% Get the populations
P=exp(-h_bar*diag(H)/(k_bol*temperature));

% Normalise populations
P=P/sum(P);

% Population differences
dP=-diff(P);

end

% Consistency enforcement
function grumble(isotope,field,temperature)
if ~ischar(isotope)
    error('isotope must be a character string.');
end
if (~isnumeric(field)||(~isreal(field)))||(~isscalar(field))
    error('field must be a real scalar.');
end
if (~isnumeric(temperature)||(~isreal(temperature)))||...
   (~isscalar(temperature))||(temperature==0)
    error('temperature must be a non-zero real scalar.');
end
end

% We are not hypocrites in our sleep. 
%
% William Hazlitt

