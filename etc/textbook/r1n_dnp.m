% A simple model for nuclear longitudnal relaxation
% rate dies to the presence of an unpaired electron
% in cryogenic DNP settings. 
%
%         [literature reference goes here]
%
% Syntax:
%
%       R1n=r1n_dnp(B0,T,g,T1e,T1n_bulk,r,bet)
% 
% Parameters:
%
%     B0       - main magnet field, Tesla
%
%     T        - absolute temperature, Kelvin
%
%     g        - electron g-factor, Bohr 
%                magneton units
%
%     T1e      - electron longitudinal rela-
%                xation time, seconds
%
%     T1n_bulk - nuclear longitudinal relaxa-
%                tion time far away from the 
%                electron
%
%     r        - electron-nuclear distance,
%                Angstrom
%
%     bet      - angle between the magnet field
%                and the electron-nuclear direc-
%                tion, radians
%
% Outputs:
%
%     R1n      - nuclear relaxation rate, Hz
%
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function R1n=r1n_dnp(B0,T,g,T1e,T1n_bulk,r,bet)

% Fundamental constants
mu0=4*pi*1e-7;    % Vacuum permeability
muB=9.274010e-24; % Bohr magneton
kB =1.380649e-23; % Boltzmann constant

% Nuclear longitudinal relaxation rate, Hz
sech_sq=sech(g*muB*B0/(2*kB*T))^2;
geom_dd=(1-3*cos(bet)^2)/(r/1e10)^3;
R1n=(((mu0/(4*pi))*(g*muB/B0)*geom_dd)^2)*sech_sq/T1e+1/T1n_bulk;

end

% Ludwig Boltzmann, who spent much of his life studying statistical
% mechanics, died in 1906, by his own hand. Paul Ehrenfest, carrying
% on the work, died similarly in 1933. Now it is our turn to study 
% statistical mechanics. Perhaps it will be wise to approach the sub-
% ject cautiously.
%
% David Goodstein, "States of Matter"

