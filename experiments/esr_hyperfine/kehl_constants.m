% defines the necessary constants
%
% no input parameters
%
% output parameters:
% constants: the Map with the constants
% February 2024 A. Kehl  (akehl@gwdg.de)

function constants=kehl_constants()
    constants=containers.Map;

    % MU_B*ge/H (convert A from MHz to G)
    constants("CONST1")=2.81;

    % Boltzmann Constant    [Hz/K]
    constants("K_B")=2.0836618E10;

    % Planck Constant       [J s]
    constants("H")=6.62607E-34;

    % Bohr Magneton         [J/T]
    constants("MU_B")=9.27401E-24;

    % Nuclear Magneton      [J/T]
    constants("MU_N")=5.050783699E-27;

    % gyromag. ratio e- /2pi    [Hz/T]
    constants("GE")=28024.95266E6;

    % gyromag. ratio 1H /2pi    [Hz/T]
    constants("GN_1H")=42.576E6;

    % gyromag. ratio 2D /2pi    [Hz/T]
    constants("GN_2D")=6.536E6;

    % gyromag. ratio 19F /2pi   [Hz/T]
    constants("GN_19F")=40.078E6;

    % gyromag. ratio 14N /2pi   [Hz/T]
    constants("GN_14N")=3.077E6;

    % gyromag. ratio 17O /2pi   [Hz/T]
    constants("GN_17O")=-5.772E6;
end

