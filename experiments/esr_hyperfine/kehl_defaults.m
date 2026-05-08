% Default simulation parameters for the Kehl ENDOR context. Syntax:
%
%      parameters=kehl_defaults(parameters)
%
% Parameters:
%
%   parameters       - Kehl ENDOR context parameter structure.
%
% Outputs:
%
%   parameters       - parameter structure with missing defaults populated.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_defaults.m>

function parameters=kehl_defaults(parameters)
    if nargin<1
        parameters=struct();
    end

    % Check consistency
    grumble(parameters);

    % Is the powder pattern used?
    if ~isfield(parameters,'powder')
        parameters.powder=false;
    end

    % Is the EPR spectrum calculated in the frequency, or field domain?
    if ~isfield(parameters,'freqDomain')
        parameters.freqDomain=false;
    end

    % Is the B term included?
    if ~isfield(parameters,'Bterm')
        parameters.Bterm=false;
    end

    % How many steps are calculated for RF-pulse B-term propagation?
    if ~isfield(parameters,'N_stepRF')
        parameters.N_stepRF=100;
    end

    % Is relaxation included?
    if ~isfield(parameters,'Relax')
        parameters.Relax=false;
    end

    % Are temperature effects used?
    if ~isfield(parameters,'temp_eff')
        parameters.temp_eff=false;
    end

    % Use temperature effects if relaxation is used
    if parameters.Relax==true
        if parameters.temp_eff==false
            warning('The temperature needs to be considered if relaxation is used.');
        end
        parameters.temp_eff=true;
    end

    % Set the temperature
    if ~isfield(parameters,'T')
        parameters.T=10;
    end

    % Choose the EPR line for single-crystal calculation
    if ~isfield(parameters,'sel_I')
        parameters.sel_I=1;
    end

    % Choose the CP condition for CP ENDOR
    if ~isfield(parameters,'sel_CP')
        parameters.sel_CP=1;
    end

    % Set the powder grid resolution
    if ~isfield(parameters,'Nang')
        parameters.Nang=20;
    end

    % Set the excitation bandwidth factor
    if ~isfield(parameters,'nwidth')
        parameters.nwidth=20;
    end

    % Is Lorentzian line broadening used?
    if ~isfield(parameters,'Lorentzian')
        parameters.Lorentzian=false;
    end

    % Lorentzian line broadening, Hz
    if ~isfield(parameters,'lw_L')
        parameters.lw_L=20*1e3;
    end

    % Is Gaussian line broadening used?
    if ~isfield(parameters,'Gaussian')
        parameters.Gaussian=false;
    end

    % Gaussian line broadening, Hz
    if ~isfield(parameters,'lw_G')
        parameters.lw_G=20*1e3;
    end
end
function grumble(parameters)
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
end

