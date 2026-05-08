%KEHL_CONTEXT_FIELDS Generic field and sweep parameters for Kehl ENDOR context.
%
%   Spinach architecture migration May 2026 Talos

function parameters=kehl_context_fields(parameters)

% Check consistency
grumble(parameters);


    % Append generic field and frequency data
    parameters.mw_freq_hz=parameters.mw_freq_ghz*1e9;
    parameters.field_t=parameters.static_field_g*1e-4;
    parameters.field_step_t=parameters.field_step_g*1e-4;
    parameters.endor_res_hz=parameters.endor_res_mhz*1e6;
    parameters.endor_range_hz=parameters.endor_range_mhz*1e6;
    parameters.pulse_times_s=parameters.pulse_times_ns*1e-9;

    % Append frequency-domain EPR sweep data when present
    if isfield(parameters,'epr_freq_min_ghz')
        parameters.epr_freq_min_hz=parameters.epr_freq_min_ghz*1e9;
        parameters.epr_freq_range_hz=parameters.epr_freq_range_ghz*1e9;
        parameters.epr_freq_step_hz=parameters.epr_freq_step_ghz*1e9;
    end

    % Append optional direct RF and CP sweep starts
    if isfield(parameters,'rf_start_mhz')
        parameters.rf_start_hz=parameters.rf_start_mhz*1e6;
    end
    if isfield(parameters,'cp_start_mhz')
        parameters.cp_start_hz=parameters.cp_start_mhz*1e6;
    end
    if isfield(parameters,'cp_range_mhz')
        parameters.cp_range_hz=parameters.cp_range_mhz*1e6;
    end

    % Default shaped-pulse multiplicity flag
    if isfield(parameters,'pulse_file')&&~isfield(parameters,'multipulses')
        parameters.multipulses=false;
    end
end

% Consistency enforcement
function grumble(parameters)
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
end
