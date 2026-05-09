% Kehl ENDOR context wrapper for Spinach-style pulse sequences. Syntax:
%
%      varargout=endor_kehl_context(spin_system,sequence,parameters,assumptions)
%
% Parameters:
%
%   spin_system      - Spinach spin system structure.
%   sequence         - pulse-sequence function handle, or Kehl sequence name.
%   parameters       - Kehl ENDOR context parameter structure.
%   assumptions      - Spinach assumption set, defaults to labframe.
%
% Outputs:
%
%   varargout        - ENDOR amplitude, broadened spectrum, frequency/time axis, optional field axis, and Larmor angular-frequency data.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=endor_kehl_context.m>

function varargout=endor_kehl_context(spin_system,sequence,parameters,assumptions)
    if nargin<4
        assumptions='labframe';
    end

    % Set default simulation parameters
    parameters=kehl_defaults(parameters);

    % Check consistency
    grumble(spin_system,sequence,parameters,assumptions);

    % Normalise the pulse-sequence selector
    pulse_sequence=sequence_handle(sequence);

    % Store context assumptions
    parameters.assumptions=assumptions;

    % Build sequence-agnostic Kehl context data from Spinach parameters
    parameters.constants=kehl_context_constants(spin_system);
    parameters=kehl_context_fields(parameters);
    parameters=kehl_context_spin_data(spin_system,parameters);

    if nargin>=4&&~isempty(assumptions)
        spin_system=assume(spin_system,assumptions);
    end

    % Let the pulse sequence append its own derived parameters
    parameters=pulse_sequence(spin_system,parameters,'parameters',[]);

    % Select EPR orientations using sequence-independent context data
    parameters.paramsEPR=kehl_prep_epr(parameters);
    if parameters.freqDomain==false
        parameters.epr=kehl_ori_field(spin_system,parameters);
    else
        parameters.epr=kehl_ori_freq(spin_system,parameters);
    end

    % Build Spinach relaxation superoperator when requested
    H=[];
    if parameters.Relax==true
        [spin_system,R]=kehl_context_relaxation(spin_system,parameters);
    else
        R=mprealloc(spin_system,1);
    end

    % Run the requested pulse sequence using the Spinach experiment signature
    endor_amp=pulse_sequence(spin_system,parameters,H,R);
    endor_amp_conv=kehl_line_broaden(endor_amp,parameters);

    x_coords=parameters.paramsENDOR('x_coords');
    x_coords=x_coords(:,1)';
    v_L=parameters.paramsENDOR('v_L');

    if isfield(parameters,'y_coords')
        varargout={endor_amp,endor_amp_conv,x_coords,parameters.y_coords,v_L};
    else
        varargout={endor_amp,endor_amp_conv,x_coords,v_L};
    end
end

function pulse_sequence=sequence_handle(sequence)
    if isa(sequence,'function_handle')
        pulse_sequence=sequence;
        return
    end
    name=char(sequence);
    name=lower(name);
    if strcmp(name,'timedomain')
        name='time';
    end
    if ~startsWith(name,'endor_kehl_')
        name=['endor_kehl_' name];
    end
    pulse_sequence=str2func(name);
end

function [spin_system,R]=kehl_context_relaxation(spin_system,parameters)

    % Configure Spinach Lindblad T1/T2 relaxation on the existing system
    spin_system.rlx.theories={'lindblad'};
    spin_system.rlx.lind_r1_rates=kehl_relaxation_rates(...
        spin_system,parameters,'r1');
    spin_system.rlx.lind_r2_rates=kehl_relaxation_rates(...
        spin_system,parameters,'r2');
    spin_system.rlx.keep='diagonal';
    spin_system.rlx.equilibrium='zero';
    spin_system.rlx.dfs='ignore';

    % Calculate Spinach relaxation in the existing Zeeman-Liouville basis
    R=relaxation(spin_system);

end

function rates=kehl_relaxation_rates(spin_system,parameters,kind)

    % Preallocate one relaxation rate per spin
    rates=zeros(1,spin_system.comp.nspins);

    % Convert Kehl relaxation times to Spinach rates
    for n=1:spin_system.comp.nspins
        if kehl_is_electron(spin_system.comp.isotopes{n})
            if strcmp(kind,'r1')
                rates(n)=kehl_rate(parameters.T1e);
            else
                rates(n)=kehl_rate(parameters.T2e);
            end
        else
            if strcmp(kind,'r1')
                rates(n)=kehl_rate(parameters.T1n);
            else
                rates(n)=kehl_rate(parameters.T2n);
            end
        end
    end

end

function rate=kehl_rate(relax_time)
    if isinf(relax_time)
        rate=0;
    else
        rate=1/relax_time;
    end
end

function grumble(spin_system,sequence,parameters,assumptions)
    if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
        error('spin_system must be a Spinach spin system structure.');
    end
    if ~strcmp(spin_system.bas.formalism,'zeeman-liouv')
        error('spin_system must use zeeman-liouv formalism.');
    end
    if (~isa(sequence,'function_handle'))&&(~ischar(sequence))&&(~isstring(sequence))
        error('sequence must be a function handle, character string, or string scalar.');
    end
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
    if ~isfield(parameters,'endor_spins')
        error('parameters.endor_spins must list ENDOR nuclei by Spinach spin index.');
    end
    if ~isfield(parameters,'mw_freq')
        error('parameters.mw_freq must be specified.');
    end
    if ~isfield(parameters,'static_field')
        error('parameters.static_field must be specified.');
    end
    if ~isfield(parameters,'field_step')
        error('parameters.field_step must be specified.');
    end
    if isa(sequence,'function_handle')
        sequence_name=func2str(sequence);
    else
        sequence_name=char(sequence);
    end
    if strcmp(sequence_name,'endor_kehl_time')||strcmp(sequence_name,'time')||...
            strcmp(sequence_name,'timedomain')
        if ~isfield(parameters,'time_res')
            error('parameters.time_res must be specified for time-domain ENDOR.');
        end
        if ~isfield(parameters,'time_range')
            error('parameters.time_range must be specified for time-domain ENDOR.');
        end
    else
        if ~isfield(parameters,'endor_res')
            error('parameters.endor_res must be specified.');
        end
        if ~isfield(parameters,'endor_range')
            error('parameters.endor_range must be specified.');
        end
    end
    if ~isfield(parameters,'pulse_times')
        error('parameters.pulse_times must be specified.');
    end
    if parameters.freqDomain==true
        if ~isfield(parameters,'epr_freq_min')
            error('parameters.epr_freq_min must be specified for frequency-domain EPR.');
        end
        if ~isfield(parameters,'epr_freq_range')
            error('parameters.epr_freq_range must be specified for frequency-domain EPR.');
        end
        if ~isfield(parameters,'epr_freq_step')
            error('parameters.epr_freq_step must be specified for frequency-domain EPR.');
        end
    end
    if (~ischar(assumptions))&&(~isstring(assumptions))
        error('assumptions must be a character string or a string scalar.');
    end
end

