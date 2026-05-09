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
%   varargout        - ENDOR amplitude, broadened spectrum, frequency axis, optional field axis, and Larmor-frequency data.
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
    parameters=pulse_sequence(spin_system,parameters,'parameters',[],[]);

    % Select EPR orientations using sequence-independent context data
    parameters.paramsEPR=kehl_prep_epr(spin_system,parameters);
    if parameters.freqDomain==false
        parameters.epr=kehl_ori_field(spin_system,parameters);
    else
        parameters.epr=kehl_ori_freq(spin_system,parameters);
    end

    % Build Spinach relaxation superoperator when requested
    H=[];
    if parameters.Relax==true
        R=kehl_context_relaxation(spin_system,parameters);
    else
        R=[];
    end
    K=[];

    % Run the requested pulse sequence using the Spinach experiment signature
    endor_amp=pulse_sequence(spin_system,parameters,H,R,K);
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

function R=kehl_context_relaxation(spin_system,parameters)

    % Build a spherical-tensor Liouville copy for Spinach T1/T2 relaxation
    rlx_spin_system=spin_system;
    rlx_spin_system.bas.formalism='sphten-liouv';
    rlx_spin_system.bas.approximation='none';
    rlx_spin_system.rlx.theories={'t1_t2'};
    rlx_spin_system.rlx.r1_rates=kehl_relaxation_rates(spin_system,parameters,'r1');
    rlx_spin_system.rlx.r2_rates=kehl_relaxation_rates(spin_system,parameters,'r2');
    rlx_spin_system.rlx.keep='diagonal';
    rlx_spin_system.rlx.equilibrium='zero';
    rlx_spin_system.rlx.dfs='ignore';
    rlx_spin_system=basis(rlx_spin_system,rlx_spin_system.bas);

    % Calculate Spinach relaxation in spherical-tensor Liouville space
    R_sphten=relaxation(rlx_spin_system);

    % Transform the superoperator into the Kehl Zeeman-Liouville basis
    P=sphten2zeeman(rlx_spin_system);
    if size(P,1)~=size(P,2)
        error('Kehl ENDOR relaxation requires a complete spherical-tensor basis.');
    end
    R=sparse(P*R_sphten/P);

end

function rates=kehl_relaxation_rates(spin_system,parameters,kind)

    % Preallocate one relaxation rate per spin
    rates=cell(1,spin_system.comp.nspins);

    % Convert Kehl relaxation times to Spinach rates
    for n=1:spin_system.comp.nspins
        if kehl_is_electron(spin_system.comp.isotopes{n})
            if strcmp(kind,'r1')
                rates{n}=kehl_rate(parameters.T1e);
            else
                rates{n}=kehl_rate(parameters.T2e);
            end
        else
            if strcmp(kind,'r1')
                rates{n}=kehl_rate(parameters.T1n);
            else
                rates{n}=kehl_rate(parameters.T2n);
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
    if (~isa(sequence,'function_handle'))&&(~ischar(sequence))&&(~isstring(sequence))
        error('sequence must be a function handle, character string, or string scalar.');
    end
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
    if ~isfield(parameters,'endor_spins')
        error('parameters.endor_spins must list ENDOR nuclei by Spinach spin index.');
    end
    if ~isfield(parameters,'mw_freq_ghz')
        error('parameters.mw_freq_ghz must be specified.');
    end
    if ~isfield(parameters,'static_field_g')
        error('parameters.static_field_g must be specified.');
    end
    if ~isfield(parameters,'field_step_g')
        error('parameters.field_step_g must be specified.');
    end
    if ~isfield(parameters,'endor_res_mhz')
        error('parameters.endor_res_mhz must be specified.');
    end
    if ~isfield(parameters,'endor_range_mhz')
        error('parameters.endor_range_mhz must be specified.');
    end
    if ~isfield(parameters,'pulse_times_ns')
        error('parameters.pulse_times_ns must be specified.');
    end
    if parameters.freqDomain==true
        if ~isfield(parameters,'epr_freq_min_ghz')
            error('parameters.epr_freq_min_ghz must be specified for frequency-domain EPR.');
        end
        if ~isfield(parameters,'epr_freq_range_ghz')
            error('parameters.epr_freq_range_ghz must be specified for frequency-domain EPR.');
        end
        if ~isfield(parameters,'epr_freq_step_ghz')
            error('parameters.epr_freq_step_ghz must be specified for frequency-domain EPR.');
        end
    end
    if (~ischar(assumptions))&&(~isstring(assumptions))
        error('assumptions must be a character string or a string scalar.');
    end
end

