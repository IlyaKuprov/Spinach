%ENDOR_KEHL_CONTEXT Spinach-style ENDOR context for the Kehl pulse sequences.
%
%   [ENDOR,ENDOR_LB,X,V]=ENDOR_KEHL_CONTEXT(SPIN_SYSTEM,PULSE_SEQUENCE,...
%   PARAMETERS,ASSUMPTIONS) builds Kehl experimental parameters from the
%   physical fields in PARAMETERS, prepares the Spinach-backed spin
%   operators, EPR orientation selection, ENDOR sweep axes, and line
%   broadening, then calls a pulse sequence function with Spinach's
%   standard (spin_system,parameters,H,R,K) experiment signature.
%
%   For CP ENDOR, the return signature is
%   [ENDOR,ENDOR_LB,X,Y,V].


function varargout=endor_kehl_context(spin_system,pulse_sequence,parameters,assumptions)
if nargin<4
    assumptions='labframe';
end

% Set default simulation parameters
parameters=kehl_defaults(parameters);

% Check consistency
grumble(spin_system,pulse_sequence,parameters,assumptions);

% Build legacy Kehl data structures from Spinach parameters
constants=kehl_constants();
expt=context_experiment(parameters,constants,pulse_sequence);
spinSys=context_spin_system(spin_system,parameters,constants);
spinOps=kehl_spin_ops(spin_system,parameters.endor_spins,spinSys('N_SpinSys'));

if nargin>=4 && ~isempty(assumptions)
    spin_system=assume(spin_system,assumptions);
end

if isfield(parameters,'time_domain') && parameters.time_domain
    t=expt('t');
    expt('t_start')=t(6);
    paramsENDOR=kehl_prep_time(constants,spinSys,expt);
else
    paramsENDOR=kehl_prep_endor(constants,spinSys,expt);
end

paramsEPR=kehl_prep_epr(spinSys,expt);
if parameters.freqDomain==false
    epr=kehl_ori_field(constants,spinSys,spinOps,paramsEPR,paramsENDOR,parameters,expt);
else
    epr=kehl_ori_freq(constants,spinSys,spinOps,paramsEPR,paramsENDOR,parameters,expt);
end

localpar=parameters;
localpar.constants=constants;
localpar.spinSys=spinSys;
localpar.spinOps=spinOps;
localpar.expt=expt;
localpar.paramsENDOR=paramsENDOR;
localpar.paramsEPR=paramsEPR;
localpar.epr=epr;
localpar.assumptions=assumptions;

is_cp=isfield(parameters,'cp') && parameters.cp;
if is_cp
    [y_coords,localpar]=cp_axis(localpar);
end

H=[];
R=[];
K=[];
endor_amp=pulse_sequence(spin_system,localpar,H,R,K);
endor_amp_conv=kehl_line_broaden(endor_amp,parameters,localpar.expt('range_EN'));

x_coords=paramsENDOR('x_coords');
x_coords=x_coords(:,1)';
v_L=paramsENDOR('v_L');

if is_cp
    varargout={endor_amp,endor_amp_conv,x_coords,y_coords,v_L};
else
    varargout={endor_amp,endor_amp_conv,x_coords,v_L};
end
end

function spinSys=context_spin_system(spin_system,parameters,constants)
    isotopes=spin_system.comp.isotopes;
    electron_idx=find(cellfun(@is_electron,isotopes),1);
    if isempty(electron_idx)
        error('spin_system must contain an electron spin.');
    end

    endor_spins=parameters.endor_spins;
    epr_spins=infer_epr_spins(spin_system,electron_idx,endor_spins);
    n_spin_systems=1;

    [~,electron_mult]=spin(isotopes{electron_idx});
    n_endor=numel(endor_spins);

    spinSys=containers.Map;
    spinSys('S')=(electron_mult-1)/2;
    spinSys('g')=electron_g_matrix(spin_system,electron_idx);
    spinSys('g_iso')=trace(spinSys('g'))/3;
    spinSys('Ni_ENDOR')=n_endor;
    spinSys('Nuclei')=isotopes(endor_spins);
    spinSys('N_SpinSys')=n_spin_systems;
    spinSys('I')=spin_numbers(isotopes(endor_spins));

    A=zeros(3*n_endor,3);
    Q=zeros(3*n_endor,3);
    Q_used=false;
    for n=1:n_endor
        spin_idx=endor_spins(n);
        A(3*n-2:3*n,:)=coupling_matrix(spin_system,electron_idx,spin_idx);
        Q_block=coupling_matrix(spin_system,spin_idx,spin_idx);
        Q(3*n-2:3*n,:)=Q_block;
        Q_used=Q_used||any(Q_block(:));
    end
    spinSys('A')=A;
    spinSys('Q')=Q;
    spinSys('Q_used')=Q_used;

    CS=zeros(3*n_endor,3);
    CS_used=false;
    for n=1:n_endor
        spin_idx=endor_spins(n);
        CS_block=nuclear_cs_matrix(spin_system,spin_idx);
        CS(3*n-2:3*n,:)=CS_block;
        CS_used=CS_used||any(CS_block(:));
    end
    if CS_used
        spinSys('CS')=CS;
    end
    spinSys('CS_used')=CS_used;

    pairs=infer_dipolar_pairs(spin_system,endor_spins);
    if ~isempty(pairs)
        D=zeros(3*size(pairs,1),3);
        for n=1:size(pairs,1)
            D(3*n-2:3*n,:)=coupling_matrix(spin_system,pairs(n,1),pairs(n,2));
        end
        spinSys('D')=D;
        spinSys('D_used')=true;
    else
        spinSys('D')=zeros(0,3);
        spinSys('D_used')=false;
    end

    spinSys('EPR_Nucs_used')=~isempty(epr_spins);
    if ~isempty(epr_spins)
        n_epr=numel(epr_spins);
        spinSys('EPR_Nuclei')=isotopes(epr_spins);
        spinSys('Ni_EPR')=n_epr;
        spinSys('I_EPR')=spin_numbers(isotopes(epr_spins));
        g_N_EPR=zeros(n_epr,1);
        A_EPR=zeros(3*n_epr,3);
        Q_EPR=zeros(3*n_epr,3);
        EPR_Q_used=false;
        for n=1:n_epr
            spin_idx=epr_spins(n);
            g_N_EPR(n)=kehl_nuc_gamma(constants,isotopes{spin_idx});
            A_EPR(3*n-2:3*n,:)=coupling_matrix(spin_system,electron_idx,spin_idx);
            Q_block=coupling_matrix(spin_system,spin_idx,spin_idx);
            Q_EPR(3*n-2:3*n,:)=Q_block;
            EPR_Q_used=EPR_Q_used||any(Q_block(:));
        end
        spinSys('g_N_EPR')=g_N_EPR;
        spinSys('A_EPR')=A_EPR;
        spinSys('Q_EPR')=Q_EPR;
        spinSys('EPR_Q_used')=EPR_Q_used;
    else
        spinSys('EPR_Q_used')=false;
    end
end

function [y_coords,localpar]=cp_axis(localpar)
    expt=localpar.expt;
    if localpar.powder==false
        expt('Npts_CP')=1;
        expt('range_CP')=0;
        y_coords=1;
    else
        if expt('Npts_CP')>1
            step_CP=expt('range_CP')/(expt('Npts_CP')-1);
            y_coords=zeros(expt('Npts_CP'));
            for n=1:expt('Npts_CP')
                y_coords(n)=expt('start_CP')+(n-1)*step_CP;
            end
            y_coords=y_coords(:,1)';
        else
            y_coords=expt('start_CP');
            y_coords=y_coords(:,1)';
        end
        localpar.paramsENDOR('y_coords')=y_coords;
    end
    localpar.expt=expt;
end

function expt=context_experiment(parameters,constants,pulse_sequence)
    % Create the primary experiment parameter map
    expt=kehl_exp_create(parameters.mw_freq_ghz,parameters.static_field_g,...
                         parameters.field_step_g,parameters.endor_res_mhz,...
                         parameters.endor_range_mhz,parameters.pulse_times_ns);

    % Set microwave and radiofrequency nutation fields
    [rf_auto,rf_nutations]=rf_field_spec(parameters);
    field_fun=context_field_fun(parameters,pulse_sequence);
    expt=field_fun(constants,expt,rf_auto,rf_nutations);

    % Set the frequency-domain EPR sweep
    if isfield(parameters,'epr_freq_min_ghz')
        expt=kehl_exp_freq(expt,parameters.epr_freq_min_ghz,...
                           parameters.epr_freq_range_ghz,...
                           parameters.epr_freq_step_ghz);
    end

    % Set the ENDOR radiofrequency pulse flip angle
    if isfield(parameters,'rf_flip_angle_deg')
        expt=kehl_exp_angle(expt,parameters.rf_flip_angle_deg);
    end

    % Set the shaped pulse profile file
    if isfield(parameters,'pulse_file')
        multipulses=false;
        if isfield(parameters,'multipulses')
            multipulses=parameters.multipulses;
        end
        expt=kehl_exp_pulse(expt,parameters.pulse_file,multipulses);
    end

    % Set the CP sweep dimension when present
    if isfield(parameters,'cp_start_mhz')
        expt=kehl_cp_def(expt,parameters.cp_start_mhz,...
                         parameters.cp_npoints,parameters.cp_range_mhz);
    end
end

function [rf_auto,rf_nutations]=rf_field_spec(parameters)
    rf_nutations=[];
    if isfield(parameters,'rf_nutation_freqs')
        rf_nutations=parameters.rf_nutation_freqs;
    end
    if isfield(parameters,'rf_field_from_pulses')
        rf_auto=parameters.rf_field_from_pulses;
    else
        rf_auto=isempty(rf_nutations);
    end
end

function field_fun=context_field_fun(parameters,pulse_sequence)
    if isfield(parameters,'cp') && parameters.cp==true
        field_fun=@kehl_cp_fields;
        return
    end
    if isfield(parameters,'time_domain') && parameters.time_domain==true
        field_fun=@kehl_time_fields;
        return
    end
    seq_name=func2str(pulse_sequence);
    if contains(seq_name,'davies')
        field_fun=@kehl_davies_fields;
    elseif contains(seq_name,'spinlock')
        field_fun=@kehl_spinlock_fields;
    elseif contains(seq_name,'tensor')
        field_fun=@kehl_tensor_fields;
    elseif contains(seq_name,'cp')
        field_fun=@kehl_cp_fields;
    else
        field_fun=@kehl_mims_fields;
    end
end

function M=electron_g_matrix(spin_system,spin_idx)
    M=spin_system.inter.zeeman.matrix{spin_idx}*...
      spin_system.tols.freeg/spin_system.inter.basefrqs(spin_idx);
end

function M=nuclear_cs_matrix(spin_system,spin_idx)
    M=(spin_system.inter.zeeman.matrix{spin_idx}/...
       spin_system.inter.basefrqs(spin_idx)-eye(3,3))*1e12;
end

function M=coupling_matrix(spin_system,row,col)
    M=matrix_from_cell(spin_system.inter.coupling.matrix,row,col)/(2*pi);
end

function pairs=infer_dipolar_pairs(spin_system,endor_spins)
    pairs=zeros(numel(endor_spins)*(numel(endor_spins)-1)/2,2);
    pair_count=0;
    for n=1:numel(endor_spins)
        for k=n+1:numel(endor_spins)
            spin_a=endor_spins(n);
            spin_b=endor_spins(k);
            if norm(coupling_matrix(spin_system,spin_a,spin_b),'fro')>0
                pair_count=pair_count+1;
                pairs(pair_count,:)=[spin_a,spin_b];
            end
        end
    end
    pairs=pairs(1:pair_count,:);
end

function epr_spins=infer_epr_spins(spin_system,electron_idx,endor_spins)
    epr_spins=zeros(1,spin_system.comp.nspins);
    spin_count=0;
    for n=1:spin_system.comp.nspins
        if (n~=electron_idx) && (~ismember(n,endor_spins)) &&...
           (spin_system.comp.mults(n)>1)
            has_hfc=norm(coupling_matrix(spin_system,electron_idx,n),'fro')>0;
            has_nqi=norm(coupling_matrix(spin_system,n,n),'fro')>0;
            if has_hfc || has_nqi
                spin_count=spin_count+1;
                epr_spins(spin_count)=n;
            end
        end
    end
    epr_spins=epr_spins(1:spin_count);
end

function M=matrix_from_cell(cells,row,col)
    M=zeros(3,3);
    if isempty(cells)
        return
    end
    if size(cells,1)>=row && size(cells,2)>=col && ~isempty(cells{row,col})
        M=cells{row,col};
    elseif size(cells,1)>=col && size(cells,2)>=row && ~isempty(cells{col,row})
        M=cells{col,row};
    end
end

function qnums=spin_numbers(labels)
    qnums=zeros(numel(labels),1);
    for n=1:numel(labels)
        [~,multiplicity]=spin(labels{n});
        qnums(n)=(multiplicity-1)/2;
    end
end

function tf=is_electron(label)
    tf=strcmp(label,'E')||~isempty(regexp(label,'^E\d+$','once'));
end

% Consistency enforcement
function grumble(spin_system,pulse_sequence,parameters,assumptions)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isa(pulse_sequence,'function_handle')
    error('pulse_sequence must be a function handle.');
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
if isfield(parameters,'cp') && parameters.cp==true
    if ~isfield(parameters,'cp_start_mhz')
        error('parameters.cp_start_mhz must be specified for CP ENDOR.');
    end
    if ~isfield(parameters,'cp_npoints')
        error('parameters.cp_npoints must be specified for CP ENDOR.');
    end
    if ~isfield(parameters,'cp_range_mhz')
        error('parameters.cp_range_mhz must be specified for CP ENDOR.');
    end
end
if (~ischar(assumptions))&&(~isstring(assumptions))
    error('assumptions must be a character string or a string scalar.');
end
end

