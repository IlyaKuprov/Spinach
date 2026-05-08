%ENDOR_KEHL_CONTEXT Spinach-style ENDOR context for the Kehl pulse sequences.
%
%   [ENDOR,ENDOR_LB,X,V]=ENDOR_KEHL_CONTEXT(SPIN_SYSTEM,SEQUENCE,...
%   PARAMETERS,ASSUMPTIONS) builds Kehl experimental parameters from the
%   physical fields in PARAMETERS, prepares the Spinach-backed spin
%   operators, EPR orientation selection, ENDOR sweep axes, and line
%   broadening, then calls the specified pulse-sequence kernel.
%
%   For CP ENDOR, the return signature is
%   [ENDOR,ENDOR_LB,X,Y,V].


function varargout=endor_kehl_context(spin_system,sequence,parameters,assumptions)
if nargin<4
    assumptions='labframe';
end

% Set default simulation parameters
parameters=kehl_defaults(parameters);

% Check consistency
grumble(spin_system,sequence,parameters,assumptions);

% Normalise the pulse-sequence selector
seq_name=sequence_name(sequence);

% Build legacy Kehl data structures from Spinach parameters
constants=containers.Map;
constants('CONST1')=2.81;
constants('K_B')=2.0836618E10;
constants('H')=6.62607E-34;
constants('MU_B')=9.27401E-24;
constants('MU_N')=5.050783699E-27;
constants('GE')=28024.95266E6;
constants('GN_1H')=42.576E6;
constants('GN_2D')=6.536E6;
constants('GN_19F')=40.078E6;
constants('GN_14N')=3.077E6;
constants('GN_17O')=-5.772E6;
expt=context_experiment(parameters,constants,seq_name);
spinSys=context_spin_system(spin_system,parameters,constants);
spinOps=kehl_spin_ops(spin_system,parameters.endor_spins,spinSys('N_SpinSys'));

if nargin>=4 && ~isempty(assumptions)
    spin_system=assume(spin_system,assumptions);
end

if strcmp(seq_name,'time')
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

is_cp=strcmp(seq_name,'cp');
if is_cp
    [y_coords,localpar]=cp_axis(localpar);
end

switch seq_name
    case 'mims'
        if parameters.Relax==true
            endor_amp=kehl_mims_rlx(localpar.constants,localpar.spinOps,...
                                    localpar.spinSys,localpar.expt,localpar,...
                                    localpar.paramsENDOR,localpar.epr);
        else
            endor_amp=kehl_mims_calc(localpar.constants,localpar.spinOps,...
                                     localpar.spinSys,localpar.expt,localpar,...
                                     localpar.paramsENDOR,localpar.epr);
        end
    case 'time'
        if parameters.Relax==true
            endor_amp=kehl_time_rlx(localpar.constants,localpar.spinOps,...
                                    localpar.spinSys,localpar.expt,localpar,...
                                    localpar.paramsENDOR,localpar.epr);
        else
            endor_amp=kehl_time_calc(localpar.constants,localpar.spinOps,...
                                     localpar.spinSys,localpar.expt,localpar,...
                                     localpar.paramsENDOR,localpar.epr);
        end
    case 'davies'
        if parameters.Relax==true
            endor_amp=kehl_davies_rlx(localpar.constants,localpar.spinOps,...
                                      localpar.spinSys,localpar.expt,localpar,...
                                      localpar.paramsENDOR,localpar.epr);
        else
            endor_amp=kehl_davies_calc(localpar.constants,localpar.spinOps,...
                                       localpar.spinSys,localpar.expt,localpar,...
                                       localpar.paramsENDOR,localpar.epr);
        end
    case 'spinlock'
        if parameters.Relax==true
            endor_amp=kehl_spinlock_rlx(localpar.constants,localpar.spinOps,...
                                        localpar.spinSys,localpar.expt,localpar,...
                                        localpar.paramsENDOR,localpar.epr);
        else
            endor_amp=kehl_spinlock_calc(localpar.constants,localpar.spinOps,...
                                         localpar.spinSys,localpar.expt,localpar,...
                                         localpar.paramsENDOR,localpar.epr);
        end
    case 'tensor'
        if parameters.Relax==true
            endor_amp=kehl_tensor_rlx(localpar.constants,localpar.spinOps,...
                                      localpar.spinSys,localpar.expt,localpar,...
                                      localpar.paramsENDOR,localpar.epr);
        else
            endor_amp=kehl_tensor_calc(localpar.constants,localpar.spinOps,...
                                       localpar.spinSys,localpar.expt,localpar,...
                                       localpar.paramsENDOR,localpar.epr);
        end
    case 'cp'
        if parameters.Relax==true
            endor_amp=kehl_cp_calc_rlx(localpar.constants,localpar.spinOps,...
                                       localpar.spinSys,localpar.expt,localpar,...
                                       localpar.paramsENDOR,localpar.epr);
        else
            endor_amp=kehl_cp_calc(localpar.constants,localpar.spinOps,...
                                   localpar.spinSys,localpar.expt,localpar,...
                                   localpar.paramsENDOR,localpar.epr);
        end
    otherwise
        error('unknown Kehl ENDOR pulse sequence.');
end
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

function seq_name=sequence_name(sequence)
    if isa(sequence,'function_handle')
        seq_name=func2str(sequence);
    else
        seq_name=char(sequence);
    end
    seq_name=lower(seq_name);
    seq_name=strrep(seq_name,'endor_kehl_','');
    seq_name=strrep(seq_name,'kehl_','');
    seq_name=strrep(seq_name,'_calc','');
    seq_name=strrep(seq_name,'_rlx','');
    if strcmp(seq_name,'timedomain')
        seq_name='time';
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

function expt=context_experiment(parameters,constants,seq_name)
    % Create the primary experiment parameter map
    expt=containers.Map;
    expt('FreqMeas')=parameters.mw_freq_ghz*1e9;
    expt('Field')=parameters.static_field_g*1e-4;
    expt('deltaField')=parameters.field_step_g*1e-4;
    expt('res_EN')=parameters.endor_res_mhz*1e6;
    expt('range_EN')=parameters.endor_range_mhz*1e6;
    expt('t')=parameters.pulse_times_ns*1e-9;

    % Get radiofrequency and microwave nutation field policy
    rf_nutations=[];
    if isfield(parameters,'rf_nutation_freqs')
        rf_nutations=parameters.rf_nutation_freqs;
    end
    if isfield(parameters,'rf_field_from_pulses')
        rf_auto=parameters.rf_field_from_pulses;
    else
        rf_auto=isempty(rf_nutations);
    end

    % Set sequence-specific nutation fields
    t=expt('t');
    switch seq_name
        case 'davies'
            if rf_auto==false
                if size(rf_nutations,1)==3
                    expt('prep')=rf_nutations(1)*2*pi*1e6;
                    expt('oneE')=rf_nutations(2)*2*pi*1e6;
                    expt('oneN')=rf_nutations(3)*2*pi*1e3;
                    expt('pulsewidth')=expt('oneE')/(2*pi*constants('CONST1')*1e10);
                elseif ~isempty(rf_nutations)
                    error('parameters.rf_nutation_freqs has incompatible dimensions.');
                end
            else
                expt('prep')=2*pi/(2*t(1));
                expt('oneE')=2*pi/(4*t(5));
                expt('oneN')=2*pi/(2*t(3));
                expt('pulsewidth')=expt('oneE')/(2*pi*constants('CONST1')*1e10);
            end
        case 'spinlock'
            if rf_auto==false
                if size(rf_nutations,1)==3
                    expt('SL')=rf_nutations(1)*2*pi*1e6;
                    expt('oneE')=rf_nutations(2)*2*pi*1e6;
                    expt('oneN')=rf_nutations(3)*2*pi*1e3;
                    expt('pulsewidth')=expt('SL')/(2*pi*constants('CONST1')*1e10);
                elseif size(rf_nutations,1)==1
                    expt('SL')=rf_nutations(1)*2*pi*1e6;
                    expt('oneE')=2*pi/(4*t(5));
                    expt('oneN')=2*pi/(2*t(3));
                    expt('pulsewidth')=expt('SL')/(2*pi*constants('CONST1')*1e10);
                elseif ~isempty(rf_nutations)
                    error('parameters.rf_nutation_freqs has incompatible dimensions.');
                end
            else
                expt('SL')=2*pi/(4*t(5));
                expt('oneE')=2*pi/(4*t(5));
                expt('oneN')=2*pi/(2*t(3));
                expt('pulsewidth')=expt('SL')/(2*pi*constants('CONST1')*1e10);
            end
        case 'tensor'
            if rf_auto==false
                if size(rf_nutations,1)==2
                    expt('oneE')=rf_nutations(1)*2*pi*1e6;
                    expt('oneN')=rf_nutations(2)*2*pi*1e3;
                    expt('pulsewidth')=expt('oneE')/(2*pi*constants('CONST1')*1e10);
                elseif ~isempty(rf_nutations)
                    error('parameters.rf_nutation_freqs has incompatible dimensions.');
                end
            else
                expt('oneE')=2*pi/(2*t(2));
                expt('oneN')=2*pi/(2*t(1));
                expt('pulsewidth')=expt('oneE')/(2*pi*constants('CONST1')*1e10);
            end
        case 'time'
            if rf_auto==false
                if size(rf_nutations,1)==2
                    expt('oneE')=rf_nutations(1)*2*pi*1e6;
                    expt('oneN')=rf_nutations(2)*2*pi*1e3;
                    expt('pulsewidth')=expt('oneE')/(2*pi*constants('CONST1')*1e10);
                elseif ~isempty(rf_nutations)
                    error('parameters.rf_nutation_freqs has incompatible dimensions.');
                end
            else
                expt('oneE')=2*pi/(4*t(1));
                if isKey(expt,'ang')
                    expt('oneN')=(expt('ang')/180)*2*pi/(4*t(5));
                else
                    expt('oneN')=2*pi/(4*t(5));
                end
                expt('pulsewidth')=expt('oneE')/(2*pi*constants('CONST1')*1e10);
            end
        case 'cp'
            if rf_auto==false
                if size(rf_nutations,2)==5
                    expt('prep')=rf_nutations(1)*2*pi*1e6;
                    expt('SL')=rf_nutations(2)*2*pi*1e6;
                    expt('CP')=rf_nutations(3)*2*pi*1e3;
                    expt('oneE')=rf_nutations(4)*2*pi*1e6;
                    expt('oneN')=rf_nutations(5)*2*pi*1e3;
                    expt('pulsewidth')=expt('prep')/(2*pi*constants('CONST1')*1e10);
                elseif size(rf_nutations,2)==3
                    expt('prep')=rf_nutations(1)*2*pi*1e6;
                    expt('SL')=rf_nutations(2)*2*pi*1e6;
                    expt('CP')=rf_nutations(3)*2*pi*1e3;
                    expt('oneE')=2*pi/(4*t(7));
                    expt('oneN')=2*pi/(2*t(5));
                    expt('pulsewidth')=expt('prep')/(2*pi*constants('CONST1')*1e10);
                elseif size(rf_nutations,2)==2
                    if t(1)==0
                        expt('prep')=2*pi/(4*t(7));
                        expt('SL')=rf_nutations(1)*2*pi*1e6;
                        expt('CP')=rf_nutations(2)*2*pi*1e3;
                        expt('pulsewidth')=expt('SL')/(2*pi*constants('CONST1')*1e10);
                    else
                        expt('prep')=2*pi/(4*t(1));
                        expt('SL')=rf_nutations(1)*2*pi*1e6;
                        expt('CP')=rf_nutations(2)*2*pi*1e3;
                        expt('pulsewidth')=expt('prep')/(2*pi*constants('CONST1')*1e10);
                    end
                    expt('oneE')=2*pi/(4*t(7));
                    expt('oneN')=2*pi/(2*t(5));
                elseif ~isempty(rf_nutations)
                    error('parameters.rf_nutation_freqs has incompatible dimensions.');
                end
            else
                expt('prep')=2*pi/(4*t(7));
                expt('SL')=2*pi/(4*t(7));
                expt('CP')=2*pi/(2*t(5));
                expt('oneE')=2*pi/(4*t(7));
                expt('oneN')=2*pi/(2*t(5));
                expt('pulsewidth')=expt('oneE')/(2*pi*constants('CONST1')*1e10);
            end
        otherwise
            if rf_auto==false
                if size(rf_nutations,1)==2
                    expt('oneE')=rf_nutations(1)*2*pi*1e6;
                    expt('oneN')=rf_nutations(2)*2*pi*1e3;
                    expt('pulsewidth')=expt('oneE')/(2*pi*constants('CONST1')*1e10);
                elseif ~isempty(rf_nutations)
                    error('parameters.rf_nutation_freqs has incompatible dimensions.');
                end
            else
                expt('oneE')=2*pi/(4*t(1));
                if isKey(expt,'ang')
                    expt('oneN')=(expt('ang')/180)*2*pi/(2*t(5));
                else
                    expt('oneN')=2*pi/(2*t(5));
                end
                expt('pulsewidth')=expt('oneE')/(2*pi*constants('CONST1')*1e10);
            end
    end

    % Set the frequency-domain EPR sweep
    if isfield(parameters,'epr_freq_min_ghz')
        expt('FreqMin')=parameters.epr_freq_min_ghz*1e9;
        expt('FreqRange')=parameters.epr_freq_range_ghz*1e9;
        expt('FreqSteps')=parameters.epr_freq_step_ghz*1e9;
    end

    % Set the ENDOR radiofrequency pulse flip angle
    if isfield(parameters,'rf_flip_angle_deg')
        expt('ang')=parameters.rf_flip_angle_deg;
    end

    % Set the shaped pulse profile file
    if isfield(parameters,'pulse_file')
        multipulses=false;
        if isfield(parameters,'multipulses')
            multipulses=parameters.multipulses;
        end
        expt('pulse')=parameters.pulse_file;
        expt('3pulses')=multipulses;
    end

    % Set the CP sweep dimension when present
    if isfield(parameters,'cp_start_mhz')
        expt('start_CP')=parameters.cp_start_mhz*1e6;
        expt('Npts_CP')=parameters.cp_npoints;
        expt('range_CP')=parameters.cp_range_mhz*1e6;
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
function grumble(spin_system,sequence,parameters,assumptions)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if (~isa(sequence,'function_handle'))&&(~ischar(sequence))&&(~isstring(sequence))
    error('sequence must be a function handle, character string, or string scalar.');
end
seq_name=sequence_name(sequence);
if ~ismember(seq_name,{'mims','time','davies','spinlock','tensor','cp'})
    error('sequence must specify mims, time, davies, spinlock, tensor, or cp.');
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
if strcmp(seq_name,'cp')
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
% Localised single-call Kehl helper functions


% --- Localised from kehl_defaults.m ---
% supplies default simulation parameters for Kehl ENDOR experiments
%
% input parameters:
% parameters: structure containing ENDOR context and simulation parameters
%
% output parameters:
% parameters: structure with missing simulation parameters populated
%
% Febuary 2024 A. Kehl (akehl@gwdg.de)
% Spinach architecture migration May 2026 Talos
%

function parameters=kehl_defaults(parameters)
if nargin<1
    parameters=struct();
end

% Check consistency
kehl_defaults_grumble(parameters);

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

function kehl_defaults_grumble(parameters)
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
end


% --- Localised from kehl_prep_time.m ---
% calculates/sets necessary parameters for the time domain ENDOR calculation from
% the experimental values
%
% input parameters:
% constants:the Map containing the constants
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
%
% output parameters:
% paramsENDOR: the Map containing the ENDOR parameters
%
% February 2024 A. Kehl (akehl@gwdg.de)

% Larmor frequencies

function paramsENDOR=kehl_prep_time(constants,spinSys,expt)

    % Check consistency
    kehl_prep_time_grumble(constants,spinSys,expt);
    Nuclei=spinSys("Nuclei");
    v_L=zeros(size(Nuclei,2),1);
    for i=1:size(Nuclei,2)
        v_L(i)=kehl_nuc_gamma(constants,Nuclei{i})*expt("Field");
    end

    % number of points and values for x-axis

    % input in us, factor 10e-6 due to usual conversion MHz to Hz
    range_EN=expt("range_EN")*1e-12;

    % input in us, factor 10e-6 due to usual conversion MHz to Hz
    res_EN=expt("res_EN")*1e-12;
    Npts_EN=(round(range_EN/res_EN));

    % ENDOR start x-axis
    start_EN=expt("t_start");

    % X-Axis steps
    step_EN=range_EN/(Npts_EN-1);
    x_coords=zeros(Npts_EN);
    for ii=1:Npts_EN

        %-v_L;
        x_coords(ii)=start_EN+(ii-1)*step_EN;
    end

    paramsENDOR=containers.Map;

    paramsENDOR("v_L")=v_L;
    paramsENDOR("start_EN")=start_EN;
    paramsENDOR("step_EN")=step_EN;
    paramsENDOR("range_EN")=expt("range_EN");
    paramsENDOR("Npts_EN")=Npts_EN;
    paramsENDOR("x_coords")=x_coords;



end

function kehl_prep_time_grumble(constants,spinSys,expt)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
end


% --- Localised from kehl_prep_endor.m ---
% calculates/sets necessary parameters for the ENDOR calculation from
% the experimental values
%
% input parameters:
% constants:the Map containing the constants
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
%
% output parameters:
% paramsENDOR: the Map containing the ENDOR parameters
%
% February 2024 A. Kehl (akehl@gwdg.de)

% Larmor frequencies

function paramsENDOR=kehl_prep_endor(constants,spinSys,expt)

    % Check consistency
    kehl_prep_endor_grumble(constants,spinSys,expt);
    Nuclei=spinSys("Nuclei");
    v_L=zeros(size(Nuclei,2),1);
    for i=1:size(Nuclei,2)
        v_L(i)=kehl_nuc_gamma(constants,Nuclei{i})*expt("Field");
    end

    % number of points and values for x-axis
    Npts_EN=(round(expt("range_EN")/expt("res_EN")));
    if isKey(expt,"RF_start")

        % ENDOR start x-axis
        start_EN=expt("RF_start");
    else

        % ENDOR start x-axis
        start_EN=v_L(1)-expt("range_EN")/2;
    end

    % X-Axis steps
    step_EN=expt("range_EN")/(Npts_EN-1);
    x_coords=zeros(Npts_EN);
    for ii=1:Npts_EN

        %-v_L;
        x_coords(ii)=start_EN+(ii-1)*step_EN;
    end

    paramsENDOR=containers.Map;

    paramsENDOR("v_L")=v_L;
    paramsENDOR("start_EN")=start_EN;
    paramsENDOR("step_EN")=step_EN;
    paramsENDOR("range_EN")=expt("range_EN");
    paramsENDOR("Npts_EN")=Npts_EN;
    paramsENDOR("x_coords")=x_coords;


end

function kehl_prep_endor_grumble(constants,spinSys,expt)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
end


% --- Localised from kehl_prep_epr.m ---
% calculates/sets necessary parameters for the EPR calculation from
% the experimental values
%
% input parameters:
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
%
% output parameters:
% paramsEPR: the Map containing the EPR parameters
%
% February 2024 A. Kehl (akehl@gwdg.de)



function paramsEPR=kehl_prep_epr(spinSys,expt)

    % Check consistency
    kehl_prep_epr_grumble(spinSys,expt);
    g_iso=spinSys("g_iso");
    obsField=expt("Field");

    gBepr=obsField*g_iso;
    geff=g_iso;
    B=gBepr/geff;
    veff=geff*obsField*9.27401*1e-24/(6.62607*1e-34);

    % EPR simulation range
    fieldCenter=gBepr/g_iso;
    fieldmin=fieldCenter-0.260;
    fieldmax=fieldCenter+0.260;

    % EPR x-axis definition

    % no of points in Field dimension
    Npts=(round((fieldmax-fieldmin)/expt("deltaField"))+1);
    field=linspace(fieldmin,fieldmax,Npts);

    paramsEPR=containers.Map;
    paramsEPR("Npts")=Npts;
    paramsEPR("fieldAxis")=field;

end

function kehl_prep_epr_grumble(spinSys,expt)
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
end


% --- Localised from kehl_ori_field.m ---
% calculates the selected EPR orientations in the field domain and the
% corresponding effective spin parameter values
%
% input parameters:
% constants: the Map containing the constants
% spinSys: the Map describing the spin system
% spinOps: the Map containing the spin operators
% paramsEPR: the Map containing the EPR parameters
% parameters: structure containing simulation parameters
% expt: the Map containing the experimental parameters
%
% output parameters:
% EPR: Map containing the information from the EPR experiment
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

% get parameters from Maps

function EPR=kehl_ori_field(constants,spinSys,spinOps,paramsEPR,paramsENDOR,parameters,expt)

    % Check consistency
    kehl_ori_field_grumble(constants,spinSys,spinOps,paramsEPR,paramsENDOR,parameters,expt);
    Ntheta=parameters.Nang;
    Nphimax=parameters.Nang;
    g=spinSys("g");
    A=spinSys("A");
    Q=spinSys("Q");
    if spinSys("EPR_Nucs_used")==true
        Ni_EPR=spinSys("Ni_EPR");
        A_EPR=spinSys("A_EPR");
        Q_EPR=spinSys("Q_EPR");
    else
        Ni_EPR=0;
    end
    if spinSys("CS_used")==true
        CS=spinSys("CS");
    end
    if spinSys("D_used")==true
        D=spinSys("D");
    end
    Ni_ENDOR=spinSys("Ni_ENDOR");
    g2=g*g;

    % initialize arrays
    tmp_epr=zeros(paramsEPR("Npts"),1);
    epr_amp=zeros(paramsEPR("Npts"),1);

    or=0;
    geff_sel=[];
    B_sel=[];
    HF_zz_sel=[];
    HF_zy_sel=[];
    HF_zx_sel=[];
    NQI_zz_sel=[];
    NQI_sel=[];
    CS_zz_sel=[];
    D_zz_sel=[];
    S_sel=[];
    offsets_sel=[];


    % for powder pattern
    if parameters.powder==true
        if isKey(expt,'exciteWidth')
            W1=expt('exciteWidth');
        else
            W1=expt("pulsewidth")*constants("CONST1")*1e10;
        end
        nor=0;
    % loop through orientations
    for ii=1:Ntheta
        theta=ii*pi/Ntheta;
        Nphi=round(sin(theta)*Nphimax)*1;
        for jj=1:Nphi
            phi=(jj-1)*pi*2/(Nphi);

            % direction cosine vector
            dc=[cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
            % Effective g-value for a given theta, and phi combination
            geff=(dc*g2*dc')^.5;

            % effective B field for given theta and phi
            B=(expt("Field")*spinSys("g_iso"))/geff(1);
            Beff=expt('FreqMeas')*constants('H')/(constants('MU_B')*geff);

            % HF ENDOR
            HF_zz=zeros(1,Ni_ENDOR);
            HF_zx=zeros(1,Ni_ENDOR);
            HF_zy=zeros(1,Ni_ENDOR);

            for m=1:Ni_ENDOR

                 HF_zz(m)=(sin(theta))^2*(cos(phi))^2*A(3*m-2,1)...
+(sin(theta))^2*(sin(phi))^2*A(3*m-1,2)...
+(cos(theta))^2*A(3*m,3)...
+2*(sin(theta))^2*sin(phi)*cos(phi)*A(3*m-2,2) ...
+2*sin(theta)*cos(theta)*cos(phi)*A(3*m-2,3)...
+2*sin(theta)*cos(theta)*sin(phi)*A(3*m-1,3);
                    %A11+A22+A33+A12+A13+A23

                 if parameters.Bterm==1
                        HF_zx(m)=(A(3*m-2,1)*sin(theta)*cos(phi)...
+A(3*m-1,1)*sin(theta)*sin(phi)...
+A(3*m,1)*cos(theta))*cos(theta)*cos(phi)...
+(A(3*m-2,2)*sin(theta)*cos(phi)...
+A(3*m-1,2)*sin(theta)*sin(phi)...
+A(3*m,2)*cos(theta))*cos(theta)*sin(phi)...
-(A(3*m-2,3)*sin(theta)*cos(phi)...
+A(3*m-1,3)*sin(theta)*sin(phi)...
+A(3*m,3)*cos(theta))*sin(theta);

                        HF_zy(m)=-A(3*m-2,1)*sin(theta)*cos(phi)*sin(phi)...
-A(3*m-1,1)*sin(theta)*sin(phi)^2 ...
-A(3*m,1)*cos(theta)*sin(phi)...
+A(3*m-2,2)*sin(theta)*cos(phi)^2 ...
+A(3*m-1,2)*sin(theta)*sin(phi)*cos(phi)...
+A(3*m,2)*cos(theta)*cos(phi);
                end
            end

            % NQI ENDOR
            NQI_zz=zeros(1,Ni_ENDOR);
            NQI=zeros(Ni_ENDOR,3,3);

            for m=1:Ni_ENDOR
                if spinSys("Q_used")==true

                    NQI_zz(m)=(sin(theta))^2*(cos(phi))^2*Q(3*m-2,1)...
+(sin(theta))^2*(sin(phi))^2*Q(3*m-1,2)...
+(cos(theta))^2*Q(3*m,3)...
+2*(sin(theta))^2*sin(phi)*cos(phi)*Q(3*m-2,2)...
+2*sin(theta)*cos(theta)*cos(phi)*Q(3*m-2,3)...
+2*sin(theta)*cos(theta)*sin(phi)*Q(3*m-1,3);
                    % Q11,Q22,Q33,Q12,Q13,Q23

                    % Rotation matrix into lab system
                    R1=zeros(3);
                    R1(1,1)=cos(theta)*cos(phi);
                    R1(1,2)=cos(theta)*sin(phi);
                    R1(1,3)=-sin(theta);
                    R1(2,1)=-sin(phi);
                    R1(2,2)=cos(phi);
                    R1(2,3)=0;
                    R1(3,1)=sin(theta)*cos(phi);
                    R1(3,2)=sin(theta)*sin(phi);
                    R1(3,3)=cos(theta);

                    Q_ENDOR=spinSys("Q");
                    qq2=Q_ENDOR((m-1)*3+1:(m-1)*3+3,:);
                    Y2=R1*qq2*R1';
                    NQI(m,:,:)=Y2;

                else
                    NQI_zz(m)=0;
                    NQI (m,:,:)=0;
                end
            end

            % CS ENDOR
            CS_zz=zeros(1,Ni_ENDOR);
            if spinSys("CS_used")==true
                for m=1:Ni_ENDOR
                    CS_zz(m)=(sin(theta))^2*(cos(phi))^2*CS(3*m-2,1)...
+(sin(theta))^2*(sin(phi))^2*CS(3*m-1,2)...
+(cos(theta))^2*CS(3*m,3)...
+2*(sin(theta))^2*sin(phi)*cos(phi)*CS(3*m-2,2)...
+2*sin(theta)*cos(theta)*cos(phi)*CS(3*m-2,3)...
+2*sin(theta)*cos(theta)*sin(phi)*CS(3*m-1,3);
                end
            end

            % dipolar SSC ENDOR
            D_zz=zeros(1,Ni_ENDOR);

             if spinSys("D_used")==true
                 D_zz=zeros(size(D,1)/3);
                 for m=1:size(D,1)/3
                    D_zz(m)=(sin(theta))^2*(cos(phi))^2*D(3*m-2,1)...
+(sin(theta))^2*(sin(phi))^2*D(3*m-1,2)...
+(cos(theta))^2*D(3*m,3)...
+2*(sin(theta))^2*sin(phi)*cos(phi)*D(3*m-2,2)...
+2*sin(theta)*cos(theta)*cos(phi)*D(3*m-2,3)...
+2*sin(theta)*cos(theta)*sin(phi)*D(3*m-1,3);
                 end
             end

            % HF EPR
            HF_zz_EPR=zeros(1,Ni_EPR);
            if Ni_EPR>0
                for m=1:Ni_EPR
                    HF_zz_EPR(m)=(sin(theta))^2*(cos(phi))^2*A_EPR(3*m-2,1)...
+(sin(theta))^2*(sin(phi))^2*A_EPR(3*m-1,2)...
+(cos(theta))^2*A_EPR(3*m,3)...
+2*(sin(theta))^2*sin(phi)*cos(phi)*A_EPR(3*m-2,2)...
+2*sin(theta)*cos(theta)*cos(phi)*A_EPR(3*m-2,3)...
+2*sin(theta)*cos(theta)*sin(phi)*A_EPR(3*m-1,3);
                end
            end


            fieldAxis=paramsEPR('fieldAxis');

            if Ni_EPR>0
                I_EPR=spinSys("I_EPR");
                mI=kehl_ori_field_BuildSpace(I_EPR);
                for i=1:Ni_EPR
                        % EPR resonance with hyperfine coupling only
                        E=Beff+mI*(HF_zz_EPR(i)/constants("CONST1")*10^(-10))';
                end
                bin=(round((E-fieldAxis(1))/expt("deltaField"))+1);
            else
                E=Beff;
                bin=(round((E-fieldAxis(1))/expt("deltaField"))+1);
            end

            for p=1:length(bin)

                % Add resonances to EPR spectrum
                tmp_epr(bin(p))=tmp_epr(bin(p))+1;
            end

            offsets=kehl_offsets(constants,spinSys,spinOps,paramsENDOR,B,geff,HF_zz,NQI_zz);

            DeltaB=0;
            scalefactor=0;

            for l=1:length(E)
                nor=nor+1;

                % Magnetic field offset in T
                DeltaB=fieldAxis(bin(l))-expt("Field");

                if isKey(expt,"pulse")
                    scalefactor=kehl_ori_field_pulsescale(expt,DeltaB,constants("CONST1"));
                else

                    if abs(DeltaB)<=abs(parameters.nwidth*(W1/(constants("CONST1")*1e10)))

                        %exp('res_EN')
                        if abs(HF_zz(m))>1
                            Scale=(DeltaB^2-(W1/(constants("CONST1")*1e10))^2)/(DeltaB^2+(W1/(constants("CONST1")*1e10))^2);
                            Scale=(1-Scale)/2;
                        else
                            Scale=0;
                        end
                    else
                       Scale=0;
                    end
                    scalefactor=Scale;
                end


            %Select only those parameters, for which scalefactor > 0
            if scalefactor>1e-3
                or=or+1;
                geff_sel(or)=geff(1);
                B_sel(or)=Beff;
                HF_zz_sel(or,:)=HF_zz(:);
                HF_zy_sel(or,:)=HF_zy(:);
                HF_zx_sel(or,:)=HF_zx(:);
                NQI_zz_sel(or,:)=NQI_zz(:);
                NQI_sel(or,:,:,:)=NQI(:,:,:);
                CS_zz_sel(or,:)=CS_zz(:);
                D_zz_sel(or,:)=D_zz(:);

                S_sel(or)=scalefactor;
                offsets_sel(or,:)=offsets(:);
            end
            end
        end
        epr_amp=epr_amp+tmp_epr;
    end

    % only one orientation for single crystal calculation
    elseif parameters.powder==false
        or=or+1;
        geff=spinSys("g_iso");

        %effective B field for given theta and phi
        Beff=expt('FreqMeas')*constants('H')/(constants('MU_B')*geff);

        HF_zz=zeros(1,Ni_ENDOR);
        HF_zx=zeros(1,Ni_ENDOR);
        HF_zy=zeros(1,Ni_ENDOR);
        NQI_zz=zeros(1,Ni_ENDOR);
        NQI=zeros(Ni_ENDOR,3,3);


        CS_zz=zeros(1,Ni_ENDOR);

        for m=1:Ni_ENDOR

            HF_zz(m)=A(3*m,3);
            if parameters.Bterm==1
                HF_zy(m)=A(3*m,2);
                HF_zx(m)=A(3*m,1);
            end

            if spinSys("Q_used")==true
                NQI_zz(m)=Q(3*m,3);
                NQI(m,:,:)=Q;
            end

            % CS ENDOR
            if spinSys("CS_used")==true
                CS_zz(m)=CS(3*m,3);
            end

        end

        % dipolar SSC ENDOR
        D_zz=zeros(1,Ni_ENDOR);

        if spinSys("D_used")==true
             D_zz=zeros(1,size(D,1)/3);
             for m=1:size(D,1)/3
                 D_zz(m)=D(3*m,3);
             end
         end


        offsets_tmp=kehl_offsets(constants,spinSys,spinOps,paramsENDOR,Beff,geff,HF_zz,NQI_zz);

        sel_I=parameters.sel_I;
        if spinSys('N_SpinSys')>1
            s=size(offsets_tmp,2)/Ni_ENDOR;
            for n=1:Ni_ENDOR
                offsets(n)=offsets_tmp(sel_I+(n-1)*s);
            end
        else
            size(offsets_tmp,2)
            s=size(offsets_tmp,2)/Ni_ENDOR;
            offsets=offsets_tmp(((sel_I-1)*s+1):(sel_I*s));
        end

        scalefactor=1;
        geff_sel(or)=geff(1);
        B_sel(or)=Beff;
        HF_zz_sel(or,:)=HF_zz(:);
        HF_zy_sel(or,:)=HF_zy(:);
        HF_zx_sel(or,:)=HF_zx(:);
        NQI_zz_sel(or,:)=NQI_zz(:);
        NQI_sel(or,:,:,:)=NQI(:,:,:);
        CS_zz_sel(or,:)=CS_zz(:);
        D_zz_sel(or,:)=D_zz(:);

        S_sel(or)=scalefactor;
        offsets_sel(or,:)=offsets(:);
    end
    % build output Map
    EPR=containers.Map;
    EPR("geff_sel")=geff_sel;
    EPR("B_sel")=B_sel;
    EPR("HF_zz_sel")=HF_zz_sel;
    EPR("HF_zy_sel")=HF_zy_sel;
    EPR("HF_zx_sel")=HF_zx_sel;
    EPR("CS_zz_sel")=CS_zz_sel;
    EPR("D_zz_sel")=D_zz_sel;

    EPR("NQI_zz_sel")=NQI_zz_sel;
    EPR("NQI_sel")=NQI_sel;

    EPR("S_sel")=S_sel;
    EPR("EPR_amp")=epr_amp;

    EPR("offsets")=offsets_sel;
end

function scalefactor=kehl_ori_field_pulsescale(expt,DeltaB,CONST1)
    % calculates the scalefactor in dependence of the mw pulse's excitation
    % profile
    %
    % input parameters:
    % expt: the Map containing the experimental parameters
    % DeltaB: offset of the actual field from the effective resonance field
    % CONST1: constant to convert between field and frequency
    %
    % output parameters:
    % scalefactor: scalefactor for the specific orientation
    %
    % February 2024 A. Kehl (akehl@gwdg.de)
    %

    data=expt("pulse");
    pulse=fopen(data);

    pulse_data=textscan(pulse,'%f %f %f %f');
    fclose('all');
    pulse_x=pulse_data{2};
    dx=1000/(pulse_x(2)-pulse_x(1));

    n=1000;


    pulse_y=pulse_data{3}+1i*pulse_data{4};

    Y=abs(fftshift(fft(pulse_y,(n+1))));
    y=Y.^3;
    y=y/max(y);
    Y=Y/max(Y);
    if isKey(expt,"3pulses") && expt("3pulses")==true
       Y=y;
    end

    %Delta omega/2pi in MHz
    DeltaOM=DeltaB*CONST1*10^4*2*pi;

    binDOM=round(DeltaOM)+(dx/2+1);
    if binDOM>0 && binDOM<(n+1)
        scalefactor=Y(binDOM);
    else
        scalefactor=0;
    end

end

function M=kehl_ori_field_BuildSpace(S)
    % builds up the space of basis vectors in columns,
    % script by M. Bennati
    %
    % input parameters:
    % S: spin quantum number
    %
    % output parameters:
    % M: basis vectors in columns
    %
    % February 2024 A. Kehl (akehl@gwdg.de)
    %

    % Initialize the space to a null
    M=[];

    % Initialize the dimension to one
    dim=1;

    % Loop over each spin
    for j=1:length(S)

      % Initialize the holder
      temp=[];

      % Loop over the z components
      for k=-S(j):S(j)

        %   adding block of the old space
        temp=[temp; M k*ones(dim,1)];

      %   plus a column of the z component
      end

      % Update M
      M=temp;

      % The dimensionality is the number
      dim=size(M,1);
                                            %   of rows in M
    end
end

function kehl_ori_field_grumble(constants,spinSys,spinOps,paramsEPR,paramsENDOR,parameters,expt)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if ~isa(paramsEPR,'containers.Map')
    error('paramsEPR must be a containers.Map object.');
end
if ~isa(paramsENDOR,'containers.Map')
    error('paramsENDOR must be a containers.Map object.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
end


% --- Localised from kehl_ori_freq.m ---
% calculates the selected EPR orientations in the frequency domain and the
% corresponding effective spin parameter values
%
% input parameters:
% constants: the Map containing the constants
% spinSys: the Map describing the spin system
% spinOps: the Map containing the spin operators
% paramsEPR: the Map containing the EPR parameters
% parameters: structure containing simulation parameters
% expt: the Map containing the experimental parameters
%
% output parameters:
% EPR: Map containing the information from the EPR experiment
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

% get parameters from Maps

function EPR=kehl_ori_freq(constants,spinSys,spinOps,paramsEPR,paramsENDOR,parameters,expt)

    % Check consistency
    kehl_ori_freq_grumble(constants,spinSys,spinOps,paramsEPR,paramsENDOR,parameters,expt);
    Ntheta=parameters.Nang;
    Nphimax=parameters.Nang;
    g=spinSys("g");
    A=spinSys("A");
    Q=spinSys("Q");
    if spinSys("EPR_Nucs_used")==true
        Ni_EPR=spinSys("Ni_EPR");
        A_EPR=spinSys("A_EPR");
        Q_EPR=spinSys("Q_EPR");
        g_N_EPR=spinSys("g_N_EPR");
        I_EPR=spinSys("I_EPR");
        spinOpsEPR=kehl_spin_ops(spinSys("S"),I_EPR,Ni_EPR,Ni_EPR);
        SzEPR=spinOpsEPR("Sz");
        SxEPR=spinOpsEPR("Sx");
        IzEPR=spinOpsEPR("Iz");
        IxEPR=spinOpsEPR("Ix");
        IyEPR=spinOpsEPR("Iy");
    else
        Ni_EPR=0;
        spinOpsEPR=kehl_spin_ops(spinSys("S"),zeros(0,1),0,1);
        SzEPR=spinOpsEPR("Sz");
        SxEPR=spinOpsEPR("Sx");
    end
    if spinSys("CS_used")==true
        CS=spinSys("CS");
    end
    if spinSys("D_used")==true
        D=spinSys("D");
    end
    Ni_ENDOR=spinSys("Ni_ENDOR");
    g2=g*g;

    % initialize arrays
    tmp_epr=zeros(paramsEPR("Npts"),1);
    epr_amp=zeros(paramsEPR("Npts"),1);

    or=0;
    geff_sel=[];
    B_sel=[];
    HF_zz_sel=[];
    HF_zy_sel=[];
    HF_zx_sel=[];
    NQI_zz_sel=[];
    NQI_sel=[];
    CS_zz_sel=[];
    D_zz_sel=[];
    offsets_sel=[];

    S_sel=[];

    nor=0;

    % loop over orientations for powder pattern
    if parameters.powder==true
        if isKey(expt,'exciteWidth')
            W1=expt('exciteWidth');
        else
            W1=expt("pulsewidth")*constants("CONST1")*1e10;
        end

        for ii=1:Ntheta
            theta=ii*pi/Ntheta;
            Nphi=round(sin(theta)*Nphimax)*1;
            for jj=1:Nphi
                phi=(jj-1)*pi*2/(Nphi);

                % direction cosine vector
                dc=[cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
                % Effective g-value for a given theta, and phi combination
                geff=(dc*g2*dc')^.5;

                % effective B field for given theta and phi
                Beff=(expt("Field")*spinSys("g_iso"))/geff(1);
                B=expt("Field");

                % Resonance frequency for geff at ObsField in GHz
                veff=geff*expt("Field")*9.27401*10^-24/(6.62607*10^-34);

                % Rotation matrix into lab system
                R1=zeros(3);

                R1(1,1)=cos(theta)*cos(phi);
                R1(1,2)=cos(theta)*sin(phi);
                R1(1,3)=-sin(theta);
                R1(2,1)=-sin(phi);
                R1(2,2)=cos(phi);
                R1(2,3)=0;
                R1(3,1)=sin(theta)*cos(phi);
                R1(3,2)=sin(theta)*sin(phi);
                R1(3,3)=cos(theta);

                grot=R1*g*R1';


                % ENDOR values
                HF_zz=zeros(1,Ni_ENDOR);
                HF_zx=zeros(1,Ni_ENDOR);
                HF_zy=zeros(1,Ni_ENDOR);
                NQI_zz=zeros(1,Ni_ENDOR);
                NQI=zeros(Ni_ENDOR,3,3);
                CS_zz=zeros(1,Ni_ENDOR);
                D_zz=zeros(1,Ni_ENDOR);

                for m=1:Ni_ENDOR
                     % HF ENDOR
                     HF_zz(m)=(sin(theta))^2*(cos(phi))^2*A(3*m-2,1)...
+(sin(theta))^2*(sin(phi))^2*A(3*m-1,2)...
+(cos(theta))^2*A(3*m,3)...
+2*(sin(theta))^2*sin(phi)*cos(phi)*A(3*m-2,2) ...
+2*sin(theta)*cos(theta)*cos(phi)*A(3*m-2,3)...
+2*sin(theta)*cos(theta)*sin(phi)*A(3*m-1,3);
                        %A11+A22+A33+A12+A13+A23
                     if parameters.Bterm==1
                            HF_zx(m)=(A(3*m-2,1)*sin(theta)*cos(phi)...
+A(3*m-1,1)*sin(theta)*sin(phi)...
+A(3*m,1)*cos(theta))*cos(theta)*cos(phi)...
+(A(3*m-2,2)*sin(theta)*cos(phi)...
+A(3*m-1,2)*sin(theta)*sin(phi)...
+A(3*m,2)*cos(theta))*cos(theta)*sin(phi)...
-(A(3*m-2,3)*sin(theta)*cos(phi)...
+A(3*m-1,3)*sin(theta)*sin(phi)...
+A(3*m,3)*cos(theta))*sin(theta);

                            HF_zy(m)=-A(3*m-2,1)*sin(theta)*cos(phi)*sin(phi)...
-A(3*m-1,1)*sin(theta)*sin(phi)^2 ...
-A(3*m,1)*cos(theta)*sin(phi)...
+A(3*m-2,2)*sin(theta)*cos(phi)^2 ...
+A(3*m-1,2)*sin(theta)*sin(phi)*cos(phi)...
+A(3*m,2)*cos(theta)*cos(phi);
                     end

                     % NQI ENDOR
                     if spinSys("Q_used")==true

                        NQI_zz(m)=(sin(theta))^2*(cos(phi))^2*Q(3*m-2,1)...
+(sin(theta))^2*(sin(phi))^2*Q(3*m-1,2)...
+(cos(theta))^2*Q(3*m,3)...
+2*(sin(theta))^2*sin(phi)*cos(phi)*Q(3*m-2,2)...
+2*sin(theta)*cos(theta)*cos(phi)*Q(3*m-2,3)...
+2*sin(theta)*cos(theta)*sin(phi)*Q(3*m-1,3);
                        % Q11,Q22,Q33,Q12,Q13,Q23

                        % Q_ENDOR = spinSys("Q");
                        qq2=Q((m-1)*3+1:(m-1)*3+3,:);
                        Y2=R1*qq2*R1';
                        NQI(m,:,:)=Y2;
                     end

                     % CS ENDOR
                     if spinSys("CS_used")==true
                        CS_zz(m)=(sin(theta))^2*(cos(phi))^2*CS(3*m-2,1)...
+(sin(theta))^2*(sin(phi))^2*CS(3*m-1,2)...
+(cos(theta))^2*CS(3*m,3)...
+2*(sin(theta))^2*sin(phi)*cos(phi)*CS(3*m-2,2)...
+2*sin(theta)*cos(theta)*cos(phi)*CS(3*m-2,3)...
+2*sin(theta)*cos(theta)*sin(phi)*CS(3*m-1,3);
                     end

                end
                % nulcear dipolar-dipolar coupling ENDOR
                 if spinSys("D_used")==true
                     D_zz=zeros(1,size(D,1)/3);
                     for m=1:size(D,1)/3
                        D_zz(m)=(sin(theta))^2*(cos(phi))^2*D(3*m-2,1)...
+(sin(theta))^2*(sin(phi))^2*D(3*m-1,2)...
+(cos(theta))^2*D(3*m,3)...
+2*(sin(theta))^2*sin(phi)*cos(phi)*D(3*m-2,2)...
+2*sin(theta)*cos(theta)*cos(phi)*D(3*m-2,3)...
+2*sin(theta)*cos(theta)*sin(phi)*D(3*m-1,3);
                     end
                 end

                % EPR vlaues

                HF_EPR=zeros(3,3,Ni_EPR);
                NQI_EPR=zeros(3,3,Ni_EPR);
                % CS_EPR = zeros(3,3,Ni_EPR);

                if Ni_EPR>0
                    for m=1:Ni_EPR
                        % HF EPR
                        hf2=A_EPR((m-1)*3+1:(m-1)*3+3,:);
                        X2=R1*hf2*R1';
                        HF_EPR(:,:,m)=X2;


                         % Q EPR
                        if spinSys("EPR_Q_used")==true
                            qq2=Q_EPR((m-1)*3+1:(m-1)*3+3,:);
                            Y2=R1*qq2*R1';
                            NQI_EPR(:,:,m)=Y2;
                        end

                    end
                end

               H_EZ_EPR=constants("MU_B")*grot(3,3)*B/constants("H")*SzEPR;
               H_NZ_EPR=zeros(size(SzEPR));
               H_HF_EPR=zeros(size(SzEPR));
               H_NQI_EPR=zeros(size(SzEPR));

               if spinSys('EPR_Nucs_used')>0
                   for mm=1:Ni_EPR
                        H_NZ_EPR=H_NZ_EPR+g_N_EPR(mm)*IzEPR{mm};
                        H_HF_EPR=H_HF_EPR+IzEPR{mm}*HF_EPR(3,3,mm)*SzEPR+IxEPR{mm}*(HF_EPR(1,3,mm)^2+HF_EPR(2,3,mm)^2)^0.5*SzEPR;
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(1,1,mm)*IxEPR{mm}*IxEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(1,2,mm)*IxEPR{mm}*IyEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(1,3,mm)*IxEPR{mm}*IzEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(2,1,mm)*IyEPR{mm}*IxEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(2,2,mm)*IyEPR{mm}*IyEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(2,3,mm)*IyEPR{mm}*IzEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(3,1,mm)*IzEPR{mm}*IxEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(3,2,mm)*IzEPR{mm}*IyEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(3,3,mm)*IzEPR{mm}*IzEPR{mm};
                   end
               end

               H_S_EPR=H_EZ_EPR-H_NZ_EPR+H_HF_EPR+H_NQI_EPR;


               % Diagonalize and calculate eigenvalues of the EPR Hamiltonian

               % E ... Hamiltonian Elements, V... eigenvectors
               [V_EPR,E_EPR]=eig(H_S_EPR);

               % round necessary for following calculation
               V_EPR=round((V_EPR),9);
               E_EPR=real(diag(E_EPR));


               % Calculate transitions between all elements
                trans_prob_EPR=[];
                freq_EPR=[];
                q=1;
                for x=1:length(V_EPR)
                    for y=x+1:length(V_EPR)
                        % Electron transition probability in the eigenbasis
                        trans_prob_EPR(q)=abs(round((V_EPR(:,x))'*SxEPR*(V_EPR(:,y)),9))^2;

                        %in Hz, Transition Frequency
                        freq_EPR(q)=abs(E_EPR(x)-E_EPR(y));

                        %Select EPR Transitions with frequency threshhold value (1 GHz)
                        if (freq_EPR(q)<1e9) || (trans_prob_EPR(q)<0.1)
                            trans_prob_EPR(q)=0;
                        end
                        q=q+1;
                    end
                end

                freq_EPR=freq_EPR(logical(trans_prob_EPR))   ;
                trans_prob_EPR=trans_prob_EPR(logical(trans_prob_EPR));

                % EPR_Resonances =[freq_EPR' trans_prob_EPR' grot(3,3)*ones(length(freq_EPR),1)];

                % EPR Frequency Spectrum
                for p=1:length(freq_EPR)

                    %scale resonance to freq axis bin
                    bin_freq=round((freq_EPR-expt("FreqMin"))/expt("FreqSteps"))+1;
                    if trans_prob_EPR(p)>0
                        tmp_epr(bin_freq(p))=tmp_epr(bin_freq(p))+trans_prob_EPR(p);
                        DeltaOm=freq_EPR(p)-expt("FreqMeas");

                        if isKey(expt,"pulse")
                            SF=kehl_ori_freq_pulsescale(expt,DeltaOm)/2;
                        else
                            SF=((DeltaOm^2)-(W1^2))/((DeltaOm^2)+(W1^2));
                            SF=(1-SF)/2;
                            nor=nor+1;
                        end

                        if SF>1e-3

                            % including transition probability!
                            scalefactor=SF;
                        else
                            scalefactor=0;
                        end
                    end


                    %Select only those parameters, for which scalefactor > 0

                    offsets=kehl_offsets(constants,spinSys,spinOps,paramsENDOR,B,geff,HF_zz,NQI_zz);


                    if (scalefactor>0) && (trans_prob_EPR(p)>0)

                        or=or+1;
                        geff_sel(or)=geff(1);
                        B_sel(or)=Beff;
                        HF_zz_sel(or,:)=HF_zz(:);
                        HF_zy_sel(or,:)=HF_zy(:);
                        HF_zx_sel(or,:)=HF_zx(:);
                        NQI_zz_sel(or,:)=NQI_zz(:);
                        NQI_sel(or,:,:,:)=NQI(:,:,:);
                        CS_zz_sel(or,:)=CS_zz(:);
                        D_zz_sel(or,:)=D_zz(:);
                        S_sel(or)=scalefactor;
                        offsets_sel(or,:)=offsets(:);
                    end
                end
            end
            epr_amp=epr_amp+tmp_epr;
        end

    % single orientation for single crystal
    elseif parameters.powder==false
        or=or+1;
        geff=spinSys("g_iso");

        %effective B field for given theta and phi
        B=(expt("Field")*spinSys("g_iso"))/geff(1);

        HF_zz=zeros(1,Ni_ENDOR);
        HF_zx=zeros(1,Ni_ENDOR);
        HF_zy=zeros(1,Ni_ENDOR);
        NQI_zz=zeros(1,Ni_ENDOR);
        NQI=zeros(Ni_ENDOR,3,3);


        CS_zz=zeros(1,Ni_ENDOR);

        for m=1:Ni_ENDOR

            HF_zz(m)=A(3*m,3);
            if parameters.Bterm==1
                HF_zy(m)=A(3*m,2);
                HF_zx(m)=A(3*m,1);
            end

            if spinSys("Q_used")==true
                NQI_zz(m)=Q(3*m,3);
                NQI(m,:,:)=Q;
            end

            % CS ENDOR
            if spinSys("CS_used")==true
                CS_zz(m)=CS(3*m,3);
            end
        end

        % nuclear dipolar-dipolar coupling ENDOR
        D_zz=zeros(1,Ni_ENDOR);
        if spinSys("D_used")==true
             D_zz=zeros(1,size(D,1)/3);
             for m=1:size(D,1)/3
                 D_zz(m)=D(3*m,3);
             end
         end

        offsets_tmp=kehl_offsets(constants,spinSys,spinOps,paramsENDOR,B,geff,HF_zz,NQI_zz);

        sel_I=parameters.sel_I;
        if spinSys('N_SpinSys')>1
            s=size(offsets_tmp,2)/Ni_ENDOR;
            for n=1:Ni_ENDOR
                offsets(n)=offsets_tmp(sel_I+(n-1)*s);
            end
        else
            s=size(offsets_tmp,2)/Ni_ENDOR;
            offsets=offsets_tmp(((sel_I-1)*s+1):(sel_I*s));
        end

        scalefactor=1;
        geff_sel(or)=geff(1);
        B_sel(or)=B;
        HF_zz_sel(or,:)=HF_zz(:);
        HF_zy_sel(or,:)=HF_zy(:);
        HF_zx_sel(or,:)=HF_zx(:);
        NQI_zz_sel(or,:)=NQI_zz(:);
        NQI_sel(or,:,:,:)=NQI(:,:,:);
        CS_zz_sel(or,:)=CS_zz(:);
        D_zz_sel(or,:)=D_zz(:);
        S_sel(or)=scalefactor;
        offsets_sel(or,:)=offsets(:);
    end

    % build output Map
    EPR=containers.Map;
    EPR("geff_sel")=geff_sel;
    EPR("B_sel")=B_sel;
    EPR("HF_zz_sel")=HF_zz_sel;
    EPR("HF_zy_sel")=HF_zy_sel;
    EPR("HF_zx_sel")=HF_zx_sel;
    EPR("CS_zz_sel")=CS_zz_sel;
    EPR("D_zz_sel")=D_zz_sel;

    EPR("NQI_zz_sel")=NQI_zz_sel;
    EPR("NQI_sel")=NQI_sel;

    EPR("S_sel")=S_sel;
    EPR("EPR_amp")=epr_amp;

    EPR("offsets")=offsets_sel;
end

function scalefactor=kehl_ori_freq_pulsescale(expt,DeltaOm)
    % calculates the scalefactor in dependence of the mw pulse's excitation
    % profile
    %
    % input parameters:
    % expt: the Map containing the experimental parameters
    % DeltaOm: offset of the actual freq from the resonance freq
    %
    % output parameters:
    % scalefactor: scalefactor for the specific orientation
    %
    % February 2024 A. Kehl (akehl@gwdg.de)
    %
    data=expt("pulse");
    %%
    pulse=fopen(data);

    pulse_data=textscan(pulse,'%f %f %f %f');
    fclose('all');
    pulse_x=pulse_data{2};
    dx=1000/(pulse_x(2)-pulse_x(1));

    n=1000;


    pulse_y=pulse_data{3}+1i*pulse_data{4};

    Y=abs(fftshift(fft(pulse_y(1:70),(n+1))));
    y=Y.^3;
    y=y/max(y);
    Y=Y/max(Y);
    if isKey(expt,"3pulses") && expt("3pulses")==true
       Y=y;
    end
%%

    %Delta omega in MHz
    DeltaOM=DeltaOm*2*pi*1e-6;
    binDOM=round(DeltaOM)+(dx/2+1);
    if binDOM>0 && binDOM<(n+1)
        scalefactor=Y(binDOM);
    else
        scalefactor=0;
    end
end

function kehl_ori_freq_grumble(constants,spinSys,spinOps,paramsEPR,paramsENDOR,parameters,expt)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if ~isa(paramsEPR,'containers.Map')
    error('paramsEPR must be a containers.Map object.');
end
if ~isa(paramsENDOR,'containers.Map')
    error('paramsENDOR must be a containers.Map object.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
end


% --- Localised from kehl_line_broaden.m ---
% convolutes the input data with a Lorentian and/or a
% Gaussian line broadening by applying a fft, then convoluting with an
% exponential function and applying a fft again.
%
% input parameters:
% data: data to be convoluted
% parameters: structure, should contain fields ('Lorentian'/'Gaussian' and
%       'lw_L'/'lw_G')
% sw: width of the data points
%
% output parameters:
% data_conv: convoluted data
%
% April 2024 A. Kehl (akehl@gwdg.de)
%

function data_conv=kehl_line_broaden(data,parameters,sw)

    % Check consistency
    kehl_line_broaden_grumble(data,parameters,sw);
    if parameters.Lorentzian==1

        % Lorentian
        Deltaend_L=parameters.lw_L*pi/(sw);
        endintens1=ifft(data(:));
        for i=1:size(data,2)

            % Lorentzian line shape
            endintens(i)=endintens1(i)*exp(-Deltaend_L*i);
        end
        endintens(1)=0.37*endintens(1);
        endamp_2d_L_conv(:)=real(fft(endintens));

        if parameters.Gaussian==1
            % Gaussian
            lw=parameters.lw_G;
            Deltaend_G=(lw*pi/(sw*sqrt(2*log(2))))^2;
            endintens1=ifft(endamp_2d_L_conv(:)+1000);
            for i=1:size(data,2)

                % Gaussian line shape
                endintens(i)=endintens1(i)*exp(-Deltaend_G*i^2/2);
            end
            endintens(1)=0.5*endintens(1);
            data_conv(:)=-abs(fft(endintens));

        else
            data_conv=endamp_2d_L_conv(:);
        end
    elseif parameters.Gaussian==1
        % Gaussian
        lw=parameters.lw_G;
        Deltaend_G=(lw*pi/(sw*sqrt(2*log(2))))^2;
        data=data+1000;
        endintens1=ifft(data(:));

        %size(data,2)
        for i=1:length(data)

            % Gaussian line shape
            endintens(i)=endintens1(i)*exp(-Deltaend_G*i^2/2);
        end
        endintens(1)=0.5*endintens(1);
        data_conv(:)=-abs(fft(endintens(:)));

    else
        data_conv=data;
    end
end

function kehl_line_broaden_grumble(data,parameters,sw)
if (~isnumeric(data))&&(~ischar(data))&&(~isstring(data))
    error('data must be numeric, a character string, or a string scalar.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isnumeric(sw)
    error('sw must be numeric.');
end
end


% --- Localised from kehl_mims_rlx.m ---
% performs the actual ENDOR calculation (with relaxation)
%
% input parameters:
% constants: the Map containing the constants
% spinOps: the Map containing the spin operators
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
% parameters: structure containing simulation parameters
% paramsENDOR: the Map containing the ENDOR parameters
% EPR: the Map with results from EPR calculation
%
% output parameters:
% endor_amp: endor amplitude values
%
% February 2024 A. Kehl (akehl@gwdg.de)
%


function endor_amp=kehl_mims_rlx(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)

    % Check consistency
    kehl_mims_rlx_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR);
    spinOps_D=kehl_diag_ops(spinOps,spinSys('Ni_ENDOR'));

    % get values from Maps
    Sz=spinOps("Sz");
    Sx=spinOps("Sx");
    Sy=spinOps("Sy");
    Iz=spinOps("Iz");
    Ix=spinOps("Ix");
    Iy=spinOps("Iy");

    Sz_D=spinOps_D("Sz");
    Sx_D=spinOps_D("Sx");
    Sy_D=spinOps_D("Sy");
    Iz_D=spinOps_D("Iz");
    Ix_D=spinOps_D("Ix");
    Iy_D=spinOps_D("Iy");

    t=expt("t");
    Nint=100;


    N_spinSys=spinSys("N_SpinSys");

    Ni_ENDOR=spinSys("Ni_ENDOR");

    if N_spinSys>1
       nucs=Ni_ENDOR;
        for i=1:Ni_ENDOR-1
            Ix{i+1}=Ix{1};
            Iy{i+1}=Iy{1};
            Iz{i+1}=Iz{1};
        end
    else
       nucs=1;
    end



    geff_sel=EPR("geff_sel");
    B_sel=EPR("B_sel");
    HF_zz_sel=EPR("HF_zz_sel");
    HF_zy_sel=EPR("HF_zy_sel");
    HF_zx_sel=EPR("HF_zx_sel");

    NQI_zz_sel=EPR("NQI_zz_sel");
    NQI_sel=EPR("NQI_sel");
    CS_zz_sel=EPR("CS_zz_sel");
    D_zz_sel=EPR("D_zz_sel");

    S_sel=EPR("S_sel");


    offsets_sel=EPR("offsets");
    Npts_EN=paramsENDOR("Npts_EN");

    I=spinSys("I");
    v_L=paramsENDOR("v_L");


    endor_amp=zeros(1,Npts_EN);


     if length(B_sel)==0
        % No resonance orientations were found
        return
     end

    % loop to repeat the calculation for every orientation
    parfor j=1:length(B_sel)
        % set parameters for this orientation
        geff=geff_sel(j);
        B=B_sel(j);

        HF_zz=HF_zz_sel(j,:);

        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);
        NQI=zeros(Ni_ENDOR,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);
        CS_zz=CS_zz_sel(j,:)*1e-12;
        D_zz=D_zz_sel(j,:);

        S=S_sel(j);
        offsets=offsets_sel(j,:);

        endor_amp_tmp=zeros(1,Npts_EN);


        if abs(HF_zz(1))>1/spinSys("T2e")*0.1

            % loop over all offsets (spin manifolds)
            for i=1 : size(offsets,2)

                v_off_S=offsets(i);
                off_1=offsets(1);

                rho0=kehl_rho0(constants,paramsENDOR,B,geff,spinOps,spinSys,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

                %electron T1
                RT1e=kehl_relax_t1(rho0,Sx_D,spinSys("T1e"));
                RT1n=zeros(size(RT1e));
                RT2e=kehl_relax_t2(Sx_D,spinSys("T2e"));
                RT2n=zeros(size(RT2e));

                if N_spinSys==1
                    for mm=1:Ni_ENDOR

                       %nuclear T1
                       RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{mm},spinSys("T1n"));
                       RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{mm},spinSys("T2dq"));
                       RT2n=RT2n+kehl_relax_t2(Ix_D{mm},spinSys("T2n"));
                    end
                else
                    RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{1},spinSys("T1n"));
                    RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{1},spinSys("T2dq"));
                    RT2n=RT2n+kehl_relax_t2(Ix_D{1},spinSys("T2n"));
                end
                R=RT1e+RT1n+RT2e+RT2n;


                start_EN=paramsENDOR("start_EN");
                step_EN=paramsENDOR("step_EN");

                oneE=expt("oneE");
                oneN=expt("oneN");

                Hfree_p=2*pi*v_off_S*Sz;
                if spinSys('N_SpinSys')>1
                    if parameters.powder==1
                        mn=round(i/(2*I(1)+1));
                    else
                        mn=i;
                    end
                    % HF
                    Hfree_p=Hfree_p-2*pi*v_L(mn)*Iz{1}+2*pi*v_L(mn)*Iz{1}*CS_zz(mn)+2*pi*HF_zz(mn)*(Sz*Iz{1})+2*pi*HF_zy(mn)*(Sz*Iy{1})+2*pi*HF_zx(mn)*(Sz*Ix{1});
                    % NQI
                    if parameters.Bterm==true
                        Hfree_p=Hfree_p+NQI(mn,1,1)*Ix{1}*Ix{1};
                        Hfree_p=Hfree_p+NQI(mn,1,2)*Ix{1}*Iy{1};
                        Hfree_p=Hfree_p+NQI(mn,1,3)*Ix{1}*Iz{1};
                        Hfree_p=Hfree_p+NQI(mn,2,1)*Iy{1}*Ix{1};
                        Hfree_p=Hfree_p+NQI(mn,2,2)*Iy{1}*Iy{1};
                        Hfree_p=Hfree_p+NQI(mn,2,3)*Iy{1}*Iz{1};
                        Hfree_p=Hfree_p+NQI(mn,3,1)*Iz{1}*Ix{1};
                        Hfree_p=Hfree_p+NQI(mn,3,2)*Iz{1}*Iy{1};
                        Hfree_p=Hfree_p+NQI(mn,3,3)*Iz{1}*Iz{1};
                    else
                        Hfree_p=Hfree_p+pi*NQI_zz(mn)*(3*Iz{1}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                    end
                else
                    for mm=1:Ni_ENDOR
                        % HF
                        Hfree_p=Hfree_p-2*pi*v_L(mm)*Iz{mm}+2*pi*HF_zz(mm)*(Sz*Iz{mm})+2*pi*HF_zy(mm)*(Sz*Iy{mm})+2*pi*HF_zx(mm)*(Sz*Ix{mm})+2*pi*v_L(mm)*CS_zz(mm)*Iz{mm};
                        % NQI
                        if parameters.Bterm==true
                            Hfree_p=Hfree_p+NQI(mm,1,1)*Ix{mm}*Ix{mm};
                            Hfree_p=Hfree_p+NQI(mm,1,2)*Ix{mm}*Iy{mm};
                            Hfree_p=Hfree_p+NQI(mm,1,3)*Ix{mm}*Iz{mm};
                            Hfree_p=Hfree_p+NQI(mm,2,1)*Iy{mm}*Ix{mm};
                            Hfree_p=Hfree_p+NQI(mm,2,2)*Iy{mm}*Iy{mm};
                            Hfree_p=Hfree_p+NQI(mm,2,3)*Iy{mm}*Iz{mm};
                            Hfree_p=Hfree_p+NQI(mm,3,1)*Iz{mm}*Ix{mm};
                            Hfree_p=Hfree_p+NQI(mm,3,2)*Iz{mm}*Iy{mm};
                            Hfree_p=Hfree_p+NQI(mm,3,3)*Iz{mm}*Iz{mm};
                        else
                            Hfree_p=Hfree_p+pi*NQI_zz(mm)*(3*Iz{mm}*Iz{mm}-I(mm)*(I(mm)+1)*eye(size(Hfree_p)));
                        end
                    end

                    if spinSys("D_used")==true
                        for mm=2:(size(D_zz,2)+1)
                            dipC=2*pi*D_zz(mm-1);

                            HD=zeros(size(Hfree_p));

                            HD=HD+Ix{1}*Ix{mm};
                            HD=HD+Ix{1}*Iy{mm};
                            HD=HD+Ix{1}*Iz{mm};
                            HD=HD+Iy{1}*Ix{mm};
                            HD=HD+Iy{1}*Iy{mm};
                            HD=HD+Iy{1}*Iz{mm};
                            HD=HD+Iz{1}*Ix{mm};
                            HD=HD+Iz{1}*Iy{mm};
                            HD=HD+Iz{1}*Iz{mm};

                            HDip=dipC*(3*Iz{1}*Iz{mm}-HD)/2;
                            Hfree_p=Hfree_p+HDip;
                        end
                    end
                end

                % Integration step for the Signal to account for oscillation
                if v_off_S==0
                    t9=1/(off_1*Nint);
                else
                    t9=1/(v_off_S*Nint);
                end

                % loop over rf frequencies (x-axis)
                for a=1:Npts_EN
                % for a=1   %%% for testing

                    % Radiofrequency
                    v_RF=(start_EN+step_EN*(a-1));

                    Hcorr=zeros(size(Hfree_p));
                    HRF=Hfree_p;

                    if parameters.Bterm==false
                        if spinSys('N_SpinSys')>1
                            m=1;
                            Hcorr=Hcorr+2*pi*v_RF*Iz{m};

                            HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Iy{m};
                        else
                            for mm=1:Ni_ENDOR
                                Hcorr=Hcorr+2*pi*v_RF*Iz{mm};

                                HRF=HRF+2*pi*v_RF*Iz{mm}+oneN*Iy{mm};
                            end
                        end
                    end

                    Hfree=Hfree_p+Hcorr;

                    Hnonsel=full(hilb2liouv(sparse(Hfree+oneE*Sx),'comm'));

                    Hfree=full(hilb2liouv(sparse(Hfree),'comm'));


                    if parameters.Bterm==false

                        U5=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*full(hilb2liouv(sparse(HRF),'comm'))),t(5)));

                        U1=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(1)));
                        U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)));
                        U3=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(3)));
                        U4=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(4)));

                        U6=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(6)));
                        U7=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(7)));
                        U8=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)+t(3)/2));
                        U9=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t9));
                    else
                        U5=kehl_rf_bterm_rlx(parameters,expt,v_RF,Hfree_p,Iy,t(5),Ni_ENDOR,N_spinSys,R,spinOps('spin_system'));

                        U1=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(1)));
                        U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)));
                        U3=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(3)));
                        U4=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(4)));

                        U6=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(6)));
                        U7=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(7)));
                        U8=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)+t(3)/2));
                        U9=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t9));

                    end


                    % Evolve the densitymatrix
                    rho=kehl_mat_to_lbra(rho0)';

                    rho=U1*rho;
                    rho=U2*rho;
                    rho=U3*rho;
                    rho=U4*rho;
                    rho=U5*rho;
                    rho=U6*rho;
                    rho=U7*rho;
                    rho=U8*rho;


                    value_Sy=0;
                    for b=1
                       rho=U9*rho;
                       rho_f=kehl_lket_to_mat(rho);
                       value_Sy=value_Sy+(real(trace(rho_f*Sy)));
                    end
                    endor_amp_tmp(a)=endor_amp_tmp(a)+(value_Sy*S/(Nint*size(offsets,2)));
                end

            end
        end
        endor_amp=endor_amp+endor_amp_tmp;
    end
end

function kehl_mims_rlx_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isa(paramsENDOR,'containers.Map')
    error('paramsENDOR must be a containers.Map object.');
end
if ~isa(EPR,'containers.Map')
    error('EPR must be a containers.Map object.');
end
end


% --- Localised from kehl_mims_calc.m ---
% performs the actual ENDOR calculation (no relaxation)
%
% input parameters:
% constants: the Map containing the constants
% spinOps: the Map containing the spin operators
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
% parameters: structure containing simulation parameters
% paramsENDOR: the Map containing the ENDOR parameters
% EPR: the Map with results from EPR calculation
%
% output parameters:
% endor_amp: endor amplitude values
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

% get values from Maps

function endor_amp=kehl_mims_calc(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)

    % Check consistency
    kehl_mims_calc_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR);
    Sz=spinOps("Sz");
    Sx=spinOps("Sx");
    Sy=spinOps("Sy");
    Iz=spinOps("Iz");
    Ix=spinOps("Ix");
    Iy=spinOps("Iy");
    t=expt("t");
    Nint=8;

    N_spinSys=spinSys("N_SpinSys");

    Ni_ENDOR=spinSys("Ni_ENDOR");

    if N_spinSys>1
       nucs=Ni_ENDOR;
        for i=1:Ni_ENDOR-1
            Ix{i+1}=Ix{1};
            Iy{i+1}=Iy{1};
            Iz{i+1}=Iz{1};
        end
    else
       nucs=1;
    end


    geff_sel=EPR("geff_sel");
    B_sel=EPR("B_sel");
    HF_zz_sel=EPR("HF_zz_sel");
    HF_zy_sel=EPR("HF_zy_sel");
    HF_zx_sel=EPR("HF_zx_sel");

    NQI_zz_sel=EPR("NQI_zz_sel");
    NQI_sel=EPR("NQI_sel");
    CS_zz_sel=EPR("CS_zz_sel");
    D_zz_sel=EPR("D_zz_sel");


    S_sel=EPR("S_sel");


    offsets_sel=EPR("offsets");
    Npts_EN=paramsENDOR("Npts_EN");

    I=spinSys("I");
    v_L=paramsENDOR("v_L");


    endor_amp=zeros(1,Npts_EN);

    if length(B_sel)==0
        % No resonance orientations were found
        return
    end

    % loop to repeat the calculation for every orientation
    parfor j=1:length(B_sel)
        % set parameters for this orientation
        geff=geff_sel(j);
        B=B_sel(j);

        HF_zz=HF_zz_sel(j,:);
        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);
        NQI=zeros(Ni_ENDOR,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);

        CS_zz=CS_zz_sel(j,:)*1e-12;
        D_zz=D_zz_sel(j,:);


        S=S_sel(j);
        offsets=offsets_sel(j,:);

        endor_amp_tmp=zeros(1,Npts_EN);

        % loop over all offsets (spin manifolds)
        for i=1 : size(offsets,2)

            v_off_S=offsets(i);
            off_1=offsets(1);

            [rho0]=kehl_rho0(constants,paramsENDOR,B,geff,spinOps,spinSys,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            oneE=expt("oneE");
            oneN=expt("oneN");

            Hfree_p=2*pi*v_off_S*Sz;
            if spinSys('N_SpinSys')>1
                if parameters.powder==1
                    mn=round(i/(2*I(1)+1));
                else
                    mn=i;
                end
                % HF
                Hfree_p=Hfree_p-2*pi*v_L(mn)*Iz{1}+2*pi*v_L(mn)*Iz{1}*CS_zz(mn)+2*pi*HF_zz(mn)*(Sz*Iz{1})+2*pi*HF_zy(mn)*(Sz*Iy{1})+2*pi*HF_zx(mn)*(Sz*Ix{1});
                % NQI
                if parameters.Bterm==true
                    Hfree_p=Hfree_p+NQI(mn,1,1)*Ix{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,1,2)*Ix{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,1,3)*Ix{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(mn,2,1)*Iy{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,2,2)*Iy{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,2,3)*Iy{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(mn,3,1)*Iz{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,3,2)*Iz{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,3,3)*Iz{1}*Iz{1};
                else
                    Hfree_p=Hfree_p+pi*NQI_zz(mn)*(3*Iz{1}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                end
            else
                for mm=1:Ni_ENDOR
                    % HF
                    Hfree_p=Hfree_p-2*pi*v_L(mm)*Iz{mm}+2*pi*HF_zz(mm)*(Sz*Iz{mm})+2*pi*HF_zy(mm)*(Sz*Iy{mm})+2*pi*HF_zx(mm)*(Sz*Ix{mm})+2*pi*v_L(mm)*CS_zz(mm)*Iz{mm};
                    % NQI
                    if parameters.Bterm==true
                        HNQI=zeros(size(Hfree_p));
                        HNQI=HNQI+Ix{mm}*NQI(mm,1,1)*Ix{mm};
                        HNQI=HNQI+Ix{mm}*NQI(mm,1,2)*Iy{mm};
                        HNQI=HNQI+Ix{mm}*NQI(mm,1,3)*Iz{mm};
                        HNQI=HNQI+Iy{mm}*NQI(mm,2,1)*Ix{mm};
                        HNQI=HNQI+Iy{mm}*NQI(mm,2,2)*Iy{mm};
                        HNQI=HNQI+Iy{mm}*NQI(mm,2,3)*Iz{mm};
                        HNQI=HNQI+Iz{mm}*NQI(mm,3,1)*Ix{mm};
                        HNQI=HNQI+Iz{mm}*NQI(mm,3,2)*Iy{mm};
                        HNQI=HNQI+Iz{mm}*NQI(mm,3,3)*Iz{mm};
                        Hfree_p=Hfree_p+HNQI;
                    else
                        Hfree_p=Hfree_p+pi*NQI_zz(mm)*(3*Iz{mm}*Iz{mm}-I(mm)*(I(mm)+1)*eye(size(Hfree_p)));

                    end

                end

                if spinSys("D_used")==true
                    for mm=2:(size(D_zz,2)+1)
                        dipC=2*pi*D_zz(mm-1);

                        HD=zeros(size(Hfree_p));

                        HD=HD+Ix{1}*Ix{mm};
                        HD=HD+Ix{1}*Iy{mm};
                        HD=HD+Ix{1}*Iz{mm};
                        HD=HD+Iy{1}*Ix{mm};
                        HD=HD+Iy{1}*Iy{mm};
                        HD=HD+Iy{1}*Iz{mm};
                        HD=HD+Iz{1}*Ix{mm};
                        HD=HD+Iz{1}*Iy{mm};
                        HD=HD+Iz{1}*Iz{mm};

                        HDip=dipC*(3*Iz{1}*Iz{mm}-HD)/2;

                        Hfree_p=Hfree_p+HDip;
                    end
                end
            end

            % mw pulses

            %Hfree_p +
            Hnonsel_p=oneE*Sx;

            % Integration step for the Signal to account for oscillation
            if v_off_S==0
                t9=1/abs(off_1*Nint);
            else
                t9=1/abs(v_off_S*Nint);
            end


            % Calculate the propagators
            U1_p=full(propagator(spinOps('spin_system'),sparse(Hnonsel_p),t(1)));
            U2_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(2)));
            if t(3)==t(1)
                U3_p=U1_p;
            else
                U3_p=full(propagator(spinOps('spin_system'),sparse(Hnonsel_p),t(3)));
            end

            U4_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(4)));

            if t(4)==t(6)
                U6_p=U4_p;
            else

                U6_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(6)));
            end

            if t(7)==t(1)
                U7_p=U1_p;
            else
                U7_p=full(propagator(spinOps('spin_system'),sparse(Hnonsel_p),t(7)));
            end
            U8_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(2)+t(3)/2));
            U9_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t9));


            % loop over rf frequencies (x-axis)
            for a=1:Npts_EN
            % for a=1   %%% for testing

                % Radiofrequency
                v_RF=(start_EN+step_EN*(a-1));

                Hcorr=zeros(size(Hfree_p));
                HRF=Hfree_p;

                if parameters.Bterm==false
                    if spinSys('N_SpinSys')>1
                        m=1;
                        HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Iy{m};
                    else
                        for mm=1:Ni_ENDOR
                            HRF=HRF+2*pi*v_RF*Iz{mm}+oneN*Iy{mm};
                        end
                    end
                end

                if spinSys('N_SpinSys')>1
                    m=1;
                    Hcorr=Hcorr+2*pi*v_RF*Iz{m};
                else

                    for mm=1:Ni_ENDOR
                        Hcorr=Hcorr+2*pi*v_RF*Iz{mm};

                    end
                end


                if parameters.Bterm==false
                    U5=full(propagator(spinOps('spin_system'),sparse(HRF),t(5)));

                    U1=U1_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(1)));
                    U2=U2_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(2)));
                    if t(3)==t(1)
                        U3=U1;
                    else
                        U3=U3_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(3)));
                    end
                    U4=U4_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(4)));

                    U4=U4_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(4)));

                    if t(4)==t(6)
                        U6=U4;
                    else
                        U6=U6_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(6)));
                    end

                    if t(7)==t(1)
                        U7=U1;
                    else
                        U7=U7_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(7)));
                    end
                    U8=U8_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(2)+t(3)/2));
                    U9=U9_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t9));
                else
                    Hfree=Hfree_p;

                    U5=kehl_rf_bterm(parameters,expt,v_RF,Hfree,Iy,t(5),Ni_ENDOR,N_spinSys,spinOps('spin_system'));

                    U1=U1_p;
                    U2=U2_p;
                    U3=U3_p;
                    U4=U4_p;

                    U6=U6_p;
                    U7=U7_p;
                    U8=U8_p;
                    U9=U9_p;
                end

                % Evolve the densitymatrix
                rho=rho0;
                rho=U1*rho*U1';
                rho=U2*rho*U2';
                rho=U3*rho*U3';
                rho=U4*rho*U4';
                rho=U5*rho*U5';

                if parameters.Bterm==true
                    rho=diag(diag(rho));
                end

                rho=U6*rho*U6';
                rho=U7*rho*U7';
                rho=U8*rho*U8';


                value_Sy=0;
                value_Sy=value_Sy+(real(trace(rho*Sy)));
                for b=1:Nint*10
                   rho=U9*rho*U9';
                   value_Sy=value_Sy+(real(trace(rho*Sy)));
                end

                % endor_amp_tmp(a) = endor_amp_tmp(a)+abs(value_Sy*S/(Nint*size(offsets,2)));

                if spinSys('N_SpinSys')>1
                    if parameters.powder==0
                        s=1;
                    else
                        s=size(offsets,2)/Ni_ENDOR;
                    end
                else
                    s=size(offsets,2);
                end
                endor_amp_tmp(a)=endor_amp_tmp(a)+(value_Sy*S/(Nint*s));

            end
        end
        endor_amp=endor_amp+endor_amp_tmp;
    end
end

function kehl_mims_calc_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isa(paramsENDOR,'containers.Map')
    error('paramsENDOR must be a containers.Map object.');
end
if ~isa(EPR,'containers.Map')
    error('EPR must be a containers.Map object.');
end
end


% --- Localised from kehl_time_rlx.m ---
% performs the actual ENDOR calculation (with relaxation)
% CAUTION: this is a beta version and not fully tested
%
% input parameters:
% constants: the Map containing the constants
% spinOps: the Map containing the spin operators
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
% parameters: structure containing simulation parameters
% paramsENDOR: the Map containing the ENDOR parameters
% EPR: the Map with results from EPR calculation
%
% output parameters:
% endor_amp: endor amplitude values
%
% April 2024 A. Kehl (akehl@gwdg.de)
%

% get values from Maps

function endor_amp=kehl_time_rlx(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)

    % Check consistency
    kehl_time_rlx_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR);
    spinOps_D=kehl_diag_ops(spinOps,spinSys('Ni_ENDOR'));

    Sz=spinOps("Sz");
    Sx=spinOps("Sx");
    Sy=spinOps("Sy");
    Iz=spinOps("Iz");
    Ix=spinOps("Ix");
    Iy=spinOps("Iy");

    Sz_D=spinOps_D("Sz");
    Sx_D=spinOps_D("Sx");
    Sy_D=spinOps_D("Sy");
    Iz_D=spinOps_D("Iz");
    Ix_D=spinOps_D("Ix");
    Iy_D=spinOps_D("Iy");

    t=expt("t");
    Nint=8;

    N_spinSys=spinSys("N_SpinSys");

    Ni_ENDOR=spinSys("Ni_ENDOR");

    if N_spinSys>1
       nucs=Ni_ENDOR;
        for i=1:Ni_ENDOR-1
            Ix{i+1}=Ix{1};
            Iy{i+1}=Iy{1};
            Iz{i+1}=Iz{1};
        end
    else
       nucs=1;
    end

    geff_sel=EPR("geff_sel");
    B_sel=EPR("B_sel");
    HF_zz_sel=EPR("HF_zz_sel");
    HF_zy_sel=EPR("HF_zy_sel");
    HF_zx_sel=EPR("HF_zx_sel");

    NQI_zz_sel=EPR("NQI_zz_sel");
    NQI_sel=EPR("NQI_sel");
    CS_zz_sel=EPR("CS_zz_sel");
    D_zz_sel=EPR("D_zz_sel");

    S_sel=EPR("S_sel");


    offsets_sel=EPR("offsets");
    Npts_EN=paramsENDOR("Npts_EN");

    I=spinSys("I");
    v_L=paramsENDOR("v_L");


    endor_amp=zeros(1,Npts_EN);
    if length(B_sel)==0
        % No resonance orientations were found
        return
    end
    warning('This is a beta version and has not been fully tested.');


    % loop to repeat the calculation for every orientation
    parfor j=1:length(B_sel)
        endor_amp_tmp=zeros(1,Npts_EN);

        % set parameters for this orientation
        geff=geff_sel(j);
        B=B_sel(j);

        HF_zz=HF_zz_sel(j,:);

        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);
        NQI=zeros(Ni_ENDOR,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);
        CS_zz=CS_zz_sel(j,:)*1e-12;
        D_zz=D_zz_sel(j,:);

        S=S_sel(j);
        offsets=offsets_sel(j,:);


        % loop over all offsets (spin manifolds)
        for i=1 : size(offsets,2)

            v_off_S=offsets(i);
            off_1=offsets(1);

            rho0=kehl_rho0(constants,paramsENDOR,B,geff,spinOps,spinSys,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            %electron T1
            RT1e=kehl_relax_t1(rho0,Sx_D,spinSys("T1e"));
            RT1n=zeros(size(RT1e));
            RT2e=kehl_relax_t2(Sx_D,spinSys("T2e"));
            RT2n=zeros(size(RT2e));

            if N_spinSys==1
                for mn=1:Ni_ENDOR

                   %nuclear T1
                   RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{mn},spinSys("T1n"));
                   RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{mn},spinSys("T2dq"));
                   RT2n=RT2n+kehl_relax_t2(Ix_D{mn},spinSys("T2n"));
                end
            else
                RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{1},spinSys("T1n"));
                RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{1},spinSys("T2dq"));
                RT2n=RT2n+kehl_relax_t2(Ix_D{1},spinSys("T2n"));
            end
            R=RT1e+RT1n+RT2e+RT2n;


            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            oneE=expt("oneE");
            oneN=expt("oneN");

            Hfree_p=2*pi*v_off_S*Sz;
            if spinSys('N_SpinSys')>1
                nuc=round(i/(2*I(1)+1)) ;
                % HF
                Hfree_p=Hfree_p-2*pi*v_L(1)*Iz{1}+2*pi*v_L(1)*Iz{1}*CS_zz(nuc)+2*pi*HF_zz(nuc)*(Sz*Iz{1})+2*pi*HF_zy(nuc)*(Sz*Iy{1})+2*pi*HF_zx(nuc)*(Sz*Ix{1});
                % NQI
                if parameters.Bterm==true
                    Hfree_p=Hfree_p+NQI(1,1,1)*Ix{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(1,1,2)*Ix{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(1,1,3)*Ix{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(1,2,1)*Iy{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(1,2,2)*Iy{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(1,2,3)*Iy{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(1,3,1)*Iz{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(1,3,2)*Iz{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(1,3,3)*Iz{1}*Iz{1};
                else
                    Hfree_p=Hfree_p+pi*NQI_zz(1)*(3*Iz{1}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                end
            else
                for mn=1:Ni_ENDOR
                    % HF
                    Hfree_p=Hfree_p-2*pi*v_L(mn)*Iz{mn}+2*pi*HF_zz(mn)*(Sz*Iz{mn})+2*pi*HF_zy(mn)*(Sz*Iy{mn})+2*pi*HF_zx(mn)*(Sz*Ix{mn})+2*pi*v_L(mn)*CS_zz(mn)*Iz{mn};
                    % NQI
                    if parameters.Bterm==true
                        Hfree_p=Hfree_p+NQI(mn,1,1)*Ix{mn}*Ix{mn};
                        Hfree_p=Hfree_p+NQI(mn,1,2)*Ix{mn}*Iy{mn};
                        Hfree_p=Hfree_p+NQI(mn,1,3)*Ix{mn}*Iz{mn};
                        Hfree_p=Hfree_p+NQI(mn,2,1)*Iy{mn}*Ix{mn};
                        Hfree_p=Hfree_p+NQI(mn,2,2)*Iy{mn}*Iy{mn};
                        Hfree_p=Hfree_p+NQI(mn,2,3)*Iy{mn}*Iz{mn};
                        Hfree_p=Hfree_p+NQI(mn,3,1)*Iz{mn}*Ix{mn};
                        Hfree_p=Hfree_p+NQI(mn,3,2)*Iz{mn}*Iy{mn};
                        Hfree_p=Hfree_p+NQI(mn,3,3)*Iz{mn}*Iz{mn};
                    else
                        Hfree_p=Hfree_p+pi*NQI_zz(mn)*(3*Iz{mn}*Iz{mn}-I(mn)*(I(mn)+1)*eye(size(Hfree_p)));
                    end
                end

                if spinSys("D_used")==true
                    for mn=2:(size(D_zz,2)+1)
                        dipC=2*pi*D_zz(mn-1);

                        HD=zeros(size(Hfree_p));

                        HD=HD+Ix{1}*Ix{mn};
                        HD=HD+Ix{1}*Iy{mn};
                        HD=HD+Ix{1}*Iz{mn};
                        HD=HD+Iy{1}*Ix{mn};
                        HD=HD+Iy{1}*Iy{mn};
                        HD=HD+Iy{1}*Iz{mn};
                        HD=HD+Iz{1}*Ix{mn};
                        HD=HD+Iz{1}*Iy{mn};
                        HD=HD+Iz{1}*Iz{mn};

                        HDip=dipC*(3*Iz{1}*Iz{mn}-HD)/2;
                        Hfree_p=Hfree_p+HDip;
                    end
                end
            end

            % Integration step for the Signal to account for oscillation
            if v_off_S==0
                t11=1/(off_1*Nint);
            else
                t11=1/(v_off_S*Nint);
            end

            % Radiofrequency
            v_RF=v_L(1);

                Hcorr=zeros(size(Hfree_p));
                HRF=Hfree_p;

                if parameters.Bterm==false
                    if spinSys('N_SpinSys')>1
                        m=1;
                        Hcorr=Hcorr+2*pi*v_RF*Iz{m};

                        HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Iy{m};
                    else
                        for mn=1:Ni_ENDOR
                            Hcorr=Hcorr+2*pi*v_RF*Iz{mn};

                            HRF=HRF+2*pi*v_RF*Iz{mn}+oneN*Iy{mn};
                        end
                    end
                end

                Hfree=Hfree_p+Hcorr;

                Hnonsel=full(hilb2liouv(sparse(Hfree+oneE*Sx),'comm'));

                Hfree=full(hilb2liouv(sparse(Hfree),'comm'));



                if parameters.Bterm==false

                    U5=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*full(hilb2liouv(sparse(HRF),'comm'))),t(5)));
                    if t(5)==t(7)
                        U7=U5;
                    else
                        U7=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*full(hilb2liouv(sparse(HRF),'comm'))),t(7)));
                    end

                    U1=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(1)));

                    if t(1)==t(3)
                        U3=U1;
                    else
                        U3=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(3)));
                    end

                    if t(1)==t(9)
                        U9=U1;
                    elseif t(3)==t(9)
                        U9=U3;
                    else
                        U9=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(9)));
                    end

                    U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)));
                    U4=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(4)));

                    if t(8)==t(4)
                        U8=U4;
                    else
                        U8=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(8)));
                    end
                    U10=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)+t(3)/2));
                else
                    U5=kehl_rf_bterm_rlx(parameters,expt,v_RF,Hfree_p,Iy,t(5),Ni_ENDOR,N_spinSys,R,spinOps('spin_system'));
                    if t(5)==t(7)
                        U7=U5;
                    else
                        U7=kehl_rf_bterm_rlx(parameters,expt,v_RF,Hfree_p,Iy,t(7),Ni_ENDOR,N_spinSys,R,spinOps('spin_system'));
                    end

                    U1=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(1)));

                    if t(1)==t(3)
                        U3=U1;
                    else
                        U3=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(3)));

                    end

                    if t(1)==t(9)
                        U9=U1;
                    elseif t(3)==t(9)
                        U9=U3;
                    else
                        U9=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(9)));
                    end

                    U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)));
                    U4=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(4)));
                    if t(8)==t(4)
                        U8=U4;
                    else
                        U8=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(8)));
                    end

                    U10=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)+t(3)/2));
                end

                % Evolve the densitymatrix
                rho=kehl_mat_to_lbra(rho0)';

                rho=U1*rho;
                rho=U2*rho;
                rho=U3*rho;
                rho=U4*rho;
                rho=U5*rho;
                rho_t=rho;


            % loop over different separation of the RF pulses (x-axis)
            for a=1:Npts_EN
                t6=t(6)+step_EN*(a-1);
                U6=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t6));

                rho=rho_t;
                rho=U6*rho;
                rho=U7*rho;
                rho=U8*rho;
                rho=U9*rho;
                rho=U9*rho;
                rho=U10*rho;


                value_Sy=0;
                for b=1
                   rho_f=kehl_lket_to_mat(rho);
                   value_Sy=value_Sy+(real(trace(rho_f*Sy)));
                end

                if abs(HF_zz(1))>1/spinSys("T2e")*0.1
                    endor_amp_tmp(a)=endor_amp_tmp(a)+(value_Sy*S/(Nint*size(offsets,2)));
                end

            end
        end
        endor_amp=endor_amp+endor_amp_tmp;
    end

end

function kehl_time_rlx_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isa(paramsENDOR,'containers.Map')
    error('paramsENDOR must be a containers.Map object.');
end
if ~isa(EPR,'containers.Map')
    error('EPR must be a containers.Map object.');
end
end


% --- Localised from kehl_time_calc.m ---
% performs the actual ENDOR calculation (no relaxation)
%
% input parameters:
% constants: the Map containing the constants
% spinOps: the Map containing the spin operators
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
% parameters: structure containing simulation parameters
% paramsENDOR: the Map containing the ENDOR parameters
% EPR: the Map with results from EPR calculation
%
% output parameters:
% endor_amp: endor amplitude values
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

% get values from Maps

function endor_amp=kehl_time_calc(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)

    % Check consistency
    kehl_time_calc_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR);
    Sz=spinOps("Sz");
    Sx=spinOps("Sx");
    Sy=spinOps("Sy");
    Iz=spinOps("Iz");
    Ix=spinOps("Ix");
    Iy=spinOps("Iy");
    t=expt("t");
    Nint=8;

    N_spinSys=spinSys("N_SpinSys");

    Ni_ENDOR=spinSys("Ni_ENDOR");

    if N_spinSys>1
       nucs=Ni_ENDOR;
        for i=1:Ni_ENDOR-1
            Ix{i+1}=Ix{1};
            Iy{i+1}=Iy{1};
            Iz{i+1}=Iz{1};
        end
    else
       nucs=1;
    end


    geff_sel=EPR("geff_sel");
    B_sel=EPR("B_sel");
    HF_zz_sel=EPR("HF_zz_sel");
    HF_zy_sel=EPR("HF_zy_sel");
    HF_zx_sel=EPR("HF_zx_sel");

    NQI_zz_sel=EPR("NQI_zz_sel");
    NQI_sel=EPR("NQI_sel");
    CS_zz_sel=EPR("CS_zz_sel");
    D_zz_sel=EPR("D_zz_sel");

    if isKey(spinSys,"CS")
       [v_cs,d_cs]=eig(spinSys("CS"));
       cs_iso=trace(d_cs)/3*1e-12;
    else
        cs_iso=0;
    end


    S_sel=EPR("S_sel");


    offsets_sel=EPR("offsets");
    Npts_EN=paramsENDOR("Npts_EN");

    I=spinSys("I");
    v_L=paramsENDOR("v_L");


    endor_amp=zeros(1,Npts_EN);

    if length(B_sel)==0
        % No resonance orientations were found
        return
    end


    % loop to repeat the calculation for every orientation
    parfor j=1:length(B_sel)
        % for j=1
        % set parameters for this orientation
        endor_amp_tmp=zeros(1,Npts_EN);

        geff=geff_sel(j);
        B=B_sel(j);

        HF_zz=HF_zz_sel(j,:);
        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);
        NQI=zeros(Ni_ENDOR,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);

        CS_zz=CS_zz_sel(j,:)*1e-12;
        D_zz=D_zz_sel(j,:);


        S=S_sel(j);
        offsets=offsets_sel(j,:);


        % loop over all offsets (spin manifolds)
        for i=1 : size(offsets,2)

            v_off_S=offsets(i);
            off_1=offsets(1);

            [rho0]=kehl_rho0(constants,paramsENDOR,B,geff,spinOps,spinSys,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            oneE=expt("oneE");
            oneN=expt("oneN");

            Hfree_p=2*pi*v_off_S*Sz;
            if spinSys('N_SpinSys')>1
                nuc=round(i/(2*I(1)+1)) ;
                % HF
                Hfree_p=Hfree_p-2*pi*v_L(1)*Iz{1}+2*pi*v_L(1)*Iz{1}*CS_zz(nuc)+2*pi*HF_zz(nuc)*(Sz*Iz{1})+2*pi*HF_zy(nuc)*(Sz*Iy{1})+2*pi*HF_zx(nuc)*(Sz*Ix{1});
                % NQI
                if parameters.Bterm==true
                    Hfree_p=Hfree_p+NQI(1,1,1)*Ix{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(1,1,2)*Ix{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(1,1,3)*Ix{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(1,2,1)*Iy{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(1,2,2)*Iy{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(1,2,3)*Iy{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(1,3,1)*Iz{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(1,3,2)*Iz{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(1,3,3)*Iz{1}*Iz{1};
                else
                    Hfree_p=Hfree_p+pi*NQI_zz(1)*(3*Iz{1}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                end
            else
                for mm=1:Ni_ENDOR
                    % HF
                    Hfree_p=Hfree_p-2*pi*v_L(mm)*Iz{mm}+2*pi*HF_zz(mm)*(Sz*Iz{mm})+2*pi*HF_zy(mm)*(Sz*Iy{mm})+2*pi*HF_zx(mm)*(Sz*Ix{mm})+2*pi*v_L(mm)*CS_zz(mm)*Iz{mm};
                    % NQI
                    if parameters.Bterm==true
                        HNQI=zeros(size(Hfree_p));
                        HNQI=HNQI+Ix{mm}*NQI(mm,1,1)*Ix{mm};
                        HNQI=HNQI+Ix{mm}*NQI(mm,1,2)*Iy{mm};
                        HNQI=HNQI+Ix{mm}*NQI(mm,1,3)*Iz{mm};
                        HNQI=HNQI+Iy{mm}*NQI(mm,2,1)*Ix{mm};
                        HNQI=HNQI+Iy{mm}*NQI(mm,2,2)*Iy{mm};
                        HNQI=HNQI+Iy{mm}*NQI(mm,2,3)*Iz{mm};
                        HNQI=HNQI+Iz{mm}*NQI(mm,3,1)*Ix{mm};
                        HNQI=HNQI+Iz{mm}*NQI(mm,3,2)*Iy{mm};
                        HNQI=HNQI+Iz{mm}*NQI(mm,3,3)*Iz{mm};
                        Hfree_p=Hfree_p+HNQI;
                    else
                        Hfree_p=Hfree_p+pi*NQI_zz(mm)*(3*Iz{mm}*Iz{mm}-I(mm)*(I(mm)+1)*eye(size(Hfree_p)));

                    end

                end

                if spinSys("D_used")==true
                    for mm=2:(size(D_zz,2)+1)
                        dipC=2*pi*D_zz(mm-1);

                        HD=zeros(size(Hfree_p));

                        HD=HD+Ix{1}*Ix{mm};
                        HD=HD+Ix{1}*Iy{mm};
                        HD=HD+Ix{1}*Iz{mm};
                        HD=HD+Iy{1}*Ix{mm};
                        HD=HD+Iy{1}*Iy{mm};
                        HD=HD+Iy{1}*Iz{mm};
                        HD=HD+Iz{1}*Ix{mm};
                        HD=HD+Iz{1}*Iy{mm};
                        HD=HD+Iz{1}*Iz{mm};

                        HDip=dipC*(3*Iz{1}*Iz{mm}-HD)/2;

                        Hfree_p=Hfree_p+HDip;
                    end
                end
            end
                v_RF=v_L(1);
                Hcorr=zeros(size(Hfree_p));
                if spinSys('N_SpinSys')>1
                    m=1;
                    Hcorr=Hcorr+2*pi*v_RF*Iz{m};
                else

                    for mm=1:Ni_ENDOR
                        Hcorr=Hcorr+2*pi*v_RF*Iz{mm};

                    end
                end

                if parameters.Bterm==false
                    Hfree=Hfree_p+Hcorr;
                else
                    Hfree=Hfree_p;
                end

                % mw pulses
                Hnonsel=Hfree+oneE*Sx;

                HRF=Hfree;

                if parameters.Bterm==false
                    if spinSys('N_SpinSys')>1
                        m=1;
                        HRF=HRF+oneN*Iy{m};
                    else
                        for mm=1:Ni_ENDOR
                            HRF=HRF+oneN*Iy{mm};
                        end
                    end
                end


                % Integration step for the Signal to account for oscillation
                if v_off_S==0
                    t11=1/abs(off_1*Nint);
                else
                    t11=1/abs(v_off_S*Nint);
                end


                if parameters.Bterm==false
                    U5=full(propagator(spinOps('spin_system'),sparse(HRF),t(5)));
                    if t(5)==t(7)
                        U7=U5;
                    else
                        U7=full(propagator(spinOps('spin_system'),sparse(HRF),t(7)));
                    end

                else
                    U5=kehl_rf_bterm(parameters,expt,v_RF,Hfree,Iy,t(5),Ni_ENDOR,N_spinSys,spinOps('spin_system'));
                    if t(5)==t(7)
                        U7=U5;
                    else
                        U7=kehl_rf_bterm(parameters,expt,v_RF,Hfree,Iy,t(7),Ni_ENDOR,N_spinSys,spinOps('spin_system'));
                    end
                end

                % Calculate the propagators
                U1=full(propagator(spinOps('spin_system'),sparse(Hnonsel),t(1)));

                if t(1)==t(3)
                    U3=U1;
                else
                    U3=full(propagator(spinOps('spin_system'),sparse(Hnonsel),t(3)));
                end

                if t(1)==t(9)
                    U9=U1;
                elseif t(3)==t(9)
                    U9=U3;
                else
                    U9=full(propagator(spinOps('spin_system'),sparse(Hfree),t(9)));
                end

                U2=full(propagator(spinOps('spin_system'),sparse(Hfree),t(2)));

                U4=full(propagator(spinOps('spin_system'),sparse(Hfree),t(4)));

                U10=full(propagator(spinOps('spin_system'),sparse(Hfree),t(2)+t(3)/2));

                U11=full(propagator(spinOps('spin_system'),sparse(Hfree),t11));


                % Evolve the densitymatrix
                rho=rho0;
                rho=U1*rho*U1';
                rho=U2*rho*U2';
                rho=U3*rho*U3';
                rho=U4*rho*U4';
                rho=U5*rho*U5';

                if parameters.Bterm==true
                    rho=diag(diag(rho));
                end

            rho_t=rho;

             % loop over different separation of the RF pulses (x-axis)
            for a=1:Npts_EN
                t6=t(6)+step_EN*(a-1);
                t8=t(8)-step_EN*(a-1);


                if parameters.Bterm==false
                    U6=full(propagator(spinOps('spin_system'),sparse(Hfree_p+Hcorr),t6));
                else
                    Hfree=Hfree_p;
                    U6=full(propagator(spinOps('spin_system'),sparse(Hfree),t6));
                end
                U8=full(propagator(spinOps('spin_system'),sparse(Hfree),t8));


                rho=rho_t;
                rho=U6*rho*U6';
                rho=U7*rho*U7';

                if parameters.Bterm==true
                    rho=diag(diag(rho));
                end

                rho=U8*rho*U8';
                rho=U9*rho*U9';
                rho=U10*rho*U10';

                value_Sy=0;
                for b=1:Nint
                   rho=U11*rho*U11';
                   value_Sy=value_Sy+(real(trace(rho*Sy)));
                end

                endor_amp_tmp(a)=endor_amp_tmp(a)+(value_Sy*S/(Nint*size(offsets,2)));
            end
        end
        endor_amp=endor_amp+endor_amp_tmp;
    end

end

function kehl_time_calc_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isa(paramsENDOR,'containers.Map')
    error('paramsENDOR must be a containers.Map object.');
end
if ~isa(EPR,'containers.Map')
    error('EPR must be a containers.Map object.');
end
end


% --- Localised from kehl_davies_rlx.m ---
% performs the actual ENDOR calculation (with relaxation)
%
% input parameters:
% constants: the Map containing the constants
% spinOps: the Map containing the spin operators
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
% parameters: structure containing simulation parameters
% paramsENDOR: the Map containing the ENDOR parameters
% EPR: the Map with results from EPR calculation
%
% output parameters:
% endor_amp: endor amplitude values
%
% February 2024 A. Kehl (akehl@gwdg.de)
%


function endor_amp=kehl_davies_rlx(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)

    % Check consistency
    kehl_davies_rlx_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR);
    spinOps_D=kehl_diag_ops(spinOps,spinSys('Ni_ENDOR'));

    % get values from Maps
    Sz=spinOps("Sz");
    Sx=spinOps("Sx");
    Sy=spinOps("Sy");
    Iz=spinOps("Iz");
    Ix=spinOps("Ix");
    Iy=spinOps("Iy");

    Sz_D=spinOps_D("Sz");
    Sx_D=spinOps_D("Sx");
    Sy_D=spinOps_D("Sy");
    Iz_D=spinOps_D("Iz");
    Ix_D=spinOps_D("Ix");
    Iy_D=spinOps_D("Iy");


    t=expt("t");
    Nint=8;


    N_spinSys=spinSys("N_SpinSys");

    Ni_ENDOR=spinSys("Ni_ENDOR");

    if N_spinSys>1
       nucs=Ni_ENDOR;
        for i=1:Ni_ENDOR-1
            Ix{i+1}=Ix{1};
            Iy{i+1}=Iy{1};
            Iz{i+1}=Iz{1};
        end
    else
       nucs=1;
    end


    geff_sel=EPR("geff_sel");
    B_sel=EPR("B_sel");
    HF_zz_sel=EPR("HF_zz_sel");
    HF_zy_sel=EPR("HF_zy_sel");
    HF_zx_sel=EPR("HF_zx_sel");

    NQI_zz_sel=EPR("NQI_zz_sel");
    NQI_sel=EPR("NQI_sel");
    CS_zz_sel=EPR("CS_zz_sel");

    S_sel=EPR("S_sel");


    offsets_sel=EPR("offsets");
    Npts_EN=paramsENDOR("Npts_EN");

    I=spinSys("I");
    v_L=paramsENDOR("v_L");


    endor_amp=zeros(1,Npts_EN);

    if length(B_sel)==0
        % No resonance orientations were found
        return
    end

    % loop to repeat the calculation for every orientation
    parfor j=1:length(B_sel)

        % set parameters for this orientation
        geff=geff_sel(j);
        B=B_sel(j);

        HF_zz=HF_zz_sel(j,:);
        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);
        NQI=zeros(Ni_ENDOR,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);
        CS_zz=CS_zz_sel(j,:)*1e-12;


        S=S_sel(j);
        offsets=offsets_sel(j,:);

        endor_amp_tmp=zeros(1,Npts_EN);

        if abs(HF_zz(1))>1/spinSys("T2e")*0.1


        % loop over all offsets (spin manifolds)
        for i=1 : size(offsets,2)

            v_off_S=offsets(i);
            off_1=offsets(1);

            rho0=kehl_rho0(constants,paramsENDOR,B,geff,spinOps,spinSys,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            %electron T1
            RT1e=kehl_relax_t1(rho0,Sx_D,spinSys("T1e"));
            RT1n=zeros(size(RT1e));
            RT2e=kehl_relax_t2(Sx_D,spinSys("T2e"));
            RT2n=zeros(size(RT2e));

            if N_spinSys==1
                for mm=1:Ni_ENDOR

                   %nuclear T1
                   RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{mm},spinSys("T1n"));
                   RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{mm},spinSys("T2dq"));
                   RT2n=RT2n+kehl_relax_t2(Ix_D{mm},spinSys("T2n"));
                end
            else
                RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{1},spinSys("T1n"));
                RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{1},spinSys("T2dq"));
                RT2n=RT2n+kehl_relax_t2(Ix_D{1},spinSys("T2n"));
            end

            R=RT1e+RT1n+RT2e+RT2n;


            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            prep=expt("prep");
            oneE=expt("oneE");
            oneN=expt("oneN");

            Hfree_p=2*pi*v_off_S*Sz;
            if spinSys('N_SpinSys')>1
                if parameters.powder==1
                    m=round(i/(2*I(1)+1));
                else
                    m=i;
                end
                % HF
                Hfree_p=Hfree_p-2*pi*v_L(m)*Iz{1}+2*pi*v_L(m)*Iz{1}*CS_zz(m)+2*pi*HF_zz(m)*(Sz*Iz{1})+2*pi*HF_zy(m)*(Sz*Iy{1})+2*pi*HF_zx(m)*(Sz*Ix{1});
                % NQI
                if parameters.Bterm==true
                    Hfree_p=Hfree_p+NQI(m,1,1)*Ix{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(m,1,2)*Ix{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(m,1,3)*Ix{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(m,2,1)*Iy{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(m,2,2)*Iy{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(m,2,3)*Iy{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(m,3,1)*Iz{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(m,3,2)*Iz{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(m,3,3)*Iz{1}*Iz{1};
                else
                    Hfree_p=Hfree_p+pi*NQI_zz(m)*(3*Iz{1}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                end
            else
                for mm=1:Ni_ENDOR
                    % HF
                    Hfree_p=Hfree_p-2*pi*v_L(mm)*Iz{mm}+2*pi*HF_zz(mm)*(Sz*Iz{mm})+2*pi*HF_zy(mm)*(Sz*Iy{mm})+2*pi*HF_zx(mm)*(Sz*Ix{mm})+2*pi*v_L(mm)*CS_zz(mm)*Iz{mm};
                    % NQI
                    if parameters.Bterm==true
                        Hfree_p=Hfree_p+NQI(mm,1,1)*Ix{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,1,2)*Ix{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,1,3)*Ix{mm}*Iz{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,1)*Iy{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,2)*Iy{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,3)*Iy{mm}*Iz{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,1)*Iz{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,2)*Iz{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,3)*Iz{mm}*Iz{mm};
                    else
                        Hfree_p=Hfree_p+pi*NQI_zz(mm)*(3*Iz{mm}*Iz{mm}-I(mm)*(I(mm)+1)*eye(size(Hfree_p)));
                    end
                end
            end

            % Integration step for the Signal to account for oscillation
            if v_off_S==0
                t9=1/(off_1*Nint);
            else
                t9=1/(v_off_S*Nint);
            end

            % loop over rf frequencies (x-axis)
            for a=1:Npts_EN
            % for a=1   %%% for testing

                % Radiofrequency
                v_RF=(start_EN+step_EN*(a-1));

                Hcorr=zeros(size(Hfree_p));
                HRF=Hfree_p;

                if parameters.Bterm==false
                    if spinSys('N_SpinSys')>1
                        m=1;
                        Hcorr=Hcorr+2*pi*v_RF*Iz{m};

                        HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Iy{m};
                    else
                        for mm=1:Ni_ENDOR
                            Hcorr=Hcorr+2*pi*v_RF*Iz{mm};

                            HRF=HRF+2*pi*v_RF*Iz{mm}+oneN*Iy{mm};
                        end
                    end
                end

                Hfree=Hfree_p+Hcorr;

                Hprep=full(hilb2liouv(sparse(Hfree+prep*Sx),'comm'));
                Hnonsel=full(hilb2liouv(sparse(Hfree+oneE*Sx),'comm'));

                Hfree=full(hilb2liouv(sparse(Hfree),'comm'));


                if parameters.Bterm==false

                    U3=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*full(hilb2liouv(sparse(HRF),'comm'))),t(3)));

                    U1=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hprep),t(1)));
                    U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)));

                    U4=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(4)));
                    U5=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(5)));
                    U6=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(6)));
                    U7=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(7)));

                    U8=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(6)+t(5)/2));
                    U9=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t9));
                else
                    U3=kehl_rf_bterm_rlx(parameters,expt,v_RF,Hfree_p,Iy,t(3),Ni_ENDOR,N_spinSys,R,spinOps('spin_system'));

                    U1=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hprep),t(1)));
                    U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)));

                    U4=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(4)));
                    U5=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(5)));
                    U6=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(6)));
                    U7=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(7)));

                    U8=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(6)+t(5)/2));
                    U9=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t9));

                end


                % Evolve the densitymatrix

                rho=kehl_mat_to_lbra(rho0)';

                rho=U1*rho;
                rho=U2*rho;
                rho=U3*rho;
                rho=U4*rho;
                rho=U5*rho;
                rho=U6*rho;
                rho=U7*rho;
                rho=U8*rho;


                value_Sy=0;
                for b=1:Nint
                   rho=U9*rho;
                   rho_f=kehl_lket_to_mat(rho);
                   value_Sy=value_Sy+real(trace(rho_f*Sy));
                end

                endor_amp_tmp(a)=endor_amp_tmp(a)+(value_Sy*S/(Nint*size(offsets,2)));
            end
        end
        end
        endor_amp=endor_amp+endor_amp_tmp;
    end
end

function kehl_davies_rlx_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isa(paramsENDOR,'containers.Map')
    error('paramsENDOR must be a containers.Map object.');
end
if ~isa(EPR,'containers.Map')
    error('EPR must be a containers.Map object.');
end
end


% --- Localised from kehl_davies_calc.m ---
% performs the actual ENDOR calculation (no relaxation)
%
% input parameters:
% constants: the Map containing the constants
% spinOps: the Map containing the spin operators
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
% parameters: structure containing simulation parameters
% paramsENDOR: the Map containing the ENDOR parameters
% EPR: the Map with results from EPR calculation
%
% output parameters:
% endor_amp: endor amplitude values
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

% get values from Maps

function endor_amp=kehl_davies_calc(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)

    % Check consistency
    kehl_davies_calc_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR);
    Sz=spinOps("Sz");
    Sx=spinOps("Sx");
    Sy=spinOps("Sy");
    Iz=spinOps("Iz");
    Ix=spinOps("Ix");
    Iy=spinOps("Iy");
    t=expt("t");
    Nint=8;

    N_spinSys=spinSys("N_SpinSys");

    Ni_ENDOR=spinSys("Ni_ENDOR");

    if N_spinSys>1
       nucs=Ni_ENDOR;
        for i=1:Ni_ENDOR-1
            Ix{i+1}=Ix{1};
            Iy{i+1}=Iy{1};
            Iz{i+1}=Iz{1};
        end
    else
       nucs=1;
    end


    geff_sel=EPR("geff_sel");
    B_sel=EPR("B_sel");
    HF_zz_sel=EPR("HF_zz_sel");
    HF_zy_sel=EPR("HF_zy_sel");
    HF_zx_sel=EPR("HF_zx_sel");

    NQI_zz_sel=EPR("NQI_zz_sel");
    NQI_sel=EPR("NQI_sel");
    CS_zz_sel=EPR("CS_zz_sel");
    D_zz_sel=EPR("D_zz_sel");

    S_sel=EPR("S_sel");


    offsets_sel=EPR("offsets");
    Npts_EN=paramsENDOR("Npts_EN");

    I=spinSys("I");
    v_L=paramsENDOR("v_L");


    endor_amp=zeros(1,Npts_EN);

    if length(B_sel)==0
        % No resonance orientations were found
        return
    end


    % loop to repeat the calculation for every orientation
    for j=1:length(B_sel)
        endor_amp_tmp=zeros(1,Npts_EN);

        % set parameters for this orientation
        geff=geff_sel(j);
        B=B_sel(j);

        HF_zz=HF_zz_sel(j,:);
        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);
        NQI=zeros(Ni_ENDOR,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);
        CS_zz=CS_zz_sel(j,:)*1e-12;
        D_zz=D_zz_sel(j,:);



        S=S_sel(j);
        offsets=offsets_sel(j,:);

        % loop over all offsets (spin manifolds)
        for i=1: size(offsets,2)

            v_off_S=offsets(i);
            if offsets(1)~=0
                off_1=offsets(1);
            else
                off_1=-HF_zz(1);
            end

            rho0=kehl_rho0(constants,paramsENDOR,B,geff,spinOps,spinSys,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            prep=expt("prep");
            oneE=expt("oneE");
            oneN=expt("oneN");

            Hfree_p=2*pi*v_off_S*Sz;
            if spinSys('N_SpinSys')>1
                if parameters.powder==1
                    mn=round(i/(2*I(1)+1));
                else
                    mn=i;
                end
                % HF
                Hfree_p=Hfree_p-2*pi*v_L(mn)*Iz{1}+2*pi*v_L(1)*Iz{1}*CS_zz(mn)+2*pi*HF_zz(mn)*(Sz*Iz{1})+2*pi*HF_zy(mn)*(Sz*Iy{1})+2*pi*HF_zx(mn)*(Sz*Ix{1});
                % NQI
                if parameters.Bterm==true
                    Hfree_p=Hfree_p+NQI(mn,1,1)*Ix{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,1,2)*Ix{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,1,3)*Ix{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(mn,2,1)*Iy{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,2,2)*Iy{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,2,3)*Iy{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(mn,3,1)*Iz{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,3,2)*Iz{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,3,3)*Iz{1}*Iz{1};
                else
                    Hfree_p=Hfree_p+pi*NQI_zz(mn)*(3*Iz{1}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                end

            else
                for mm=1:Ni_ENDOR
                    % HF
                    Hfree_p=Hfree_p-2*pi*v_L(mm)*Iz{mm}+2*pi*HF_zz(mm)*(Sz*Iz{mm})+2*pi*HF_zy(mm)*(Sz*Iy{mm})+2*pi*HF_zx(mm)*(Sz*Ix{mm})+2*pi*v_L(mm)*CS_zz(mm)*Iz{mm};
                    % NQI
                    if parameters.Bterm==true
                        Hfree_p=Hfree_p+NQI(mm,1,1)*Ix{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,1,2)*Ix{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,1,3)*Ix{mm}*Iz{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,1)*Iy{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,2)*Iy{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,3)*Iy{mm}*Iz{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,1)*Iz{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,2)*Iz{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,3)*Iz{mm}*Iz{mm};
                    else
                        Hfree_p=Hfree_p+pi*NQI_zz(mm)*(3*Iz{mm}*Iz{mm}-I(mm)*(I(mm)+1)*eye(size(Hfree_p)));
                    end


                end
                if spinSys("D_used")==true
                    for mm=2:(size(D_zz,2)+1)
                        dipC=2*pi*D_zz(mm-1);

                        HD=zeros(size(Hfree_p));

                        HD=HD+Ix{1}*Ix{mm};
                        HD=HD+Ix{1}*Iy{mm};
                        HD=HD+Ix{1}*Iz{mm};
                        HD=HD+Iy{1}*Ix{mm};
                        HD=HD+Iy{1}*Iy{mm};
                        HD=HD+Iy{1}*Iz{mm};
                        HD=HD+Iz{1}*Ix{mm};
                        HD=HD+Iz{1}*Iy{mm};
                        HD=HD+Iz{1}*Iz{mm};

                        HDip=dipC*(3*Iz{1}*Iz{mm}-HD)/2;
                        Hfree_p=Hfree_p+HDip;
                    end
                end
            end

                % mw pulses
                Hprep_p=Hfree_p+prep*Sx;
                Hnonsel_p=Hfree_p+oneE*Sx;

                % Integration step for the Signal to account for oscillation
                if v_off_S==0
                    t9=1/(off_1*Nint);
                else
                    t9=1/(v_off_S*Nint);
                end


            % Calculate the propagators
            U1_p=full(propagator(spinOps('spin_system'),sparse(Hprep_p),t(1)));
            U2_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(2)));

            if t(4)==t(2)
                U4_p=U2_p;
            else
                U4_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(4)));
            end

            U5_p=full(propagator(spinOps('spin_system'),sparse(Hnonsel_p),t(5)));
            U6_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(6)));

            if t(7)==t(1) && prep==oneE
                U7_p=U1_p;
            else
                U7_p=full(propagator(spinOps('spin_system'),sparse(Hnonsel_p),t(7)));
            end

            U8_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(6)+t(5)/2));
            U9_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t9));


            % loop over rf frequencies (x-axis)
            parfor a=1:Npts_EN
            % for a=1   %%% for testing

                % Radiofrequency
                v_RF=(start_EN+step_EN*(a-1));

                Hcorr=zeros(size(Hfree_p));
                HRF=Hfree_p;

                if parameters.Bterm==false
                    if spinSys('N_SpinSys')>1
                        m=1;
                        HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Iy{m};
                    else
                        for mm=1:Ni_ENDOR
                            HRF=HRF+2*pi*v_RF*Iz{mm}+oneN*Iy{mm};
                        end
                    end
                end

                if spinSys('N_SpinSys')>1
                    m=1;
                    Hcorr=Hcorr+2*pi*v_RF*Iz{m};
                else

                    for mm=1:Ni_ENDOR
                        Hcorr=Hcorr+2*pi*v_RF*Iz{mm};
                    end
                end


                if parameters.Bterm==false
                    U3=full(propagator(spinOps('spin_system'),sparse(HRF),t(3)));

                    U1=U1_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(1)));
                    U2=U2_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(2)));

                    if t(4)==t(2)
                        U4=U2;
                    else
                            U4=U4_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(4)));
                    end

                    U5=U5_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(5)));
                    U6=U6_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(6)));

                    if t(7)==t(1) && prep==oneE
                        U7=U1;
                    else
                        U7=U7_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(7)));
                    end


                    U8=U8_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(6)+t(5)/2));
                    U9=U9_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t9));
                else
                    Hfree=Hfree_p;

                    U3=kehl_rf_bterm(parameters,expt,v_RF,Hfree,Iy,t(3),Ni_ENDOR,N_spinSys,spinOps('spin_system'));

                    U1=U1_p;
                    U2=U2_p;
                    U4=U4_p;
                    U5=U5_p;
                    U6=U6_p;
                    U7=U7_p;
                    U8=U8_p;
                    U9=U9_p;
                end

                % Evolve the densitymatrix
                rho=rho0;
                rho=U1*rho*U1';
                rho=U2*rho*U2';
                rho=U3*rho*U3';

                if parameters.Bterm==true
                    rho=diag(diag(rho));
                end

                rho=U4*rho*U4';
                rho=U5*rho*U5';
                rho=U6*rho*U6';
                rho=U7*rho*U7';
                rho=U8*rho*U8';


                value_Sy=0;
                for b=1:Nint
                   rho=U9*rho*U9';
                   value_Sy=value_Sy+real(trace(rho*Sy));
                end
                if spinSys('N_SpinSys')>1
                    if parameters.powder==0
                        s=1;
                    else
                        s=size(offsets,2)/Ni_ENDOR;
                    end
                else
                    s=size(offsets,2);
                end
                endor_amp_tmp(a)=endor_amp_tmp(a)+(value_Sy*S/(Nint*s));
            end
        end

        endor_amp=endor_amp+endor_amp_tmp;
    end

end

function kehl_davies_calc_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isa(paramsENDOR,'containers.Map')
    error('paramsENDOR must be a containers.Map object.');
end
if ~isa(EPR,'containers.Map')
    error('EPR must be a containers.Map object.');
end
end


% --- Localised from kehl_spinlock_rlx.m ---
% performs the actual ENDOR calculation (with relaxation)
%
% input parameters:
% constants: the Map containing the constants
% spinOps: the Map containing the spin operators
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
% parameters: structure containing simulation parameters
% paramsENDOR: the Map containing the ENDOR parameters
% EPR: the Map with results from EPR calculation
%
% output parameters:
% endor_amp: endor amplitude values
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

function endor_amp=kehl_spinlock_rlx(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)

    % Check consistency
    kehl_spinlock_rlx_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR);
    spinOps_D=kehl_diag_ops(spinOps,spinSys('Ni_ENDOR'));


    % get values from Maps
    Sz=spinOps("Sz");
    Sx=spinOps("Sx");
    Sy=spinOps("Sy");
    Iz=spinOps("Iz");
    Ix=spinOps("Ix");
    Iy=spinOps("Iy");


    Sz_D=spinOps_D("Sz");
    Sx_D=spinOps_D("Sx");
    Sy_D=spinOps_D("Sy");
    Iz_D=spinOps_D("Iz");
    Ix_D=spinOps_D("Ix");
    Iy_D=spinOps_D("Iy");


    t=expt("t");
    Nint=8;


    N_spinSys=spinSys("N_SpinSys");

    Ni_ENDOR=spinSys("Ni_ENDOR");

    if N_spinSys>1
       nucs=Ni_ENDOR;
        for i=1:Ni_ENDOR-1
            Ix{i+1}=Ix{1};
            Iy{i+1}=Iy{1};
            Iz{i+1}=Iz{1};
        end
    else
       nucs=1;
    end



    geff_sel=EPR("geff_sel");
    B_sel=EPR("B_sel");
    HF_zz_sel=EPR("HF_zz_sel");
    HF_zy_sel=EPR("HF_zy_sel");
    HF_zx_sel=EPR("HF_zx_sel");

    NQI_zz_sel=EPR("NQI_zz_sel");
    NQI_sel=EPR("NQI_sel");

    CS_zz_sel=EPR("CS_zz_sel");

    S_sel=EPR("S_sel");


    offsets_sel=EPR("offsets");
    Npts_EN=paramsENDOR("Npts_EN");

    I=spinSys("I");
    v_L=paramsENDOR("v_L");


    endor_amp=zeros(1,Npts_EN);

    if length(B_sel)==0
        % No resonance orientations were found
        return
    end

    % loop to repeat the calculation for every orientation
    parfor j=1:length(B_sel)
        % set parameters for this orientation
        geff=geff_sel(j);
        B=B_sel(j);

        HF_zz=HF_zz_sel(j,:);
        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);
        NQI=zeros(Ni_ENDOR,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);

        CS_zz=CS_zz_sel(j,:)*1e-12;


        S=S_sel(j);
        offsets=offsets_sel(j,:);

        endor_amp_tmp=zeros(1,Npts_EN);

        if abs(HF_zz(1))>1/spinSys("T2e")*0.1


        % loop over all offsets (spin manifolds)
        for i=1 : size(offsets,2)

            v_off_S=offsets(i);
            off_1=offsets(1);

            rho0=kehl_rho0(constants,paramsENDOR,B,geff,spinOps,spinSys,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            %electron T1
            RT1e=kehl_relax_t1(rho0,Sx_D,spinSys("T1e"));
            RT1n=zeros(size(RT1e));
            RT2e=kehl_relax_t2(Sx_D,spinSys("T2e"));
            RT2n=zeros(size(RT2e));

            if N_spinSys==1
                for mm=1:Ni_ENDOR

                   %nuclear T1
                   RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{mm},spinSys("T1n"));
                   RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{mm},spinSys("T2dq"));
                   RT2n=RT2n+kehl_relax_t2(Ix_D{mm},spinSys("T2n"));
                end
            else
                RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{1},spinSys("T1n"));
                RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{1},spinSys("T2dq"));
                RT2n=RT2n+kehl_relax_t2(Ix_D{1},spinSys("T2n"));
            end

            R=RT1e+RT1n+RT2e+RT2n;

            % for tilted frame Rho relaxation times

            %electron T1
            RT1eRho=kehl_relax_t1(rho0,Sx_D,spinSys("T1eR"));
            RT1nRho=zeros(size(RT1eRho));

            % electron T2
            RT2eRho=kehl_relax_t2(Sx_D,spinSys("T2eR"));
            RT2nRho=zeros(size(RT2eRho));

            if N_spinSys==1
                for mm=1:Ni_ENDOR

                   %nuclear T1
                   RT1nRho=RT1eRho+kehl_relax_t1(rho0,Ix_D{mm},spinSys("T1nR"));

                   % double and zero quantum T2
                   RT2eRho=RT2eRho+kehl_relax_t2(Sx_D*Ix_D{mm},spinSys("T2dqR"));

                   % nuclear T2
                   RT2nRho=RT2nRho+kehl_relax_t2(Ix_D{mm},spinSys("T2nR"));
                end
            else

                %nuclear T1
                RT1nRho=RT1eRho+kehl_relax_t1(rho0,Ix_D{1},spinSys("T1nR"));

                % double and zero quantum T2
                RT2eRho=RT2eRho+kehl_relax_t2(Sx_D*Ix_D{1},spinSys("T2dqR"));

                % nuclear T2
                RT2nRho=RT2nRho+kehl_relax_t2(Ix_D{1},spinSys("T2nR"));
            end

            % full relaxation superoperator
            RRho=RT1eRho+RT1nRho+RT2eRho+RT2nRho;


            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            sl=expt("SL");
            oneE=expt("oneE");
            oneN=expt("oneN");

            Hfree_p=2*pi*v_off_S*Sz;
            if spinSys('N_SpinSys')>1
                if parameters.powder==1
                    mn=round(i/(2*I(1)+1));
                else
                    mn=i;
                end
                % HF
                Hfree_p=Hfree_p-2*pi*v_L(mn)*Iz{1}+2*pi*v_L(mn)*Iz{1}*CS_zz(mn)+2*pi*HF_zz(mn)*(Sz*Iz{1})+2*pi*HF_zy(mn)*(Sz*Iy{1})+2*pi*HF_zx(mn)*(Sz*Ix{1});
                % NQI
                if parameters.Bterm==true
                    Hfree_p=Hfree_p+NQI(mn,1,1)*Ix{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,1,2)*Ix{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,1,3)*Ix{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(mn,2,1)*Iy{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,2,2)*Iy{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,2,3)*Iy{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(mn,3,1)*Iz{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,3,2)*Iz{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,3,3)*Iz{1}*Iz{1};
                else
                    Hfree_p=Hfree_p+pi*NQI_zz(mn)*(3*Iz{1}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                end
            else
                for mm=1:Ni_ENDOR
                    % HF
                    Hfree_p=Hfree_p-2*pi*v_L(mm)*Iz{mm}+2*pi*HF_zz(mm)*(Sz*Iz{mm})+2*pi*HF_zy(mm)*(Sz*Iy{mm})+2*pi*HF_zx(mm)*(Sz*Ix{mm})+2*pi*v_L(mm)*CS_zz(mm)*Iz{mm};
                    % NQI
                    if parameters.Bterm==true
                        Hfree_p=Hfree_p+NQI(mm,1,1)*Ix{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,1,2)*Ix{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,1,3)*Ix{mm}*Iz{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,1)*Iy{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,2)*Iy{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,3)*Iy{mm}*Iz{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,1)*Iz{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,2)*Iz{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,3)*Iz{mm}*Iz{mm};
                    else
                        Hfree_p=Hfree_p+pi*NQI_zz(mm)*(3*Iz{mm}*Iz{mm}-I(mm)*(I(mm)+1)*eye(size(Hfree_p)));
                    end
                end
            end

            % Integration step for the Signal to account for oscillation
            if v_off_S==0
                t9=1/(off_1*Nint);
            else
                t9=1/(v_off_S*Nint);
            end


            % loop over rf frequencies (x-axis)
            for a=1:Npts_EN
            % for a=1   %%% for testing

                % Radiofrequency
                v_RF=(start_EN+step_EN*(a-1));

                Hcorr=zeros(size(Hfree_p));
                HRF=Hfree_p;

                if parameters.Bterm==false
                    if spinSys('N_SpinSys')>1
                        m=1;
                        Hcorr=Hcorr+2*pi*v_RF*Iz{m};

                        HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Iy{m};
                    else
                        for mm=1:Ni_ENDOR
                            Hcorr=Hcorr+2*pi*v_RF*Iz{mm};

                            HRF=HRF+2*pi*v_RF*Iz{mm}+oneN*Iy{mm};
                        end
                    end
                end

                Hfree=Hfree_p+Hcorr;

                HSL=Hfree+sl*Sy;

                [V,D]=eig(HSL);
                if ~issorted(diag(D))
                  [d,inds]=sort(diag(D));
                  D=diag(d);
                  V=V(:, inds);
                end
                %eigenvalues are order from negative large to positive large!
                W=V^(-1);

                HSL_t=W*HSL*W^(-1);


                HSL_t=full(hilb2liouv(sparse(HSL_t),'comm'));
                Hnonsel=full(hilb2liouv(sparse(Hfree+oneE*Sx),'comm'));

                Hfree=full(hilb2liouv(sparse(Hfree),'comm'));


                if parameters.Bterm==false

                    U3=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*full(hilb2liouv(sparse(HRF),'comm'))),t(3)));

                    U1=full(propagator(spinOps('spin_system'),1i*sparse(RRho-1i*HSL_t),t(1)));
                    U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)));

                    U4=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(4)));
                    U5=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(5)));
                    U6=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(6)));
                    U7=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(7)));

                    U8=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(6)+t(5)/2));
                    U9=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t9));
                else
                    U3=kehl_rf_bterm_rlx(parameters,expt,v_RF,Hfree_p,Iy,t(3),Ni_ENDOR,N_spinSys,R,spinOps('spin_system'));

                    U1=full(propagator(spinOps('spin_system'),1i*sparse(RRho-1i*HSL_t),t(1)));
                    U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)));

                    U4=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(4)));
                    U5=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(5)));
                    U6=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(6)));
                    U7=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(7)));

                    U8=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(6)+t(5)/2));
                    U9=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t9));

                end

                % Evolve the densitymatrix

                rho=kehl_mat_to_lbra(rho0)';

                rho_t=W*kehl_lket_to_mat(rho)*W^(-1);
                rho_t=kehl_mat_to_lbra(rho_t)';
                rho_t=U1*rho_t;
                rho=V*kehl_lket_to_mat(rho_t)*V^(-1);
                rho=kehl_mat_to_lbra(rho)';

                % rho = U1*rho;
                rho=U2*rho;
                rho=U3*rho;
                rho=U4*rho;
                rho=U5*rho;
                rho=U6*rho;
                rho=U7*rho;
                rho=U8*rho;


                value_Sy=0;
                for b=1:Nint
                   rho=U9*rho;
                   rho_f=kehl_lket_to_mat(rho);
                   value_Sy=value_Sy+real(trace(rho_f*Sy));
                end

                endor_amp_tmp(a)=endor_amp_tmp(a)+(value_Sy*S/(Nint*size(offsets,2)));

            end
        end
        end
        endor_amp=endor_amp+endor_amp_tmp;
    end
end

function kehl_spinlock_rlx_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isa(paramsENDOR,'containers.Map')
    error('paramsENDOR must be a containers.Map object.');
end
if ~isa(EPR,'containers.Map')
    error('EPR must be a containers.Map object.');
end
end


% --- Localised from kehl_spinlock_calc.m ---
% performs the actual ENDOR calculation (no relaxation)
%
% input parameters:
% constants: the Map containing the constants
% spinOps: the Map containing the spin operators
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
% parameters: structure containing simulation parameters
% paramsENDOR: the Map containing the ENDOR parameters
% EPR: the Map with results from EPR calculation
%
% output parameters:
% endor_amp: endor amplitude values
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

% get values from Maps

function endor_amp=kehl_spinlock_calc(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)

    % Check consistency
    kehl_spinlock_calc_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR);
    Sz=spinOps("Sz");
    Sx=spinOps("Sx");
    Sy=spinOps("Sy");
    Iz=spinOps("Iz");
    Ix=spinOps("Ix");
    Iy=spinOps("Iy");
    t=expt("t");
    Nint=8;

    N_spinSys=spinSys("N_SpinSys");

    Ni_ENDOR=spinSys("Ni_ENDOR");

    if N_spinSys>1
       nucs=Ni_ENDOR;
        for i=1:Ni_ENDOR-1
            Ix{i+1}=Ix{1};
            Iy{i+1}=Iy{1};
            Iz{i+1}=Iz{1};
        end
    else
       nucs=1;
    end


    geff_sel=EPR("geff_sel");
    B_sel=EPR("B_sel");
    HF_zz_sel=EPR("HF_zz_sel");
    HF_zy_sel=EPR("HF_zy_sel");
    HF_zx_sel=EPR("HF_zx_sel");

    NQI_zz_sel=EPR("NQI_zz_sel");
    NQI_sel=EPR("NQI_sel");

    CS_zz_sel=EPR("CS_zz_sel");
    D_zz_sel=EPR("D_zz_sel");

    S_sel=EPR("S_sel");


    offsets_sel=EPR("offsets");
    Npts_EN=paramsENDOR("Npts_EN");

    I=spinSys("I");
    v_L=paramsENDOR("v_L");


    endor_amp=zeros(1,Npts_EN);
    if length(B_sel)==0
        % No resonance orientations were found
        return
    end


    % loop to repeat the calculation for every orientation
    parfor j=1:length(B_sel)

        % set parameters for this orientation
        geff=geff_sel(j);
        B=B_sel(j);

        HF_zz=HF_zz_sel(j,:);
        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);

        NQI=zeros(Ni_ENDOR,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);

        CS_zz=CS_zz_sel(j,:)*1e-12;
        D_zz=D_zz_sel(j,:);


        S=S_sel(j);
        offsets=offsets_sel(j,:);

        endor_amp_tmp=zeros(1,Npts_EN);


        % loop over all offsets (spin manifolds)
        for i=1 : size(offsets,2)

            v_off_S=offsets(i);
            off_1=offsets(1);

            rho0=kehl_rho0(constants,paramsENDOR,B,geff,spinOps,spinSys,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            sl=expt("SL");
            oneE=expt("oneE");
            oneN=expt("oneN");

            Hfree_p=2*pi*v_off_S*Sz;
            if spinSys('N_SpinSys')>1
                if parameters.powder==1
                    mn=round(i/(2*I(1)+1));
                else
                    mn=i;
                end
                % HF
                Hfree_p=Hfree_p-2*pi*v_L(mn)*Iz{1}+2*pi*v_L(mn)*Iz{1}*CS_zz(mn)+2*pi*HF_zz(mn)*(Sz*Iz{1})+2*pi*HF_zy(mn)*(Sz*Iy{1})+2*pi*HF_zx(mn)*(Sz*Ix{1});
                % NQI
                if parameters.Bterm==true
                    Hfree_p=Hfree_p+NQI(mn,1,1)*Ix{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,1,2)*Ix{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,1,3)*Ix{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(mn,2,1)*Iy{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,2,2)*Iy{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,2,3)*Iy{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(mn,3,1)*Iz{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,3,2)*Iz{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,3,3)*Iz{1}*Iz{1};
                else
                    Hfree_p=Hfree_p+pi*NQI_zz(mn)*(3*Iz{1}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                end

            else
                for mm=1:Ni_ENDOR
                    % HF
                    Hfree_p=Hfree_p-2*pi*v_L(mm)*Iz{mm}+2*pi*HF_zz(mm)*(Sz*Iz{mm})+2*pi*HF_zy(mm)*(Sz*Iy{mm})+2*pi*HF_zx(mm)*(Sz*Ix{mm})+2*pi*v_L(mm)*CS_zz(mm)*Iz{mm};
                    % NQI
                    if parameters.Bterm==true
                        Hfree_p=Hfree_p+NQI(mm,1,1)*Ix{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,1,2)*Ix{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,1,3)*Ix{mm}*Iz{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,1)*Iy{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,2)*Iy{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,3)*Iy{mm}*Iz{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,1)*Iz{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,2)*Iz{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,3)*Iz{mm}*Iz{mm};
                    else
                        Hfree_p=Hfree_p+pi*NQI_zz(mm)*(3*Iz{mm}*Iz{mm}-I(mm)*(I(mm)+1)*eye(size(Hfree_p)));
                    end
                end

                if spinSys("D_used")==true
                    for mm=2:(size(D_zz,2)+1)
                        dipC=2*pi*D_zz(mm-1);

                        HD=zeros(size(Hfree_p));

                        HD=HD+Ix{1}*Ix{mm};
                        HD=HD+Ix{1}*Iy{mm};
                        HD=HD+Ix{1}*Iz{mm};
                        HD=HD+Iy{1}*Ix{mm};
                        HD=HD+Iy{1}*Iy{mm};
                        HD=HD+Iy{1}*Iz{mm};
                        HD=HD+Iz{1}*Ix{mm};
                        HD=HD+Iz{1}*Iy{mm};
                        HD=HD+Iz{1}*Iz{mm};

                        HDip=dipC*(3*Iz{1}*Iz{mm}-HD)/2;
                        Hfree_p=Hfree_p+HDip;
                    end
                end
            end

            % mw pulses
            HSL_p=Hfree_p+sl*Sy;
            Hnonsel_p=Hfree_p+oneE*Sx;

            % Integration step for the Signal to account for oscillation
            if v_off_S==0
                t9=1/(off_1*Nint);
            else
                t9=1/(v_off_S*Nint);
            end


            % Calculate the propagators
            U1_p=full(propagator(spinOps('spin_system'),sparse(HSL_p),t(1)));
            U2_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(2)));

            U4_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(4)));
            U5_p=full(propagator(spinOps('spin_system'),sparse(Hnonsel_p),t(5)));
            U6_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(6)));
            U7_p=full(propagator(spinOps('spin_system'),sparse(Hnonsel_p),t(7)));
            U8_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(6)+t(5)/2));
            U9_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t9));


            % loop over rf frequencies (x-axis)
            for a=1:Npts_EN
            % for a=1   %%% for testing

                % Radiofrequency
                v_RF=(start_EN+step_EN*(a-1));

                Hcorr=zeros(size(Hfree_p));
                HRF=Hfree_p;

                if parameters.Bterm==false
                    if spinSys('N_SpinSys')>1
                        m=1;
                        HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Iy{m};
                    else
                        for mm=1:Ni_ENDOR
                            HRF=HRF+2*pi*v_RF*Iz{mm}+oneN*Iy{mm};
                        end
                    end
                end

                if spinSys('N_SpinSys')>1
                    m=1;
                    Hcorr=Hcorr+2*pi*v_RF*Iz{m};
                else

                    for mm=1:Ni_ENDOR
                        Hcorr=Hcorr+2*pi*v_RF*Iz{mm};
                    end
                end


                if parameters.Bterm==false
                    U3=full(propagator(spinOps('spin_system'),sparse(HRF),t(3)));

                    U1=U1_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(1)));
                    U2=U2_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(2)));

                    U4=U4_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(4)));
                    U5=U5_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(5)));
                    U6=U6_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(6)));
                    U7=U7_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(7)));

                    U8=U8_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(6)+t(5)/2));
                    U9=U9_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t9));
                else
                    Hfree=Hfree_p;

                    U3=kehl_rf_bterm(parameters,expt,v_RF,Hfree,Iy,t(3),Ni_ENDOR,N_spinSys,spinOps('spin_system'));

                    U1=U1_p;
                    U2=U2_p;
                    U4=U4_p;
                    U5=U5_p;
                    U6=U6_p;
                    U7=U7_p;
                    U8=U8_p;
                    U9=U9_p;
                end

                HSL=HSL_p;
                if spinSys('N_SpinSys')>1
                    m=1;
                    HSL=HSL+2*pi*v_RF*Iz{m};
                else
                    for mm=1:Ni_ENDOR
                        HSL=HSL+2*pi*v_RF*Iz{mm};
                    end
                end

                [V,D]=eig(HSL);
                if ~issorted(diag(D))
                  [d,inds]=sort(diag(D));
                  D=diag(d);
                  V=V(:, inds);
                end
                %eigenvalues are order from negative large to positive large!
                W=V^(-1);



                % Evolve the densitymatrix
                rho=rho0;
                rho=U1*rho*U1';
                rho_t=W*rho*W^(-1);
                rho_t=diag(diag(rho_t));
                rho=V*rho_t*V^(-1);
                rho=U2*rho*U2';
                rho=U3*rho*U3';

                if parameters.Bterm==true
                    rho=diag(diag(rho));
                end

                rho=U4*rho*U4';
                rho=U5*rho*U5';
                rho=U6*rho*U6';
                rho=U7*rho*U7';
                rho=U8*rho*U8';


                value_Sy=0;
                for b=1:Nint
                   rho=U9*rho*U9';
                   value_Sy=value_Sy+real(trace(rho*Sy));
                end

                if spinSys('N_SpinSys')>1
                    if parameters.powder==0
                        s=1;
                    else
                        s=size(offsets,2)/Ni_ENDOR;
                    end
                else
                    s=size(offsets,2);
                end

                endor_amp_tmp(a)=endor_amp_tmp(a)+(value_Sy*S/(Nint*s));

            end
        end
        endor_amp=endor_amp+endor_amp_tmp;
    end
end

function kehl_spinlock_calc_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isa(paramsENDOR,'containers.Map')
    error('paramsENDOR must be a containers.Map object.');
end
if ~isa(EPR,'containers.Map')
    error('EPR must be a containers.Map object.');
end
end


% --- Localised from kehl_tensor_rlx.m ---
% performs the actual ENDOR calculation (with relaxation)
%
% input parameters:
% constants: the Map containing the constants
% spinOps: the Map containing the spin operators
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
% parameters: structure containing simulation parameters
% paramsENDOR: the Map containing the ENDOR parameters
% EPR: the Map with results from EPR calculation
%
% output parameters:
% endor_amp: endor amplitude values
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

% redefine, one spin system is not possible

function endor_amp=kehl_tensor_rlx(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)

    % Check consistency
    kehl_tensor_rlx_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR);
    spinSys('N_SpinSys')=spinSys('Ni_ENDOR');
    spinOps=kehl_spin_ops(spinSys("S"),spinSys("I"),spinSys("Ni_ENDOR"),spinSys("Ni_ENDOR"));

    spinOps_D=kehl_diag_ops(spinOps,spinSys('Ni_ENDOR'));

    % get values from Maps
    Sz=spinOps("Sz");
    Sx=spinOps("Sx");
    Sy=spinOps("Sy");
    Iz=spinOps("Iz");
    Ix=spinOps("Ix");
    Iy=spinOps("Iy");

    Sz_D=spinOps_D("Sz");
    Sx_D=spinOps_D("Sx");
    Sy_D=spinOps_D("Sy");
    Iz_D=spinOps_D("Iz");
    Ix_D=spinOps_D("Ix");
    Iy_D=spinOps_D("Iy");

    t=expt("t");
    Nint=8;


    N_spinSys=spinSys("N_SpinSys");

    Ni_ENDOR=spinSys("Ni_ENDOR");

    geff_sel=EPR("geff_sel");
    B_sel=EPR("B_sel");
    HF_zz_sel=EPR("HF_zz_sel");
    HF_zy_sel=EPR("HF_zy_sel");
    HF_zx_sel=EPR("HF_zx_sel");

    NQI_zz_sel=EPR("NQI_zz_sel");
    NQI_sel=EPR("NQI_sel");
    CS_zz_sel=EPR("CS_zz_sel");


    S_sel=EPR("S_sel");


    offsets_sel=EPR("offsets");
    Npts_EN=paramsENDOR("Npts_EN");

    I=spinSys("I");
    v_L=paramsENDOR("v_L");


    endor_amp=zeros(1,Npts_EN);

    if length(B_sel)==0
        % No resonance orientations were found
        return
    end

    % loop to repeat the calculation for every orientation
    parfor j=1:length(B_sel)

        % set parameters for this orientation
        geff=geff_sel(j);
        B=B_sel(j);
        const_R=1/size(Sz,2)*constants('GE')*B/(2*pi*constants('K_B')*parameters.T);


        HF_zz=HF_zz_sel(j,:);
        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);
        NQI=zeros(Ni_ENDOR,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);
        CS_zz=CS_zz_sel(j,:)*1e-12;

        S=S_sel(j);
        offsets=offsets_sel(j,:);

        endor_amp_tmp=zeros(1,Npts_EN);

        % loop over nuclei
        for i=1 : Ni_ENDOR
            I=spinSys('I');
            mI=-I(i):1:I(i);

            for jj=mI

                v_off_S=jj*HF_zz(i);
                off_1=-I(i)*HF_zz(i);

                [rho0]=2*Sz*Iz{1};

            RT2e=kehl_relax_t2(Sx_D,spinSys("T2e"));
            RT2n=zeros(size(RT2e));

            if N_spinSys==1
                for mm=1:Ni_ENDOR
                   RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{mm},spinSys("T2dq"));
                   RT2n=RT2n+kehl_relax_t2(Ix_D{mm},spinSys("T2n"));
                end
            else
                RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{1},spinSys("T2dq"));
                RT2n=RT2n+kehl_relax_t2(Ix_D{1},spinSys("T2n"));
            end
            R=RT2e+RT2n;
            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");


            oneN=expt("oneN");

            Hfree_p=2*pi*v_off_S*Sz;
            if spinSys('N_SpinSys')>1
                nuc=i;
                m=1;
                % HF
                Hfree_p=Hfree_p-2*pi*v_L(m)*Iz{1}+2*pi*v_L(m)*Iz{1}*CS_zz(m)+2*pi*HF_zz(m)*(Sz*Iz{1})+2*pi*HF_zy(m)*(Sz*Iy{1})+2*pi*HF_zx(m)*(Sz*Ix{1});
                % NQI
                if parameters.Bterm==true
                    Hfree_p=Hfree_p+NQI(m,1,1)*Ix{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(m,1,2)*Ix{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(m,1,3)*Ix{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(m,2,1)*Iy{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(m,2,2)*Iy{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(m,2,3)*Iy{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(m,3,1)*Iz{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(m,3,2)*Iz{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(m,3,3)*Iz{1}*Iz{1};
                else
                    Hfree_p=Hfree_p+pi*NQI_zz(m)*(3*Iz{1}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                end
            else
                for mm=1:Ni_ENDOR
                    % HF
                    Hfree_p=Hfree_p-2*pi*v_L(mm)*Iz{mm}+2*pi*HF_zz(mm)*(Sz*Iz{mm})+2*pi*HF_zy(mm)*(Sz*Iy{mm})+2*pi*HF_zx(mm)*(Sz*Ix{mm})+2*pi*v_L(mm)*CS_zz(mm)*Iz{mm};
                    % NQI
                    if parameters.Bterm==true
                        Hfree_p=Hfree_p+NQI(mm,1,1)*Ix{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,1,2)*Ix{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,1,3)*Ix{mm}*Iz{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,1)*Iy{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,2)*Iy{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,3)*Iy{mm}*Iz{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,1)*Iz{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,2)*Iz{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,3)*Iz{mm}*Iz{mm};
                    else
                        Hfree_p=Hfree_p+pi*NQI_zz(mm)*(3*Iz{mm}*Iz{mm}-I(mm)*(I(mm)+1)*eye(size(Hfree_p)));
                    end
                end
            end

            % Integration step for the Signal to account for oscillation
            if v_off_S==0
                t2=1/(off_1*Nint);
            else
                t2=1/(v_off_S*Nint);
            end


            % loop over rf frequencies (x-axis)
            for a=1:Npts_EN
            % for a=1   %%% for testing

                % Radiofrequency
                v_RF=(start_EN+step_EN*(a-1));

                Hcorr=zeros(size(Hfree_p));
                HRF=Hfree_p;

                if parameters.Bterm==false
                    if spinSys('N_SpinSys')>1
                        m=1;
                        Hcorr=Hcorr+2*pi*v_RF*Iz{m};

                        HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Iy{m};
                    else
                        for mm=1:Ni_ENDOR
                            Hcorr=Hcorr+2*pi*v_RF*Iz{mm};

                            HRF=HRF+2*pi*v_RF*Iz{mm}+oneN*Iy{mm};
                        end
                    end
                end

                Hfree=Hfree_p+Hcorr;
                Hfree=full(hilb2liouv(sparse(Hfree),'comm'));


                if parameters.Bterm==false

                    U1=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*full(hilb2liouv(sparse(HRF),'comm'))),t(1)));
                    U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t2));
                else
                    U1=kehl_rf_bterm_rlx(parameters,expt,v_RF,Hfree_p,Iy,t(1),Ni_ENDOR,N_spinSys,R,spinOps('spin_system'));
                    U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t2));
                end

                % Evolve the densitymatrix

                rho=kehl_mat_to_lbra(rho0)';
                rho=U1*rho;

                value_Sy=0;
                for b=1
                   rho=U2*rho;
                   rho_f=kehl_lket_to_mat(rho);
                   value_Sy=value_Sy+(real(trace(rho_f*Sz*Iz{1})));
                end
                if parameters.temp_eff==true
                    endor_amp_tmp(a)=endor_amp_tmp(a)+abs(value_Sy*S/(Nint*size(mI,2)))*const_R;
                else
                    endor_amp_tmp(a)=endor_amp_tmp(a)+abs(value_Sy*S/(Nint*size(mI,2)));

                end

            end
            end
        end
        endor_amp=endor_amp+endor_amp_tmp;
    end
end

function kehl_tensor_rlx_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isa(paramsENDOR,'containers.Map')
    error('paramsENDOR must be a containers.Map object.');
end
if ~isa(EPR,'containers.Map')
    error('EPR must be a containers.Map object.');
end
end


% --- Localised from kehl_tensor_calc.m ---
% performs the actual ENDOR calculation (no relaxation)
%
% input parameters:
% constants: the Map containing the constants
% spinOps: the Map containing the spin operators
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
% parameters: structure containing simulation parameters
% paramsENDOR: the Map containing the ENDOR parameters
% EPR: the Map with results from EPR calculation
%
% output parameters:
% endor_amp: endor amplitude values
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

% redefine, one spin system is not possible

function endor_amp=kehl_tensor_calc(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)

    % Check consistency
    kehl_tensor_calc_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR);
    spinSys('N_SpinSys')=spinSys("Ni_ENDOR");
    spinOps=kehl_spin_ops(spinSys("S"),spinSys("I"),spinSys("Ni_ENDOR"),spinSys("Ni_ENDOR"));

    % get values from Maps
    Sz=spinOps("Sz");
    Sx=spinOps("Sx");
    Sy=spinOps("Sy");
    Iz=spinOps("Iz");
    Ix=spinOps("Ix");
    Iy=spinOps("Iy");
    t=expt("t");
    Nint=8;

    N_spinSys=spinSys("Ni_ENDOR");
    Ni_ENDOR=spinSys("Ni_ENDOR");

    if N_spinSys>1
       nucs=Ni_ENDOR;
        for i=1:Ni_ENDOR-1
            Ix{i+1}=Ix{1};
            Iy{i+1}=Iy{1};
            Iz{i+1}=Iz{1};
        end
    else
       nucs=1;
    end



    geff_sel=EPR("geff_sel");
    B_sel=EPR("B_sel");
    HF_zz_sel=EPR("HF_zz_sel");
    HF_zy_sel=EPR("HF_zy_sel");
    HF_zx_sel=EPR("HF_zx_sel");

    NQI_zz_sel=EPR("NQI_zz_sel");
    NQI_sel=EPR("NQI_sel");
    CS_zz_sel=EPR("CS_zz_sel");
    D_zz_sel=EPR("D_zz_sel");

    S_sel=EPR("S_sel");


    offsets_sel=EPR("offsets");
    Npts_EN=paramsENDOR("Npts_EN");

    I=spinSys("I");
    v_L=paramsENDOR("v_L");




    endor_amp=zeros(1,Npts_EN);

    if length(B_sel)==0
        % No resonance orientations were found
        return
    end

    % loop to repeat the calculation for every orientation
    parfor j=1:length(B_sel)
        % set parameters for this orientation
        geff=geff_sel(j);
        B=B_sel(j);
        const_R=1/size(Sz,2)*constants('GE')*B/(2*pi*constants('K_B')*parameters.T);

        HF_zz=HF_zz_sel(j,:);
        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);
        NQI=zeros(Ni_ENDOR,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);

        CS_zz=CS_zz_sel(j,:)*1e-12;
        D_zz=D_zz_sel(j,:);


        S=S_sel(j);
        offsets=offsets_sel(j,:);

        endor_amp_tmp=zeros(1,Npts_EN);

        % loop over nuclei
        for i=1 : Ni_ENDOR
            I=spinSys('I');
            mI=-I(i):1:I(i);

            for jj=mI

                v_off_S=jj*HF_zz(i);
                off_1=-I(i)*HF_zz(i);

                [rho0]=2*Sz*Iz{i};

                start_EN=paramsENDOR("start_EN");
                step_EN=paramsENDOR("step_EN");

                oneE=expt("oneE");
                oneN=expt("oneN");

                Hfree_p=2*pi*v_off_S*Sz;
                nuc=i
                % HF
                Hfree_p=Hfree_p-2*pi*v_L(nuc)*Iz{1}+2*pi*v_L(nuc)*Iz{1}*CS_zz(nuc)+2*pi*HF_zz(nuc)*(Sz*Iz{1})+2*pi*HF_zy(nuc)*(Sz*Iy{1})+2*pi*HF_zx(nuc)*(Sz*Ix{1});
                % NQI
                if parameters.Bterm==true
                    m=1;
                    Hfree_p=Hfree_p+NQI(m,1,1)*Ix{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(m,1,2)*Ix{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(m,1,3)*Ix{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(m,2,1)*Iy{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(m,2,2)*Iy{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(m,2,3)*Iy{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(m,3,1)*Iz{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(m,3,2)*Iz{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(m,3,3)*Iz{1}*Iz{1};
                else
                    m=1;
                    Hfree_p=Hfree_p+pi*NQI_zz(m)*(3*Iz{m}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                end

                if spinSys("D_used")==true && i==1
                    for mm=2:(size(D_zz,2)+1)
                        dipC=2*pi*D_zz(mm-1);

                        HD=zeros(size(Hfree_p));

                        HD=HD+Ix{1}*Ix{mm};
                        HD=HD+Ix{1}*Iy{mm};
                        HD=HD+Ix{1}*Iz{mm};
                        HD=HD+Iy{1}*Ix{mm};
                        HD=HD+Iy{1}*Iy{mm};
                        HD=HD+Iy{1}*Iz{mm};
                        HD=HD+Iz{1}*Ix{mm};
                        HD=HD+Iz{1}*Iy{mm};
                        HD=HD+Iz{1}*Iz{mm};

                        HDip=dipC*(3*Iz{1}*Iz{mm}-HD)/2;
                        Hfree_p=Hfree_p+HDip;
                    end
                end

                % mw pulses
                Hnonsel_p=Hfree_p+oneE*Sx;

                % Integration step for the Signal to account for oscillation
                if v_off_S==0
                    t2=1/abs(off_1*Nint);
                else
                    t2=1/abs(v_off_S*Nint);
                end



                U2_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t2));


            % loop over rf frequencies (x-axis)
            for a=1:Npts_EN
            % for a=1   %%% for testing

                % Radiofrequency
                v_RF=(start_EN+step_EN*(a-1));

                Hcorr=zeros(size(Hfree_p));
                HRF=Hfree_p;

                if parameters.Bterm==false
                   HRF=HRF+2*pi*v_RF*Iz{1}+oneN*Iy{1};
                end

                Hcorr=Hcorr+2*pi*v_RF*Iz{1};


                if parameters.Bterm==false
                    U1=full(propagator(spinOps('spin_system'),sparse(HRF),t(1)));
                    U2=U2_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t2));
                else
                    Hfree=Hfree_p;
                    U1=kehl_rf_bterm(parameters,expt,v_RF,Hfree,Iy,t(1),Ni_ENDOR,N_spinSys,spinOps('spin_system'));
                    U2=U2_p;
                end


                % Evolve the densitymatrix
                rho=rho0;
                rho=U1*rho*U1';

                if parameters.Bterm==true
                    rho=diag(diag(rho));
                end


                value_Sy=0;
                value_Sy=value_Sy+(real(trace(rho*Sy)));
                for b=1:Nint*10
                   rho=U2*rho*U2';
                   value_Sy=value_Sy+(real(trace(rho*Sz*Iz{1})));
                end

                if parameters.temp_eff==true
                    endor_amp_tmp(a)=endor_amp_tmp(a)+abs(value_Sy*S/(Nint*size(mI,2)))*const_R;
                else
                    endor_amp_tmp(a)=endor_amp_tmp(a)+abs(value_Sy*S/(Nint*size(mI,2)));

                end

            end
            end
        end
        endor_amp=endor_amp+endor_amp_tmp;
    end
end

function kehl_tensor_calc_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isa(paramsENDOR,'containers.Map')
    error('paramsENDOR must be a containers.Map object.');
end
if ~isa(EPR,'containers.Map')
    error('EPR must be a containers.Map object.');
end
end


% --- Localised from kehl_cp_calc_rlx.m ---
% performs the actual CP ENDOR calculation (with relaxation)
%
% input parameters:
% constants: the Map containing the constants
% spinOps: the Map containing the spin operators
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
% parameters: structure containing simulation parameters
% paramsENDOR: the Map containing the ENDOR parameters
% EPR: the Map with results from EPR calculation
%
% output parameters:
% endor_amp: endor amplitude values
%
% February 2024 A. Kehl (akehl@gwdg.de)
%


function endor_amp=kehl_cp_calc_rlx(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)

    % Check consistency
    kehl_cp_calc_rlx_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR);
    spinOps_D=kehl_diag_ops(spinOps,spinSys('Ni_ENDOR'));

    % get parameters from Maps
    Sz=spinOps("Sz");
    Sx=spinOps("Sx");
    Sy=spinOps("Sy");
    Iz=spinOps("Iz");
    Ix=spinOps("Ix");
    Iy=spinOps("Iy");

    Sz_D=spinOps_D("Sz");
    Sx_D=spinOps_D("Sx");
    Sy_D=spinOps_D("Sy");
    Iz_D=spinOps_D("Iz");
    Ix_D=spinOps_D("Ix");
    Iy_D=spinOps_D("Iy");


    t=expt("t");
    Nint=8;

    N_spinSys=spinSys("N_SpinSys");

    Ni_ENDOR=spinSys("Ni_ENDOR");
    Npts_CP=expt("Npts_CP");

    if expt("Npts_CP")>1
        step_CP=expt("range_CP")/(expt("Npts_CP")-1);
    else
        step_CP=0;
    end

    if N_spinSys>1
       nucs=Ni_ENDOR;
        for i=1:Ni_ENDOR-1
            Ix{i+1}=Ix{1};
            Iy{i+1}=Iy{1};
            Iz{i+1}=Iz{1};
        end
    else
       nucs=1;
    end

    geff_sel=EPR("geff_sel");
    B_sel=EPR("B_sel");
    HF_zz_sel=EPR("HF_zz_sel");
    HF_zy_sel=EPR("HF_zy_sel");
    HF_zx_sel=EPR("HF_zx_sel");

    NQI_zz_sel=EPR("NQI_zz_sel");
    NQI_sel=EPR("NQI_sel");

    CS_zz_sel=EPR("CS_zz_sel");
    D_zz_sel=EPR("D_zz_sel");


    S_sel=EPR("S_sel");


    offsets_sel=EPR("offsets");
    Npts_EN=paramsENDOR("Npts_EN");

    I=spinSys("I");
    v_L=paramsENDOR("v_L");


    endor_amp=zeros(1,Npts_EN);
    % length(B_sel)

    if length(B_sel)==0
        % No resonance orientations were found
        return
    end

    % loop to repeat the calculation for every orientation
    for j=1:length(B_sel)
        % set parameters for this orientation
        geff=geff_sel(j);
        B=B_sel(j);

        HF_zz=HF_zz_sel(j,:);
        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);
        NQI=zeros(Ni_ENDOR,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);

        CS_zz=CS_zz_sel(j,:)*1e-12;
        D_zz=D_zz_sel(j,:);



        S=S_sel(j);
        offsets=offsets_sel(j,:);

            % % loop over nuclei if not all in one spinSys
        % for n = 1 : nucs
        % loop over all offsets (spin manifolds)
        for i=1 : size(offsets,2)

            v_off_S=offsets(i);
            off_1=offsets(1);


            rho0=kehl_rho0(constants,paramsENDOR,B,geff,spinOps,spinSys,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            %electron T1
            RT1e=kehl_relax_t1(rho0,Sx_D,spinSys("T1e"));
            RT1n=zeros(size(RT1e));

            % electron T2
            RT2e=kehl_relax_t2(Sx_D,spinSys("T2e"));
            RT2n=zeros(size(RT2e));

            if N_spinSys==1
                for mm=1:Ni_ENDOR

                   %nuclear T1
                   RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{mm},spinSys("T1n"));

                   % double and zero quantum T2
                   RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{mm},spinSys("T2dq"));

                   % nuclear T2
                   RT2n=RT2n+kehl_relax_t2(Ix_D{mm},spinSys("T2n"));
                end
            else

                %nuclear T1
                RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{1},spinSys("T1n"));

                % double and zero quantum T2
                RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{1},spinSys("T2dq"));

                % nuclear T2
                RT2n=RT2n+kehl_relax_t2(Ix_D{1},spinSys("T2n"));
            end

            % full relaxation superoperator
            R=RT1e+RT1n+RT2e+RT2n;

            % for tilted frame Rho relaxation times

            %electron T1
            RT1eRho=kehl_relax_t1(rho0,Sx_D,spinSys("T1eR"));
            RT1nRho=zeros(size(RT1eRho));

            % electron T2
            RT2eRho=kehl_relax_t2(Sx_D,spinSys("T2eR"));
            RT2nRho=zeros(size(RT2eRho));

            if N_spinSys==1
                for mm=1:Ni_ENDOR

                   %nuclear T1
                   RT1nRho=RT1eRho+kehl_relax_t1(rho0,Ix_D{mm},spinSys("T1nR"));

                   % double and zero quantum T2
                   RT2eRho=RT2eRho+kehl_relax_t2(Sx_D*Ix_D{mm},spinSys("T2dqR"));

                   % nuclear T2
                   RT2nRho=RT2nRho+kehl_relax_t2(Ix_D{mm},spinSys("T2nR"));
                end
            else

                %nuclear T1
                RT1nRho=RT1eRho+kehl_relax_t1(rho0,Ix_D{1},spinSys("T1nR"));

                % double and zero quantum T2
                RT2eRho=RT2eRho+kehl_relax_t2(Sx_D*Ix_D{1},spinSys("T2dqR"));

                % nuclear T2
                RT2nRho=RT2nRho+kehl_relax_t2(Ix_D{1},spinSys("T2nR"));
            end

            % full relaxation superoperator
            RRho=RT1eRho+RT1nRho+RT2eRho+RT2nRho;

            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            prep=expt("prep");
            sl=expt("SL");
            cp=expt("CP");
            oneE=expt("oneE");
            oneN=expt("oneN");

            Hfree_p=2*pi*v_off_S*Sz;
            if spinSys('N_SpinSys')>1
                if parameters.powder==1
                    mn=round(i/(2*I(1)+1));
                else
                    mn=i;
                end
                % HF
                Hfree_p=Hfree_p-2*pi*v_L(mn)*Iz{1}+2*pi*v_L(mn)*Iz{1}*CS_zz(mn)+2*pi*HF_zz(mn)*(Sz*Iz{1})+2*pi*HF_zy(mn)*(Sz*Iy{1})+2*pi*HF_zx(mn)*(Sz*Ix{1});
                % NQI
                if parameters.Bterm==true
                    Hfree_p=Hfree_p+NQI(mn,1,1)*Ix{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,1,2)*Ix{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,1,3)*Ix{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(mn,2,1)*Iy{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,2,2)*Iy{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,2,3)*Iy{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(mn,3,1)*Iz{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,3,2)*Iz{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,3,3)*Iz{1}*Iz{1};
                else
                    Hfree_p=Hfree_p+pi*NQI_zz(mn)*(3*Iz{1}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                end

                % get offset for sc CP calculation
                if parameters.powder==false
                    [v_CP,paramsENDOR]=kehl_cp_offset(parameters,paramsENDOR,expt,spinSys,v_off_S,HF_zz,NQI_zz,mn,i);
                end

            else
                for mm=1:Ni_ENDOR
                    % HF
                    Hfree_p=Hfree_p-2*pi*v_L(mm)*Iz{mm}+2*pi*HF_zz(mm)*(Sz*Iz{mm})+2*pi*HF_zy(mm)*(Sz*Iy{mm})+2*pi*HF_zx(mm)*(Sz*Ix{mm})+2*pi*v_L(mm)*CS_zz(mm)*Iz{mm};
                    % NQI
                    if parameters.Bterm==true
                        Hfree_p=Hfree_p+NQI(mm,1,1)*Ix{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,1,2)*Ix{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,1,3)*Ix{mm}*Iz{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,1)*Iy{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,2)*Iy{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,3)*Iy{mm}*Iz{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,1)*Iz{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,2)*Iz{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,3)*Iz{mm}*Iz{mm};
                    else
                        Hfree_p=Hfree_p+pi*NQI_zz(mm)*(3*Iz{mm}*Iz{mm}-I(mm)*(I(mm)+1)*eye(size(Hfree_p)));
                    end

                    % get offset for sc CP calculation
                    if parameters.powder==false
                        [v_CP,paramsENDOR]=kehl_cp_offset(parameters,paramsENDOR,expt,spinSys,v_off_S,HF_zz,NQI_zz,mm,i);
                    end
                end
                 if spinSys("D_used")==true
                    for mm=2:(size(D_zz,2)+1)
                        dipC=2*pi*D_zz(mm-1);

                        HD=zeros(size(Hfree_p));

                        HD=HD+Ix{1}*Ix{mm};
                        HD=HD+Ix{1}*Iy{mm};
                        HD=HD+Ix{1}*Iz{mm};
                        HD=HD+Iy{1}*Ix{mm};
                        HD=HD+Iy{1}*Iy{mm};
                        HD=HD+Iy{1}*Iz{mm};
                        HD=HD+Iz{1}*Ix{mm};
                        HD=HD+Iz{1}*Iy{mm};
                        HD=HD+Iz{1}*Iz{mm};

                        HDip=dipC*(3*Iz{1}*Iz{mm}-HD)/2;
                        Hfree_p=Hfree_p+HDip;
                    end
                end
            end

            % Integration step for the Signal to account for oscillation
            if v_off_S==0
                t11=1/(off_1*Nint);
            else
                t11=1/(v_off_S*Nint);
            end

            % loop over cp frequencies (y-axis)
            for c=1:Npts_CP

            if parameters.powder==true
                v_CP=expt("start_CP")+step_CP*(c-1);
            end

            % loop over rf frequencies (x-axis)
            parfor a=1:Npts_EN
            % for a=1   %%% for testing

                if abs(HF_zz(1))>1/spinSys("T2e")*0.1

                % Radiofrequency
                v_RF=(start_EN+step_EN*(a-1));

                Hcorr=zeros(size(Hfree_p));
                % Hamiltonian for RF pulse (no HF enhancement)
                HRF=Hfree_p;
                HRFb=zeros(size(Hfree_p));
                HSL=Hfree_p+sl*Sy;


                if parameters.Bterm==false
                    if spinSys('N_SpinSys')>1
                        m=1;
                        HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Ix{m};
                        HRFb=2*pi*(v_RF-v_CP)*Iz{m};
                        HSL=HSL+2*pi*v_CP*Iz{m}+cp*Ix{m};
                    else
                        for mm=1:Ni_ENDOR
                            HRF=HRF+2*pi*v_RF*Iz{mm}+oneN*Ix{mm};
                            HRFb=2*pi*(v_RF-v_CP)*Iz{mm};
                            HSL=HSL+2*pi*v_CP*Iz{mm}+cp*Ix{mm};
                        end
                    end
                end

                [V,D]=eig(HSL);
                if ~issorted(diag(D))
                  [d,inds]=sort(diag(D));
                  D=diag(d);
                  V=V(:, inds);
                end
                %eigenvalues are order from negative large to positive large!
                W=V^(-1);

                HSL_t=W*HSL*W^(-1);


                if spinSys('N_SpinSys')>1
                    m=1;
                    Hcorr=Hcorr+2*pi*v_CP*Iz{m};
                else

                    for mm=1:Ni_ENDOR
                        Hcorr=Hcorr+2*pi*v_CP*Iz{mm};
                    end
                end


                Hfree=Hfree_p+Hcorr;

                Hnonsel=full(hilb2liouv(sparse(Hfree+oneE*Sx),'comm'));
                Hprep=full(hilb2liouv(sparse(Hfree+prep*Sx),'comm'));

                Hfree=full(hilb2liouv(sparse(Hfree),'comm'));


                    if parameters.Bterm==false
                        U3=full(propagator(spinOps('spin_system'),1i*sparse(RRho-1i*full(hilb2liouv(sparse(HSL_t),'comm'))),t(3)));

                        U5=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*full(hilb2liouv(sparse(HRF),'comm'))),t(5)));

                        U1=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hprep),t(1)));
                        U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)));

                        U4=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(4)));

                        U6=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(6)));
                        U7=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(7)));
                        U8=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(8)));
                        U9=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(9)));
                        U10=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(8)+t(7)/2));
                        U11=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t11));
                    else
                        t_stepCP=1/(v_CP*parameters.N_stepRF);
                        H_SL=HSL_t;
                        U_SL=eye(size(HSL_t,1)^2,size(HSL_t,1)^2);
                        for ll=1:parameters.N_stepRF
                            if N_spinSys==1
                                for mm=1:Ni_ENDOR
                                    H_SL=H_SL+2*expt('CP')*Ix{mm}*cos(2*pi*v_CP*t_stepCP*(ll-1));
                                end
                            else
                                H_SL=H_SL+2*expt('CP')*Ix{1}*cos(2*pi*v_CP*t_stepCP*(ll-1));
                            end
                            G=RRho/t_stepCP-1i*full(hilb2liouv(sparse(H_SL),'comm'));
                            U_step=full(propagator(spinOps('spin_system'),1i*sparse(G),t_stepCP));
                            U_SL=U_step*U_SL;
                        end
                        U3=(U_SL^(t(3)*v_CP));
                        U5=kehl_rf_bterm_rlx(parameters,expt,v_RF,Hfree_p,Iy,t(5),Ni_ENDOR,N_spinSys,R,spinOps('spin_system'));

                        U1=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hprep),t(1)));
                        U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)));

                        U4=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(4)));

                        U6=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(6)));
                        U7=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(7)));
                        U8=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(8)));
                        U9=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(9)));
                        U10=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(8)+t(7)/2));
                        U11=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t11));

                    end

                % Evolve the densitymatrix
                rho=kehl_mat_to_lbra(rho0)';
                rho=U1*rho;
                rho=U2*rho;

                rho_t=W*kehl_lket_to_mat(rho)*W^(-1);
                rho_t=kehl_mat_to_lbra(rho_t)';
                rho_t=U3*rho_t;
                rho=V*kehl_lket_to_mat(rho_t)*V^(-1);
                rho=kehl_mat_to_lbra(rho)';

                rho=U4*rho;
                rho=U5*rho;
                rho=U6*rho;
                rho=U7*rho;
                rho=U8*rho;
                rho=U9*rho;
                rho=U10*rho;

                value_Sy=0;
                for b=1:Nint
                   rho=U11*rho;
                   rho_f=kehl_lket_to_mat(rho);
                   value_Sy=value_Sy+real(trace(rho_f*Sy));
                end
                    endor_amp(a)=endor_amp(a)+(value_Sy*S/(Nint*size(offsets,2)));
                end
            end
            end

        end
    end
end

function kehl_cp_calc_rlx_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isa(paramsENDOR,'containers.Map')
    error('paramsENDOR must be a containers.Map object.');
end
if ~isa(EPR,'containers.Map')
    error('EPR must be a containers.Map object.');
end
end


% --- Localised from kehl_cp_calc.m ---
% performs the actual CP ENDOR calculation (no relaxation)
%
% input parameters:
% constants: the Map containing the constants
% spinOps: the Map containing the spin operators
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
% parameters: structure containing simulation parameters
% paramsENDOR: the Map containing the ENDOR parameters
% EPR: the Map with results from EPR calculation
%
% output parameters:
% endor_amp: endor amplitude values
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

% get values from Maps

function endor_amp=kehl_cp_calc(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)

    % Check consistency
    kehl_cp_calc_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR);
    Sz=spinOps("Sz");
    Sx=spinOps("Sx");
    Sy=spinOps("Sy");
    Iz=spinOps("Iz");
    Ix=spinOps("Ix");
    Iy=spinOps("Iy");
    t=expt("t");
    Nint=8;

    N_spinSys=spinSys("N_SpinSys");

    Ni_ENDOR=spinSys("Ni_ENDOR");

    if N_spinSys>1
       nucs=Ni_ENDOR;
        for i=1:Ni_ENDOR-1
            Ix{i+1}=Ix{1};
            Iy{i+1}=Iy{1};
            Iz{i+1}=Iz{1};
        end
    else
       nucs=1;
    end

    geff_sel=EPR("geff_sel");
    B_sel=EPR("B_sel");
    HF_zz_sel=EPR("HF_zz_sel");
    HF_zy_sel=EPR("HF_zy_sel");
    HF_zx_sel=EPR("HF_zx_sel");

    NQI_zz_sel=EPR("NQI_zz_sel");
    NQI_sel=EPR("NQI_sel");

    CS_zz_sel=EPR("CS_zz_sel");
    D_zz_sel=EPR("D_zz_sel");

    S_sel=EPR("S_sel");


    offsets_sel=EPR("offsets");
    Npts_EN=paramsENDOR("Npts_EN");
    Npts_CP=expt("Npts_CP");

    if expt("Npts_CP")>1
        step_CP=expt("range_CP")/(expt("Npts_CP")-1);
    else
        step_CP=0;
    end

    I=spinSys("I");
    v_L=paramsENDOR("v_L");


    endor_amp=zeros(1,Npts_EN);
     length(B_sel)
    if length(B_sel)==0
        % No resonance orientations were found
        return
    end


    % loop to repeat the calculation for every orientation
    for j=1:length(B_sel)
        % set parameters for this orientation
        geff=geff_sel(j);
        B=B_sel(j);

        HF_zz=HF_zz_sel(j,:);
        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);

        NQI=zeros(Ni_ENDOR,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);

        CS_zz=CS_zz_sel(j,:)*1e-12;
        D_zz=D_zz_sel(j,:);


        S=S_sel(j);
        offsets=offsets_sel(j,:);

        % % loop over nuclei if not all in one spinSys
        % for n = 1:nucs

        % loop over all offsets (spin manifolds)
        for i=1 : size(offsets,2)

            v_off_S=offsets(i);
            off_1=offsets(1);

            rho0=kehl_rho0(constants,paramsENDOR,B,geff,spinOps,spinSys,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            prep=expt("prep");
            sl=expt("SL");
            cp=expt("CP");
            oneE=expt("oneE");
            oneN=expt("oneN");

            % define free Hamiltonian
            Hfree_p=2*pi*v_off_S*Sz;
            if spinSys('N_SpinSys')>1
                if parameters.powder==1
                    mn=round(i/(2*I(1)+1));
                else
                    mn=i;
                end
                % HF
                Hfree_p=Hfree_p-2*pi*v_L(mn)*Iz{1}+2*pi*v_L(mn)*Iz{1}*CS_zz(mn)+2*pi*HF_zz(mn)*(Sz*Iz{1})+2*pi*HF_zy(mn)*(Sz*Iy{1})+2*pi*HF_zx(mn)*(Sz*Ix{1});
                % NQI
                if parameters.Bterm==true
                    Hfree_p=Hfree_p+NQI(mn,1,1)*Ix{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,1,2)*Ix{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,1,3)*Ix{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(mn,2,1)*Iy{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,2,2)*Iy{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,2,3)*Iy{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(mn,3,1)*Iz{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,3,2)*Iz{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,3,3)*Iz{1}*Iz{1};
                else
                    Hfree_p=Hfree_p+pi*NQI_zz(mn)*(3*Iz{1}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                end

                % get offset for sc CP calculation
                if parameters.powder==false
                    [v_CP,paramsENDOR]=kehl_cp_offset(parameters,paramsENDOR,expt,spinSys,v_off_S,HF_zz,NQI_zz,mn,i);
                end
            else
                for mm=1:Ni_ENDOR
                    % HF
                    Hfree_p=Hfree_p-2*pi*v_L(mm)*Iz{mm}+2*pi*HF_zz(mm)*(Sz*Iz{mm})+2*pi*HF_zy(mm)*(Sz*Iy{mm})+2*pi*HF_zx(mm)*(Sz*Ix{mm})+2*pi*v_L(mm)*CS_zz(mm)*Iz{mm};
                    % NQI
                    if parameters.Bterm==true
                        Hfree_p=Hfree_p+NQI(mm,1,1)*Ix{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,1,2)*Ix{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,1,3)*Ix{mm}*Iz{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,1)*Iy{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,2)*Iy{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,3)*Iy{mm}*Iz{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,1)*Iz{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,2)*Iz{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,3)*Iz{mm}*Iz{mm};
                    else
                        Hfree_p=Hfree_p+pi*NQI_zz(mm)*(3*Iz{mm}*Iz{mm}-I(mm)*(I(mm)+1)*eye(size(Hfree_p)));
                    end

                    % get offset for sc CP calculation
                    if parameters.powder==false
                        [v_CP,paramsENDOR]=kehl_cp_offset(parameters,paramsENDOR,expt,spinSys,v_off_S,HF_zz,NQI_zz,mm,i);
                    end
                end

                if spinSys("D_used")==true
                    for mm=2:(size(D_zz,2)+1)
                        dipC=2*pi*D_zz(mm-1);

                        HD=zeros(size(Hfree_p));

                        HD=HD+Ix{1}*Ix{mm};
                        HD=HD+Ix{1}*Iy{mm};
                        HD=HD+Ix{1}*Iz{mm};
                        HD=HD+Iy{1}*Ix{mm};
                        HD=HD+Iy{1}*Iy{mm};
                        HD=HD+Iy{1}*Iz{mm};
                        HD=HD+Iz{1}*Ix{mm};
                        HD=HD+Iz{1}*Iy{mm};
                        HD=HD+Iz{1}*Iz{mm};

                        HDip=dipC*(3*Iz{1}*Iz{mm}-HD)/2;
                        Hfree_p=Hfree_p+HDip;
                    end
                end
            end

                % Hamiltonian for mw pulses
                Hprep_p=Hfree_p+prep*Sx;
                Hnonsel_p=Hfree_p+oneE*Sx;

                % Integration step for the Signal to account for oscillation
                if v_off_S==0
                    t11=1/(off_1*Nint);
                else
                    t11=1/(v_off_S*Nint);
                end


                % Calculate the propagators
                U1_p=full(propagator(spinOps('spin_system'),sparse(Hprep_p),t(1)));
                U2_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(2)));
                % SL/CP step
                U4_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(4)));
                % RF pulse
                U6_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(6)));
                U7_p=full(propagator(spinOps('spin_system'),sparse(Hnonsel_p),t(7)));
                U8_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(8)));
                U9_p=full(propagator(spinOps('spin_system'),sparse(Hnonsel_p),t(9)));
                U10_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(8)+t(7)/2));
                U11_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t11));

            % loop over cp frequencies (y-axis)
            for c=1:Npts_CP

            if parameters.powder==true
                v_CP=expt("start_CP")+step_CP*(c-1);
            end


            % loop over rf frequencies (x-axis)
            parfor a=1:Npts_EN
            % for a=1   %%% for testing

                % Radiofrequency
                v_RF=(start_EN+step_EN*(a-1));

                % Calculate the Free evolution hamiltonian
                % (not including the Bterm)

                Hcorr=zeros(size(Hfree_p));
                % Hamiltonian for RF pulse (no HF enhancement)
                HRF=Hfree_p;
                HRFb=zeros(size(Hfree_p));

                HSL=Hfree_p+sl*Sy;

                [V,D]=eig(HSL);
                if ~issorted(diag(D))
                  [d,inds]=sort(diag(D));
                  D=diag(d);
                  V=V(:, inds);
                end
                %eigenvalues are order from negative large to positive large!
                W=V^(-1);


                if parameters.Bterm==false
                    if spinSys('N_SpinSys')>1
                        m=1;
                        HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Ix{m};
                        HRFb=2*pi*(v_RF-v_CP)*Iz{m};
                        HSL=HSL+2*pi*v_CP*Iz{m}+cp*Ix{m};
                    else
                        for mm=1:Ni_ENDOR
                            HRF=HRF+2*pi*v_RF*Iz{mm}+oneN*Ix{mm};
                            HRFb=2*pi*(v_RF-v_CP)*Iz{mm};
                            HSL=HSL+2*pi*v_CP*Iz{mm}+cp*Ix{mm};
                        end
                    end
                end

                if spinSys('N_SpinSys')>1
                    m=1;
                    Hcorr=Hcorr+2*pi*v_CP*Iz{m};
                else

                    for mm=1:Ni_ENDOR
                        Hcorr=Hcorr+2*pi*v_CP*Iz{mm};
                    end
                end

                U5b=[];
                if parameters.Bterm==false
                    U3=full(propagator(spinOps('spin_system'),sparse(HSL),t(3)));
                    U5=full(propagator(spinOps('spin_system'),sparse(HRF),t(5)));
                    U5b=full(propagator(spinOps('spin_system'),sparse(HRFb),t(5)));

                    U1=U1_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(1)));
                    U2=U2_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(2)));

                    U4=U4_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(4)));

                    U6=U6_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(6)));
                    U7=U7_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(7)));
                    U8=U8_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(8)));
                    U9=U9_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(9)));
                    U10=U10_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(8)+t(7)/2));
                    U11=U11_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t11));
                else
                    Hfree=Hfree_p;
                    HSL_p=Hfree+sl*Sy;

                    t_stepCP=1/(v_RF*parameters.N_stepRF);
                    U_SL=eye(size(HSL_p));
                    H_SL=HSL_p;
                    for ll=1:parameters.N_stepRF
                        if N_spinSys==1
                            for mm=1:Ni_ENDOR
                                H_SL=H_SL+2*expt('CP')*Ix{mm}*cos(2*pi*v_RF*t_stepCP*(ll-1));
                            end
                        else
                            H_SL=H_SL+2*expt('CP')*Ix{1}*cos(2*pi*v_RF*t_stepCP*(ll-1));
                        end
                        U_step=full(propagator(spinOps('spin_system'),sparse(H_SL),t_stepCP));
                        U_SL=U_step*U_SL;
                    end
                    U3=(U_SL^(t(3)*v_RF));
                    U5=kehl_rf_bterm(parameters,expt,v_RF,Hfree,Iy,t(5),Ni_ENDOR,N_spinSys,spinOps('spin_system'));

                    U1=U1_p;
                    U2=U2_p;

                    U4=U4_p;

                    U6=U6_p;
                    U7=U7_p;
                    U8=U8_p;
                    U9=U9_p;
                    U10=U10_p;
                    U11=U11_p;
                end


                % Evolve the densitymatrix
                rho=rho0;
                rho=U1*rho*U1';
                rho=U2*rho*U2';
                if t(2)~=0
                    rho=diag(diag(rho));
                end

                rho=U3*rho*U3';

                rho_t=W*rho*W^(-1);
                rho_t=diag(diag(rho_t));
                rho=V*rho_t*V^(-1);

                rho=U4*rho*U4';
                rho=U5*rho*U5';

                if parameters.Bterm==true
                    rho=diag(diag(rho));
                else
                    rho=U5b*rho*U5b';
                end

                rho=U6*rho*U6';
                rho=U7*rho*U7';
                rho=U8*rho*U8';
                rho=U9*rho*U9';
                rho=U10*rho*U10';


                value_Sy=0;
                for b=1:Nint*10
                   rho=U11*rho*U11';
                   value_Sy=value_Sy+real(trace(rho*Sy));
                end

                if spinSys('N_SpinSys')>1
                    if parameters.powder==0
                        s=1;
                    else
                        s=size(offsets,2)/Ni_ENDOR;
                    end
                else
                    s=size(offsets,2);
                end

                endor_amp(a)=endor_amp(a)+(value_Sy*S/(Nint*s));

            end
            end
        end

    end
end

function kehl_cp_calc_grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isa(paramsENDOR,'containers.Map')
    error('paramsENDOR must be a containers.Map object.');
end
if ~isa(EPR,'containers.Map')
    error('EPR must be a containers.Map object.');
end
end

