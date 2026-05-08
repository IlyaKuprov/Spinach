%ENDOR_KEHL_CONTEXT Spinach-style ENDOR context for the Kehl pulse sequences.
%
%   [ENDOR,ENDOR_LB,X,V]=ENDOR_KEHL_CONTEXT(SPIN_SYSTEM,PULSE_SEQUENCE,...
%   PARAMETERS,ASSUMPTIONS) builds Kehl experimental parameters from the
%   physical fields in PARAMETERS, prepares the Spinach-backed spin
%   operators, EPR orientation selection, ENDOR sweep axes, and line
%   broadening, then calls PULSE_SEQUENCE with Spinach's standard
%   (spin_system,parameters,H,R,K) experiment signature.
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
pulse_sequence=sequence_handle(sequence);

% Store context assumptions
parameters.assumptions=assumptions;

% Build sequence-agnostic Kehl context data from Spinach parameters
parameters.constants=context_constants(spin_system);
parameters=context_fields(parameters);
parameters=context_spin_data(spin_system,parameters);
parameters.operator_spin_system=kehl_spin_system(parameters.operator_isotopes,...
                                                parameters.n_spin_systems);

if nargin>=4 && ~isempty(assumptions)
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

% Run the requested pulse sequence using the Spinach experiment signature
H=[];
R=[];
K=[];
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

function constants=context_constants(spin_system)
    constants=containers.Map;
    constants('H')=2*pi*spin_system.tols.hbar;
    constants('K_B')=spin_system.tols.kbol/constants('H');
    constants('MU_B')=spin_system.tols.muB;
    constants('GE')=spin_system.tols.freeg*spin_system.tols.muB/constants('H');
    constants('CONST1')=constants('GE')/1e10;
end

function parameters=context_spin_data(spin_system,parameters)
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

    parameters.electron_spin=(electron_mult-1)/2;
    parameters.electron_isotope=isotopes{electron_idx};
    parameters.g_matrix=electron_g_matrix(spin_system,electron_idx);
    parameters.g_iso=trace(parameters.g_matrix)/3;
    parameters.n_endor=n_endor;
    parameters.endor_isotopes=isotopes(endor_spins);
    parameters.n_spin_systems=n_spin_systems;
    parameters.endor_spin_numbers=spin_numbers(isotopes(endor_spins));
    parameters.operator_isotopes=[{parameters.electron_isotope} parameters.endor_isotopes(:).'];

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
    parameters.hfc_matrix=A;
    parameters.nqi_matrix=Q;
    parameters.nqi_active=Q_used;

    CS=zeros(3*n_endor,3);
    CS_used=false;
    for n=1:n_endor
        spin_idx=endor_spins(n);
        CS_block=nuclear_cs_matrix(spin_system,spin_idx);
        CS(3*n-2:3*n,:)=CS_block;
        CS_used=CS_used||any(CS_block(:));
    end
    if CS_used
        parameters.cs_matrix=CS;
    end
    parameters.cs_active=CS_used;

    pairs=infer_dipolar_pairs(spin_system,endor_spins);
    if ~isempty(pairs)
        D=zeros(3*size(pairs,1),3);
        for n=1:size(pairs,1)
            D(3*n-2:3*n,:)=coupling_matrix(spin_system,pairs(n,1),pairs(n,2));
        end
        parameters.dipolar_matrix=D;
        parameters.dipolar_active=true;
    else
        parameters.dipolar_matrix=zeros(0,3);
        parameters.dipolar_active=false;
    end

    parameters.epr_nuclei_active=~isempty(epr_spins);
    if ~isempty(epr_spins)
        n_epr=numel(epr_spins);
        parameters.epr_isotopes=isotopes(epr_spins);
        parameters.n_epr=n_epr;
        parameters.epr_spin_numbers=spin_numbers(isotopes(epr_spins));
        g_N_EPR=zeros(n_epr,1);
        A_EPR=zeros(3*n_epr,3);
        Q_EPR=zeros(3*n_epr,3);
        EPR_Q_used=false;
        for n=1:n_epr
            spin_idx=epr_spins(n);

            % Convert Spinach magnetogyric ratio to Hz/T
            [gamma,~]=spin(isotopes{spin_idx});
            g_N_EPR(n)=gamma/(2*pi);
            A_EPR(3*n-2:3*n,:)=coupling_matrix(spin_system,electron_idx,spin_idx);
            Q_block=coupling_matrix(spin_system,spin_idx,spin_idx);
            Q_EPR(3*n-2:3*n,:)=Q_block;
            EPR_Q_used=EPR_Q_used||any(Q_block(:));
        end
        parameters.epr_gamma_hz_t=g_N_EPR;
        parameters.epr_hfc_matrix=A_EPR;
        parameters.epr_nqi_matrix=Q_EPR;
        parameters.epr_nqi_active=EPR_Q_used;
    else
        parameters.epr_isotopes={};
        parameters.n_epr=0;
        parameters.epr_spin_numbers=zeros(0,1);
        parameters.epr_gamma_hz_t=zeros(0,1);
        parameters.epr_hfc_matrix=zeros(0,3);
        parameters.epr_nqi_matrix=zeros(0,3);
        parameters.epr_nqi_active=false;
    end
end

function parameters=context_fields(parameters)

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
% Localised Kehl context helper functions


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


function paramsEPR=kehl_prep_epr(spin_system,parameters)

    % Check consistency
    kehl_prep_epr_grumble(spin_system,parameters);

    % Unpack context data
    constants=parameters.constants;
        g_iso=parameters.g_iso;
    obsField=parameters.field_t;

    gBepr=obsField*g_iso;
    geff=g_iso;
    B=gBepr/geff;
    veff=geff*obsField*constants('MU_B')/constants('H');

    % EPR simulation range
    fieldCenter=gBepr/g_iso;
    fieldmin=fieldCenter-0.260;
    fieldmax=fieldCenter+0.260;

    % EPR x-axis definition

    % no of points in Field dimension
    Npts=(round((fieldmax-fieldmin)/parameters.field_step_t)+1);
    field=linspace(fieldmin,fieldmax,Npts);

    paramsEPR=containers.Map;
    paramsEPR("Npts")=Npts;
    paramsEPR("fieldAxis")=field;

end

function kehl_prep_epr_grumble(spin_system,parameters)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
end

function EPR=kehl_ori_field(spin_system,parameters)

    % Check consistency
    kehl_ori_field_grumble(spin_system,parameters);

    % Unpack context data
    constants=parameters.constants;
    paramsEPR=parameters.paramsEPR;
    paramsENDOR=parameters.paramsENDOR;
        Ntheta=parameters.Nang;
    Nphimax=parameters.Nang;
    g=parameters.g_matrix;
    A=parameters.hfc_matrix;
    Q=parameters.nqi_matrix;
    if parameters.epr_nuclei_active==true
        Ni_EPR=parameters.n_epr;
        A_EPR=parameters.epr_hfc_matrix;
        Q_EPR=parameters.epr_nqi_matrix;
    else
        Ni_EPR=0;
    end
    if parameters.cs_active==true
        CS=parameters.cs_matrix;
    end
    if parameters.dipolar_active==true
        D=parameters.dipolar_matrix;
    end
    n_endor=parameters.n_endor;
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
        if isfield(parameters,'excite_width')
            W1=parameters.excite_width;
        else
            W1=parameters.pulse_width*constants("CONST1")*1e10;
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
            B=(parameters.field_t*parameters.g_iso)/geff(1);
            Beff=parameters.mw_freq_hz*constants('H')/(constants('MU_B')*geff);

            % HF ENDOR
            HF_zz=zeros(1,n_endor);
            HF_zx=zeros(1,n_endor);
            HF_zy=zeros(1,n_endor);

            for m=1:n_endor

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
            NQI_zz=zeros(1,n_endor);
            NQI=zeros(n_endor,3,3);

            for m=1:n_endor
                if parameters.nqi_active==true

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

                    Q_ENDOR=parameters.nqi_matrix;
                    qq2=Q_ENDOR((m-1)*3+1:(m-1)*3+3,:);
                    Y2=R1*qq2*R1';
                    NQI(m,:,:)=Y2;

                else
                    NQI_zz(m)=0;
                    NQI (m,:,:)=0;
                end
            end

            % CS ENDOR
            CS_zz=zeros(1,n_endor);
            if parameters.cs_active==true
                for m=1:n_endor
                    CS_zz(m)=(sin(theta))^2*(cos(phi))^2*CS(3*m-2,1)...
+(sin(theta))^2*(sin(phi))^2*CS(3*m-1,2)...
+(cos(theta))^2*CS(3*m,3)...
+2*(sin(theta))^2*sin(phi)*cos(phi)*CS(3*m-2,2)...
+2*sin(theta)*cos(theta)*cos(phi)*CS(3*m-2,3)...
+2*sin(theta)*cos(theta)*sin(phi)*CS(3*m-1,3);
                end
            end

            % dipolar SSC ENDOR
            D_zz=zeros(1,n_endor);

             if parameters.dipolar_active==true
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
                I_EPR=parameters.epr_spin_numbers;
                mI=kehl_ori_field_BuildSpace(I_EPR);
                for i=1:Ni_EPR
                        % EPR resonance with hyperfine coupling only
                        E=Beff+mI*(HF_zz_EPR(i)/constants("CONST1")*10^(-10))';
                end
                bin=(round((E-fieldAxis(1))/parameters.field_step_t)+1);
            else
                E=Beff;
                bin=(round((E-fieldAxis(1))/parameters.field_step_t)+1);
            end

            for p=1:length(bin)

                % Add resonances to EPR spectrum
                tmp_epr(bin(p))=tmp_epr(bin(p))+1;
            end

            offsets=kehl_offsets(constants,parameters,parameters.operator_spin_system,paramsENDOR,B,geff,HF_zz,NQI_zz);

            DeltaB=0;
            scalefactor=0;

            for l=1:length(E)
                nor=nor+1;

                % Magnetic field offset in T
                DeltaB=fieldAxis(bin(l))-parameters.field_t;

                if isfield(parameters,'pulse_file')
                    scalefactor=kehl_ori_field_pulsescale(parameters,DeltaB,constants("CONST1"));
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
        geff=parameters.g_iso;

        %effective B field for given theta and phi
        Beff=parameters.mw_freq_hz*constants('H')/(constants('MU_B')*geff);

        HF_zz=zeros(1,n_endor);
        HF_zx=zeros(1,n_endor);
        HF_zy=zeros(1,n_endor);
        NQI_zz=zeros(1,n_endor);
        NQI=zeros(n_endor,3,3);


        CS_zz=zeros(1,n_endor);

        for m=1:n_endor

            HF_zz(m)=A(3*m,3);
            if parameters.Bterm==1
                HF_zy(m)=A(3*m,2);
                HF_zx(m)=A(3*m,1);
            end

            if parameters.nqi_active==true
                NQI_zz(m)=Q(3*m,3);
                NQI(m,:,:)=Q;
            end

            % CS ENDOR
            if parameters.cs_active==true
                CS_zz(m)=CS(3*m,3);
            end

        end

        % dipolar SSC ENDOR
        D_zz=zeros(1,n_endor);

        if parameters.dipolar_active==true
             D_zz=zeros(1,size(D,1)/3);
             for m=1:size(D,1)/3
                 D_zz(m)=D(3*m,3);
             end
         end


        offsets_tmp=kehl_offsets(constants,parameters,parameters.operator_spin_system,paramsENDOR,Beff,geff,HF_zz,NQI_zz);

        sel_I=parameters.sel_I;
        if parameters.n_spin_systems>1
            s=size(offsets_tmp,2)/n_endor;
            for n=1:n_endor
                offsets(n)=offsets_tmp(sel_I+(n-1)*s);
            end
        else
            size(offsets_tmp,2)
            s=size(offsets_tmp,2)/n_endor;
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

function scalefactor=kehl_ori_field_pulsescale(parameters,DeltaB,CONST1)
    % calculates the scalefactor in dependence of the mw pulse's excitation
    % profile
    %
    % input parameters:
    % parameters: Kehl ENDOR context parameters
    % DeltaB: offset of the actual field from the effective resonance field
    % CONST1: constant to convert between field and frequency
    %
    % output parameters:
    % scalefactor: scalefactor for the specific orientation
    %
    % February 2024 A. Kehl (akehl@gwdg.de)
    %

    data=parameters.pulse_file;
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
    if isfield(parameters,'multipulses')&&parameters.multipulses==true
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

function kehl_ori_field_grumble(spin_system,parameters)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
end

function EPR=kehl_ori_freq(spin_system,parameters)

    % Check consistency
    kehl_ori_freq_grumble(spin_system,parameters);

    % Unpack context data
    constants=parameters.constants;
    paramsEPR=parameters.paramsEPR;
    paramsENDOR=parameters.paramsENDOR;
        Ntheta=parameters.Nang;
    Nphimax=parameters.Nang;
    g=parameters.g_matrix;
    A=parameters.hfc_matrix;
    Q=parameters.nqi_matrix;
    if parameters.epr_nuclei_active==true
        Ni_EPR=parameters.n_epr;
        A_EPR=parameters.epr_hfc_matrix;
        Q_EPR=parameters.epr_nqi_matrix;
        g_N_EPR=parameters.epr_gamma_hz_t;
        I_EPR=parameters.epr_spin_numbers;
        epr_spin_system=kehl_spin_system([{parameters.electron_isotope} parameters.epr_isotopes(:).'],1);
        SzEPR=full(operator(epr_spin_system,'Lz',1));
        SxEPR=full(operator(epr_spin_system,'Lx',1));
        IzEPR=cell(1,Ni_EPR);
        IxEPR=cell(1,Ni_EPR);
        IyEPR=cell(1,Ni_EPR);
        for n=1:Ni_EPR
            spin_idx=n+1;
            IzEPR{n}=full(operator(epr_spin_system,'Lz',spin_idx));
            IxEPR{n}=full(operator(epr_spin_system,'Lx',spin_idx));
            IyEPR{n}=full(operator(epr_spin_system,'Ly',spin_idx));
        end
    else
        Ni_EPR=0;
        epr_spin_system=kehl_spin_system({parameters.electron_isotope},1);
        SzEPR=full(operator(epr_spin_system,'Lz',1));
        SxEPR=full(operator(epr_spin_system,'Lx',1));
    end
    if parameters.cs_active==true
        CS=parameters.cs_matrix;
    end
    if parameters.dipolar_active==true
        D=parameters.dipolar_matrix;
    end
    n_endor=parameters.n_endor;
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
        if isfield(parameters,'excite_width')
            W1=parameters.excite_width;
        else
            W1=parameters.pulse_width*constants("CONST1")*1e10;
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
                Beff=(parameters.field_t*parameters.g_iso)/geff(1);
                B=parameters.field_t;

                % Resonance frequency for geff at ObsField in GHz
                veff=geff*parameters.field_t*constants('MU_B')/constants('H');

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
                HF_zz=zeros(1,n_endor);
                HF_zx=zeros(1,n_endor);
                HF_zy=zeros(1,n_endor);
                NQI_zz=zeros(1,n_endor);
                NQI=zeros(n_endor,3,3);
                CS_zz=zeros(1,n_endor);
                D_zz=zeros(1,n_endor);

                for m=1:n_endor
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
                     if parameters.nqi_active==true

                        NQI_zz(m)=(sin(theta))^2*(cos(phi))^2*Q(3*m-2,1)...
+(sin(theta))^2*(sin(phi))^2*Q(3*m-1,2)...
+(cos(theta))^2*Q(3*m,3)...
+2*(sin(theta))^2*sin(phi)*cos(phi)*Q(3*m-2,2)...
+2*sin(theta)*cos(theta)*cos(phi)*Q(3*m-2,3)...
+2*sin(theta)*cos(theta)*sin(phi)*Q(3*m-1,3);
                        % Q11,Q22,Q33,Q12,Q13,Q23

                        % Q_ENDOR = parameters.nqi_matrix;
                        qq2=Q((m-1)*3+1:(m-1)*3+3,:);
                        Y2=R1*qq2*R1';
                        NQI(m,:,:)=Y2;
                     end

                     % CS ENDOR
                     if parameters.cs_active==true
                        CS_zz(m)=(sin(theta))^2*(cos(phi))^2*CS(3*m-2,1)...
+(sin(theta))^2*(sin(phi))^2*CS(3*m-1,2)...
+(cos(theta))^2*CS(3*m,3)...
+2*(sin(theta))^2*sin(phi)*cos(phi)*CS(3*m-2,2)...
+2*sin(theta)*cos(theta)*cos(phi)*CS(3*m-2,3)...
+2*sin(theta)*cos(theta)*sin(phi)*CS(3*m-1,3);
                     end

                end
                % nulcear dipolar-dipolar coupling ENDOR
                 if parameters.dipolar_active==true
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
                        if parameters.epr_nqi_active==true
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

               if parameters.epr_nuclei_active>0
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
                    bin_freq=round((freq_EPR-parameters.epr_freq_min_hz)/parameters.epr_freq_step_hz)+1;
                    if trans_prob_EPR(p)>0
                        tmp_epr(bin_freq(p))=tmp_epr(bin_freq(p))+trans_prob_EPR(p);
                        DeltaOm=freq_EPR(p)-parameters.mw_freq_hz;

                        if isfield(parameters,'pulse_file')
                            SF=kehl_ori_freq_pulsescale(parameters,DeltaOm)/2;
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

                    offsets=kehl_offsets(constants,parameters,parameters.operator_spin_system,paramsENDOR,B,geff,HF_zz,NQI_zz);


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
        geff=parameters.g_iso;

        %effective B field for given theta and phi
        B=(parameters.field_t*parameters.g_iso)/geff(1);

        HF_zz=zeros(1,n_endor);
        HF_zx=zeros(1,n_endor);
        HF_zy=zeros(1,n_endor);
        NQI_zz=zeros(1,n_endor);
        NQI=zeros(n_endor,3,3);


        CS_zz=zeros(1,n_endor);

        for m=1:n_endor

            HF_zz(m)=A(3*m,3);
            if parameters.Bterm==1
                HF_zy(m)=A(3*m,2);
                HF_zx(m)=A(3*m,1);
            end

            if parameters.nqi_active==true
                NQI_zz(m)=Q(3*m,3);
                NQI(m,:,:)=Q;
            end

            % CS ENDOR
            if parameters.cs_active==true
                CS_zz(m)=CS(3*m,3);
            end
        end

        % nuclear dipolar-dipolar coupling ENDOR
        D_zz=zeros(1,n_endor);
        if parameters.dipolar_active==true
             D_zz=zeros(1,size(D,1)/3);
             for m=1:size(D,1)/3
                 D_zz(m)=D(3*m,3);
             end
         end

        offsets_tmp=kehl_offsets(constants,parameters,parameters.operator_spin_system,paramsENDOR,B,geff,HF_zz,NQI_zz);

        sel_I=parameters.sel_I;
        if parameters.n_spin_systems>1
            s=size(offsets_tmp,2)/n_endor;
            for n=1:n_endor
                offsets(n)=offsets_tmp(sel_I+(n-1)*s);
            end
        else
            s=size(offsets_tmp,2)/n_endor;
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

function scalefactor=kehl_ori_freq_pulsescale(parameters,DeltaOm)
    % calculates the scalefactor in dependence of the mw pulse's excitation
    % profile
    %
    % input parameters:
    % parameters: Kehl ENDOR context parameters
    % DeltaOm: offset of the actual freq from the resonance freq
    %
    % output parameters:
    % scalefactor: scalefactor for the specific orientation
    %
    % February 2024 A. Kehl (akehl@gwdg.de)
    %
    data=parameters.pulse_file;
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
    if isfield(parameters,'multipulses')&&parameters.multipulses==true
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

function kehl_ori_freq_grumble(spin_system,parameters)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
end

function data_conv=kehl_line_broaden(data,parameters)

    % Check consistency
    kehl_line_broaden_grumble(data,parameters);

    % Get ENDOR sweep width
    sw=parameters.endor_range_hz;
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

function kehl_line_broaden_grumble(data,parameters)
if ~isnumeric(data)
    error('data must be numeric.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
end
