%ENDOR_KEHL_MIMS Mims ENDOR pulse sequence for the Kehl ENDOR context.
%
%   ENDOR_AMP=ENDOR_KEHL_MIMS(SPIN_SYSTEM,PARAMETERS,H,R,K) is called by
%   endor_kehl_context.m. H, R, and K are present for Spinach experiment
%   signature compatibility; this sequence uses orientation-selected
%   effective Hamiltonian data prepared by the context.


function endor_amp=endor_kehl_mims(spin_system,parameters,H,R,K)


% Append sequence-specific parameters when requested by the context
if nargin>=3 && ischar(H) && strcmp(H,'parameters')
    endor_amp=kehl_mims_parameters(spin_system,parameters);
    return
end
% Check consistency
grumble(spin_system,parameters,H,R,K);
if parameters.Relax==true
    endor_amp=kehl_mims_rlx(spin_system,parameters);
else
    endor_amp=kehl_mims_calc(spin_system,parameters);
end
end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if (~isempty(H))&&(~isnumeric(H))
    error('H must be empty or numeric.');
end
if (~isempty(R))&&(~isnumeric(R))
    error('R must be empty or numeric.');
end
if (~isempty(K))&&(~isnumeric(K))
    error('K must be empty or numeric.');
end
end


function parameters=kehl_mims_parameters(spin_system,parameters)

% Get pulse sequence timing and RF-field policy
constants=parameters.constants;
t=parameters.pulse_times_s;
[rf_nutations,rf_auto]=kehl_rf_policy(parameters);

% Set Mims RF nutation fields
if rf_auto==false
    if size(rf_nutations,1)==2
        parameters.electron_nutation=rf_nutations(1)*2*pi*1e6;
        parameters.nuclear_nutation=rf_nutations(2)*2*pi*1e3;
        parameters.pulse_width=parameters.electron_nutation/...
                         (2*pi*constants('CONST1')*1e10);
    else
        error('parameters.rf_nutation_freqs has incompatible dimensions.');
    end
else
    parameters.electron_nutation=2*pi/(4*t(1));
    if isfield(parameters,'rf_flip_angle_deg')
        parameters.nuclear_nutation=(parameters.rf_flip_angle_deg/180)*...
                                    2*pi/(2*t(5));
    else
        parameters.nuclear_nutation=2*pi/(2*t(5));
    end
    parameters.pulse_width=parameters.electron_nutation/...
                     (2*pi*constants('CONST1')*1e10);
end

% Append standard ENDOR sweep-axis data
parameters=kehl_endor_axis(spin_system,parameters,'endor');

end

function endor_amp=kehl_mims_rlx(spin_system,parameters)

    % Check consistency
    kehl_mims_rlx_grumble(spin_system,parameters);

    % Unpack context data
    constants=parameters.constants;
        paramsENDOR=parameters.paramsENDOR;
    EPR=parameters.epr;
    operator_spin_system=parameters.operator_spin_system;
    n_endor=parameters.n_endor;
    n_spin_systems=parameters.n_spin_systems;
    I=parameters.endor_spin_numbers;

    % Obtain electron operators directly from Spinach
    Sx=full(operator(operator_spin_system,'Lx',1));
    Sy=full(operator(operator_spin_system,'Ly',1));
    Sz=full(operator(operator_spin_system,'Lz',1));

    % Obtain nuclear operators directly from Spinach
    Ix=cell(1,operator_spin_system.comp.nspins-1);
    Iy=cell(1,operator_spin_system.comp.nspins-1);
    Iz=cell(1,operator_spin_system.comp.nspins-1);
    for n=1:numel(Ix)
        spin_idx=n+1;
        Ix{n}=full(operator(operator_spin_system,'Lx',spin_idx));
        Iy{n}=full(operator(operator_spin_system,'Ly',spin_idx));
        Iz{n}=full(operator(operator_spin_system,'Lz',spin_idx));
    end
    [Sx_D,Sy_D,Sz_D,Ix_D,Iy_D,Iz_D]=kehl_diag_ops(Sx,Sy,Sz,Ix,Iy,Iz,n_endor);

    % get values from Maps


    t=parameters.pulse_times_s;
    Nint=100;


    if n_spin_systems>1
       nucs=n_endor;
        for i=1:n_endor-1
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
        NQI=zeros(n_endor,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);
        CS_zz=CS_zz_sel(j,:)*1e-12;
        D_zz=D_zz_sel(j,:);

        S=S_sel(j);
        offsets=offsets_sel(j,:);

        endor_amp_tmp=zeros(1,Npts_EN);


        if abs(HF_zz(1))>1/parameters.T2e*0.1

            % loop over all offsets (spin manifolds)
            for i=1 : size(offsets,2)

                v_off_S=offsets(i);
                off_1=offsets(1);

                rho0=kehl_rho0(constants,paramsENDOR,B,geff,operator_spin_system,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

                %electron T1
                RT1e=kehl_relax_t1(rho0,Sx_D,parameters.T1e);
                RT1n=zeros(size(RT1e));
                RT2e=kehl_relax_t2(Sx_D,parameters.T2e);
                RT2n=zeros(size(RT2e));

                if n_spin_systems==1
                    for mm=1:n_endor

                       %nuclear T1
                       RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{mm},parameters.T1n);
                       RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{mm},parameters.T2dq);
                       RT2n=RT2n+kehl_relax_t2(Ix_D{mm},parameters.T2n);
                    end
                else
                    RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{1},parameters.T1n);
                    RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{1},parameters.T2dq);
                    RT2n=RT2n+kehl_relax_t2(Ix_D{1},parameters.T2n);
                end
                R=RT1e+RT1n+RT2e+RT2n;


                start_EN=paramsENDOR("start_EN");
                step_EN=paramsENDOR("step_EN");

                oneE=parameters.electron_nutation;
                oneN=parameters.nuclear_nutation;

                if n_spin_systems>1
                    if parameters.powder==1
                        spin_map=round(i/(2*I(1)+1));
                    else
                        spin_map=i;
                    end
                    term_map=struct();
                    use_dipolar=false;
                else
                    spin_map=1:n_endor;
                    term_map=struct();
                    use_dipolar=true;
                end
                Hfree_p=kehl_free_ham(parameters,paramsENDOR,operator_spin_system,...
                                        v_off_S,spin_map,HF_zz,HF_zy,HF_zx,...
                                        NQI,NQI_zz,CS_zz,D_zz,use_dipolar,term_map);
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
                        if n_spin_systems>1
                            m=1;
                            Hcorr=Hcorr+2*pi*v_RF*Iz{m};

                            HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Iy{m};
                        else
                            for mm=1:n_endor
                                Hcorr=Hcorr+2*pi*v_RF*Iz{mm};

                                HRF=HRF+2*pi*v_RF*Iz{mm}+oneN*Iy{mm};
                            end
                        end
                    end

                    Hfree=Hfree_p+Hcorr;

                    Hnonsel=full(hilb2liouv(sparse(Hfree+oneE*Sx),'comm'));

                    Hfree=full(hilb2liouv(sparse(Hfree),'comm'));


                    if parameters.Bterm==false

                        U5=full(propagator(operator_spin_system,1i*sparse(R-1i*full(hilb2liouv(sparse(HRF),'comm'))),t(5)));

                        U1=full(propagator(operator_spin_system,1i*sparse(R-1i*Hnonsel),t(1)));
                        U2=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(2)));
                        U3=full(propagator(operator_spin_system,1i*sparse(R-1i*Hnonsel),t(3)));
                        U4=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(4)));

                        U6=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(6)));
                        U7=full(propagator(operator_spin_system,1i*sparse(R-1i*Hnonsel),t(7)));
                        U8=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(2)+t(3)/2));
                        U9=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t9));
                    else
                        U5=kehl_rf_bterm_rlx(parameters,v_RF,Hfree_p,Iy,t(5),n_endor,n_spin_systems,R,operator_spin_system);

                        U1=full(propagator(operator_spin_system,1i*sparse(R-1i*Hnonsel),t(1)));
                        U2=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(2)));
                        U3=full(propagator(operator_spin_system,1i*sparse(R-1i*Hnonsel),t(3)));
                        U4=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(4)));

                        U6=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(6)));
                        U7=full(propagator(operator_spin_system,1i*sparse(R-1i*Hnonsel),t(7)));
                        U8=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(2)+t(3)/2));
                        U9=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t9));

                    end


                    % Evolve the densitymatrix
                    rho=hilb2liouv(rho0,'statevec');

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
                       rho_f=reshape(rho,sqrt(size(rho,1)),sqrt(size(rho,1)));
                       value_Sy=value_Sy+(real(trace(rho_f*Sy)));
                    end
                    endor_amp_tmp(a)=endor_amp_tmp(a)+(value_Sy*S/(Nint*size(offsets,2)));
                end

            end
        end
        endor_amp=endor_amp+endor_amp_tmp;
    end
end

function kehl_mims_rlx_grumble(spin_system,parameters)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
end

function endor_amp=kehl_mims_calc(spin_system,parameters)

    % Check consistency
    kehl_mims_calc_grumble(spin_system,parameters);

    % Unpack context data
    constants=parameters.constants;
        paramsENDOR=parameters.paramsENDOR;
    EPR=parameters.epr;
    operator_spin_system=parameters.operator_spin_system;
    n_endor=parameters.n_endor;
    n_spin_systems=parameters.n_spin_systems;
    I=parameters.endor_spin_numbers;

    % Obtain electron operators directly from Spinach
    Sx=full(operator(operator_spin_system,'Lx',1));
    Sy=full(operator(operator_spin_system,'Ly',1));
    Sz=full(operator(operator_spin_system,'Lz',1));

    % Obtain nuclear operators directly from Spinach
    Ix=cell(1,operator_spin_system.comp.nspins-1);
    Iy=cell(1,operator_spin_system.comp.nspins-1);
    Iz=cell(1,operator_spin_system.comp.nspins-1);
    for n=1:numel(Ix)
        spin_idx=n+1;
        Ix{n}=full(operator(operator_spin_system,'Lx',spin_idx));
        Iy{n}=full(operator(operator_spin_system,'Ly',spin_idx));
        Iz{n}=full(operator(operator_spin_system,'Lz',spin_idx));
    end
    t=parameters.pulse_times_s;
    Nint=8;

    if n_spin_systems>1
       nucs=n_endor;
        for i=1:n_endor-1
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
        NQI=zeros(n_endor,3,3);

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

            [rho0]=kehl_rho0(constants,paramsENDOR,B,geff,operator_spin_system,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            oneE=parameters.electron_nutation;
            oneN=parameters.nuclear_nutation;

            if n_spin_systems>1
                if parameters.powder==1
                    spin_map=round(i/(2*I(1)+1));
                else
                    spin_map=i;
                end
                term_map=struct();
                use_dipolar=false;
            else
                spin_map=1:n_endor;
                term_map=struct();
                use_dipolar=true;
            end
            Hfree_p=kehl_free_ham(parameters,paramsENDOR,operator_spin_system,...
                                    v_off_S,spin_map,HF_zz,HF_zy,HF_zx,...
                                    NQI,NQI_zz,CS_zz,D_zz,use_dipolar,term_map);
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
            U1_p=full(propagator(operator_spin_system,sparse(Hnonsel_p),t(1)));
            U2_p=full(propagator(operator_spin_system,sparse(Hfree_p),t(2)));
            if t(3)==t(1)
                U3_p=U1_p;
            else
                U3_p=full(propagator(operator_spin_system,sparse(Hnonsel_p),t(3)));
            end

            U4_p=full(propagator(operator_spin_system,sparse(Hfree_p),t(4)));

            if t(4)==t(6)
                U6_p=U4_p;
            else

                U6_p=full(propagator(operator_spin_system,sparse(Hfree_p),t(6)));
            end

            if t(7)==t(1)
                U7_p=U1_p;
            else
                U7_p=full(propagator(operator_spin_system,sparse(Hnonsel_p),t(7)));
            end
            U8_p=full(propagator(operator_spin_system,sparse(Hfree_p),t(2)+t(3)/2));
            U9_p=full(propagator(operator_spin_system,sparse(Hfree_p),t9));


            % loop over rf frequencies (x-axis)
            for a=1:Npts_EN
            % for a=1   %%% for testing

                % Radiofrequency
                v_RF=(start_EN+step_EN*(a-1));

                Hcorr=zeros(size(Hfree_p));
                HRF=Hfree_p;

                if parameters.Bterm==false
                    if n_spin_systems>1
                        m=1;
                        HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Iy{m};
                    else
                        for mm=1:n_endor
                            HRF=HRF+2*pi*v_RF*Iz{mm}+oneN*Iy{mm};
                        end
                    end
                end

                if n_spin_systems>1
                    m=1;
                    Hcorr=Hcorr+2*pi*v_RF*Iz{m};
                else

                    for mm=1:n_endor
                        Hcorr=Hcorr+2*pi*v_RF*Iz{mm};

                    end
                end


                if parameters.Bterm==false
                    U5=full(propagator(operator_spin_system,sparse(HRF),t(5)));

                    U1=U1_p*full(propagator(operator_spin_system,sparse(Hcorr),t(1)));
                    U2=U2_p*full(propagator(operator_spin_system,sparse(Hcorr),t(2)));
                    if t(3)==t(1)
                        U3=U1;
                    else
                        U3=U3_p*full(propagator(operator_spin_system,sparse(Hcorr),t(3)));
                    end
                    U4=U4_p*full(propagator(operator_spin_system,sparse(Hcorr),t(4)));

                    U4=U4_p*full(propagator(operator_spin_system,sparse(Hcorr),t(4)));

                    if t(4)==t(6)
                        U6=U4;
                    else
                        U6=U6_p*full(propagator(operator_spin_system,sparse(Hcorr),t(6)));
                    end

                    if t(7)==t(1)
                        U7=U1;
                    else
                        U7=U7_p*full(propagator(operator_spin_system,sparse(Hcorr),t(7)));
                    end
                    U8=U8_p*full(propagator(operator_spin_system,sparse(Hcorr),t(2)+t(3)/2));
                    U9=U9_p*full(propagator(operator_spin_system,sparse(Hcorr),t9));
                else
                    Hfree=Hfree_p;

                    U5=kehl_rf_bterm(parameters,v_RF,Hfree,Iy,t(5),n_endor,n_spin_systems,operator_spin_system);

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

                if n_spin_systems>1
                    if parameters.powder==0
                        s=1;
                    else
                        s=size(offsets,2)/n_endor;
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

function kehl_mims_calc_grumble(spin_system,parameters)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
end
