%ENDOR_KEHL_TENSOR tensor ENDOR pulse sequence for the Kehl ENDOR context.
%
%   ENDOR_AMP=ENDOR_KEHL_TENSOR(SPIN_SYSTEM,PARAMETERS,H,R,K) is called by
%   endor_kehl_context.m. H, R, and K are present for Spinach experiment
%   signature compatibility; this sequence uses orientation-selected
%   effective Hamiltonian data prepared by the context.


function endor_amp=endor_kehl_tensor(spin_system,parameters,H,R,K)


% Append sequence-specific parameters when requested by the context
if nargin>=3 && ischar(H) && strcmp(H,'parameters')
    endor_amp=kehl_tensor_parameters(spin_system,parameters);
    return
end
% Check consistency
grumble(spin_system,parameters,H,R,K);
if parameters.Relax==true
    endor_amp=kehl_tensor_rlx(spin_system,parameters);
else
    endor_amp=kehl_tensor_calc(spin_system,parameters);
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


function parameters=kehl_tensor_parameters(spin_system,parameters)

% Get pulse sequence timing and RF-field policy
constants=parameters.constants;
t=parameters.pulse_times_s;
[rf_nutations,rf_auto]=kehl_rf_policy(parameters);

% Set tensor RF nutation fields
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
    parameters.electron_nutation=2*pi/(2*t(2));
    parameters.nuclear_nutation=2*pi/(2*t(1));
    parameters.pulse_width=parameters.electron_nutation/...
                     (2*pi*constants('CONST1')*1e10);
end

% Append standard ENDOR sweep-axis data
parameters=kehl_endor_axis(spin_system,parameters,'endor');

end

function endor_amp=kehl_tensor_rlx(spin_system,parameters)

    % Check consistency
    kehl_tensor_rlx_grumble(spin_system,parameters);

    % Unpack context data
    constants=parameters.constants;
        paramsENDOR=parameters.paramsENDOR;
    EPR=parameters.epr;
    n_endor=parameters.n_endor;
    n_spin_systems=n_endor;
    I=parameters.endor_spin_numbers;
    operator_spin_system=kehl_spin_system(parameters.operator_isotopes,n_spin_systems);

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
    Nint=8;



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
        NQI=zeros(n_endor,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);
        CS_zz=CS_zz_sel(j,:)*1e-12;

        S=S_sel(j);
        offsets=offsets_sel(j,:);

        endor_amp_tmp=zeros(1,Npts_EN);

        % loop over nuclei
        for i=1 : n_endor
                    mI=-I(i):1:I(i);

            for jj=mI

                v_off_S=jj*HF_zz(i);
                off_1=-I(i)*HF_zz(i);

                [rho0]=2*Sz*Iz{1};

            RT2e=kehl_relax_t2(Sx_D,parameters.T2e);
            RT2n=zeros(size(RT2e));

            if n_spin_systems==1
                for mm=1:n_endor
                   RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{mm},parameters.T2dq);
                   RT2n=RT2n+kehl_relax_t2(Ix_D{mm},parameters.T2n);
                end
            else
                RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{1},parameters.T2dq);
                RT2n=RT2n+kehl_relax_t2(Ix_D{1},parameters.T2n);
            end
            R=RT2e+RT2n;
            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");


            oneN=parameters.nuclear_nutation;

            spin_map=1;
            Hfree_p=kehl_free_ham(parameters,paramsENDOR,operator_spin_system,...
                                    v_off_S,spin_map,HF_zz,HF_zy,HF_zx,...
                                    NQI,NQI_zz,CS_zz,[],false);
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
                Hfree=full(hilb2liouv(sparse(Hfree),'comm'));


                if parameters.Bterm==false

                    U1=full(propagator(operator_spin_system,1i*sparse(R-1i*full(hilb2liouv(sparse(HRF),'comm'))),t(1)));
                    U2=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t2));
                else
                    U1=kehl_rf_bterm_rlx(parameters,v_RF,Hfree_p,Iy,t(1),n_endor,n_spin_systems,R,operator_spin_system);
                    U2=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t2));
                end

                % Evolve the densitymatrix

                rho=hilb2liouv(rho0,'statevec');
                rho=U1*rho;

                value_Sy=0;
                for b=1
                   rho=U2*rho;
                   rho_f=reshape(rho,sqrt(size(rho,1)),sqrt(size(rho,1)));
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

function kehl_tensor_rlx_grumble(spin_system,parameters)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
end

function endor_amp=kehl_tensor_calc(spin_system,parameters)

    % Check consistency
    kehl_tensor_calc_grumble(spin_system,parameters);

    % Unpack context data
    constants=parameters.constants;
        paramsENDOR=parameters.paramsENDOR;
    EPR=parameters.epr;
    n_endor=parameters.n_endor;
    n_spin_systems=n_endor;
    I=parameters.endor_spin_numbers;
    operator_spin_system=kehl_spin_system(parameters.operator_isotopes,n_spin_systems);

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

    % get values from Maps
    t=parameters.pulse_times_s;
    Nint=8;

    n_spin_systems=n_endor;
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
        const_R=1/size(Sz,2)*constants('GE')*B/(2*pi*constants('K_B')*parameters.T);

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

        % loop over nuclei
        for i=1 : n_endor
                    mI=-I(i):1:I(i);

            for jj=mI

                v_off_S=jj*HF_zz(i);
                off_1=-I(i)*HF_zz(i);

                [rho0]=2*Sz*Iz{i};

                start_EN=paramsENDOR("start_EN");
                step_EN=paramsENDOR("step_EN");

                oneE=parameters.electron_nutation;
                oneN=parameters.nuclear_nutation;

                nuc=i;
                spin_map=nuc;
                term_map=struct('nqi',1);
                use_dipolar=(i==1);
                Hfree_p=kehl_free_ham(parameters,paramsENDOR,operator_spin_system,...
                                        v_off_S,spin_map,HF_zz,HF_zy,HF_zx,...
                                        NQI,NQI_zz,CS_zz,D_zz,use_dipolar,term_map);
                % mw pulses
                Hnonsel_p=Hfree_p+oneE*Sx;

                % Integration step for the Signal to account for oscillation
                if v_off_S==0
                    t2=1/abs(off_1*Nint);
                else
                    t2=1/abs(v_off_S*Nint);
                end



                U2_p=full(propagator(operator_spin_system,sparse(Hfree_p),t2));


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
                    U1=full(propagator(operator_spin_system,sparse(HRF),t(1)));
                    U2=U2_p*full(propagator(operator_spin_system,sparse(Hcorr),t2));
                else
                    Hfree=Hfree_p;
                    U1=kehl_rf_bterm(parameters,v_RF,Hfree,Iy,t(1),n_endor,n_spin_systems,operator_spin_system);
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

function kehl_tensor_calc_grumble(spin_system,parameters)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
end
