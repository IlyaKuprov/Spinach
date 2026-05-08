% Mims ENDOR pulse sequence for the Kehl ENDOR context. Syntax:
%
%      endor_amp=endor_kehl_mims(spin_system,parameters,H,R,K)
%
% Parameters:
%
%   spin_system      - Spinach spin system structure.
%   parameters       - Kehl ENDOR context parameter structure.
%   H,R,K            - Spinach experiment signature matrices.
%
% Outputs:
%
%   endor_amp        - simulated Mims ENDOR amplitude array.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=endor_kehl_mims.m>

function endor_amp=endor_kehl_mims(spin_system,parameters,H,R,K)

    % Append sequence-specific parameters when requested by the context
    if nargin>=3&&ischar(H)&&strcmp(H,'parameters')
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
    grumble(spin_system,parameters,[],[],[]);

    % Unpack context data
    constants=parameters.constants;
    paramsENDOR=parameters.paramsENDOR;
    EPR=parameters.epr;
    operator_spin_system=parameters.operator_spin_system;
    n_endor=parameters.n_endor;
    n_spin_systems=parameters.n_spin_systems;
    I=parameters.endor_spin_numbers;

    % Get cached operators and states
    ops=kehl_operator_basis(operator_spin_system,n_endor,n_spin_systems);
    Sx=ops.Sx;
    Sy=ops.Sy;
    Ix=ops.Ix;
    Iy=ops.Iy;
    Iz=ops.Iz;
    Sx_D=ops.Sx_D;
    Ix_D=ops.Ix_D;

    % Unpack context maps

    t=parameters.pulse_times_s;
    Nint=100;

    if n_spin_systems>1
        for spin_idx=1:n_endor-1
            Ix{spin_idx+1}=Ix{1};
            Iy{spin_idx+1}=Iy{1};
            Iz{spin_idx+1}=Iz{1};
        end
    end

    geff_sel=EPR("geff_sel");
    B_sel=EPR("B_sel");
    euler_sel=EPR("euler_sel");
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

    endor_amp=zeros(1,Npts_EN);

    if isempty(B_sel)
        % No resonance orientations were found
        return
    end

    % Loop over selected orientations
    parfor j=1:length(B_sel)
        % Select orientation-specific parameters
        geff=geff_sel(j);
        B=B_sel(j);
        euler_angles=euler_sel(j,:);

        HF_zz=HF_zz_sel(j,:);

        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);
        NQI=zeros(n_endor,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);

        S=S_sel(j);
        offsets=offsets_sel(j,:);

        endor_amp_tmp=zeros(1,Npts_EN);

        if abs(HF_zz(1))>1/parameters.T2e*0.1

            % Loop over spin-manifold offsets
            for offset_idx=1:size(offsets,2)

                v_off_S=offsets(offset_idx);
                off_1=offsets(1);

                rho0=kehl_rho0(constants,paramsENDOR,B,geff,operator_spin_system,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

                % Electron T1 relaxation
                RT1e=kehl_relax_t1(rho0,Sx_D,parameters.T1e);
                RT1n=zeros(size(RT1e));
                RT2e=kehl_relax_t2(Sx_D,parameters.T2e);
                RT2n=zeros(size(RT2e));

                if n_spin_systems==1
                    for mm=1:n_endor

                        % Nuclear T1 relaxation
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
                        spin_map=round(offset_idx/(2*I(1)+1));
                    else
                        spin_map=offset_idx;
                    end
                    term_map=struct();
                    use_dipolar=false;
                else
                    spin_map=1:n_endor;
                    term_map=struct();
                    use_dipolar=true;
                end
                Hfree_p=kehl_free_ham(parameters,paramsENDOR,operator_spin_system,...
                    v_off_S,spin_map,use_dipolar,...
                    term_map,euler_angles);
                % Integration step for the Signal to account for oscillation
                if v_off_S==0
                    t9=1/(off_1*Nint);
                else
                    t9=1/(v_off_S*Nint);
                end

                % Loop over RF frequencies
                for a=1:Npts_EN

                    % Radiofrequency
                    v_RF=(start_EN+step_EN*(a-1));

                    Hcorr=zeros(size(Hfree_p));
                    HRF=Hfree_p;

                    if parameters.Bterm==false
                        Hcorr=2*pi*v_RF*ops.Iz_rf;
                        HRF=Hfree_p+Hcorr+oneN*ops.Iy_rf;
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

function endor_amp=kehl_mims_calc(spin_system,parameters)

    % Check consistency
    grumble(spin_system,parameters,[],[],[]);

    % Unpack context data
    constants=parameters.constants;
    paramsENDOR=parameters.paramsENDOR;
    EPR=parameters.epr;
    operator_spin_system=parameters.operator_spin_system;
    n_endor=parameters.n_endor;
    n_spin_systems=parameters.n_spin_systems;
    I=parameters.endor_spin_numbers;

    % Get cached operators and states
    ops=kehl_operator_basis(operator_spin_system,n_endor,n_spin_systems);
    Sx=ops.Sx;
    Sy=ops.Sy;
    Ix=ops.Ix;
    Iy=ops.Iy;
    Iz=ops.Iz;
    t=parameters.pulse_times_s;
    Nint=8;

    if n_spin_systems>1
        for spin_idx=1:n_endor-1
            Ix{spin_idx+1}=Ix{1};
            Iy{spin_idx+1}=Iy{1};
            Iz{spin_idx+1}=Iz{1};
        end
    end

    geff_sel=EPR("geff_sel");
    B_sel=EPR("B_sel");
    euler_sel=EPR("euler_sel");
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

    endor_amp=zeros(1,Npts_EN);

    if isempty(B_sel)
        % No resonance orientations were found
        return
    end

    % Loop over selected orientations
    parfor j=1:length(B_sel)
        % Select orientation-specific parameters
        geff=geff_sel(j);
        B=B_sel(j);
        euler_angles=euler_sel(j,:);

        HF_zz=HF_zz_sel(j,:);
        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);
        NQI=zeros(n_endor,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);

        S=S_sel(j);
        offsets=offsets_sel(j,:);

        endor_amp_tmp=zeros(1,Npts_EN);

        % Loop over spin-manifold offsets
        for offset_idx=1:size(offsets,2)

            v_off_S=offsets(offset_idx);
            off_1=offsets(1);

            [rho0]=kehl_rho0(constants,paramsENDOR,B,geff,operator_spin_system,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            oneE=parameters.electron_nutation;
            oneN=parameters.nuclear_nutation;

            if n_spin_systems>1
                if parameters.powder==1
                    spin_map=round(offset_idx/(2*I(1)+1));
                else
                    spin_map=offset_idx;
                end
                term_map=struct();
                use_dipolar=false;
            else
                spin_map=1:n_endor;
                term_map=struct();
                use_dipolar=true;
            end
            Hfree_p=kehl_free_ham(parameters,paramsENDOR,operator_spin_system,...
                v_off_S,spin_map,use_dipolar,...
                term_map,euler_angles);
            % Apply microwave pulses

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

            % Loop over RF frequencies
            for a=1:Npts_EN

                % Radiofrequency
                v_RF=(start_EN+step_EN*(a-1));

                Hcorr=zeros(size(Hfree_p));
                HRF=Hfree_p;

                if parameters.Bterm==false
                    Hcorr=2*pi*v_RF*ops.Iz_rf;
                    HRF=Hfree_p+Hcorr+oneN*ops.Iy_rf;
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

