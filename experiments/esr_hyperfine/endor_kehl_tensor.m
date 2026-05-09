% ENDOR tensor pulse sequence for the Kehl ENDOR context. Syntax:
%
%      endor_amp=endor_kehl_tensor(spin_system,parameters,H,R,K)
%
% Parameters:
%
%   spin_system      - Spinach spin system structure.
%   parameters       - Kehl ENDOR context parameter structure.
%   H,R,K            - Spinach experiment signature matrices.
%
% Outputs:
%
%   endor_amp        - simulated tensor-ENDOR amplitude array.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=endor_kehl_tensor.m>

function endor_amp=endor_kehl_tensor(spin_system,parameters,H,R,K)

    % Append sequence-specific parameters when requested by the context
    if nargin>=3&&ischar(H)&&strcmp(H,'parameters')
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
    grumble(spin_system,parameters,[],[],[]);

    % Unpack context data
    constants=parameters.constants;
    paramsENDOR=parameters.paramsENDOR;
    EPR=parameters.epr;
    n_endor=parameters.n_endor;
    I=parameters.endor_spin_numbers;

    % Get cached operators and states
    ops=kehl_operator_basis(spin_system,parameters);
    Sz=ops.Sz;
    Iy=ops.Iy;
    Iz=ops.Iz;
    Sx_D=ops.Sx_D;
    Ix_D=ops.Ix_D;

    % Unpack context maps

    t=parameters.pulse_times_s;
    Nint=8;

    geff_sel=EPR("geff_sel");
    B_sel=EPR("B_sel");
    euler_sel=EPR("euler_sel");
    HF_zz_sel=EPR("HF_zz_sel");
    HF_zy_sel=EPR("HF_zy_sel");
    HF_zx_sel=EPR("HF_zx_sel");

    NQI_zz_sel=EPR("NQI_zz_sel");
    NQI_sel=EPR("NQI_sel");
    CS_zz_sel=EPR("CS_zz_sel");

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
        B=B_sel(j);
        euler_angles=euler_sel(j,:);
        const_R=1/size(Sz,2)*constants('GE')*B/(2*pi*constants('K_B')*parameters.T);

        HF_zz=HF_zz_sel(j,:);

        NQI=zeros(n_endor,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);

        S=S_sel(j);

        endor_amp_tmp=zeros(1,Npts_EN);

        % Loop over nuclei
        for spin_idx=1:n_endor
            mI=-I(spin_idx):1:I(spin_idx);

            for jj=mI

                v_off_S=jj*HF_zz(spin_idx);
                off_1=-I(spin_idx)*HF_zz(spin_idx);

                [rho0]=2*Sz*Iz{1};

                RT2e=kehl_relax_t2(Sx_D,parameters.T2e);
                RT2n=zeros(size(RT2e));

                for mm=1:n_endor
                    RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{mm},parameters.T2dq);
                    RT2n=RT2n+kehl_relax_t2(Ix_D{mm},parameters.T2n);
                end
                R=RT2e+RT2n;
                start_EN=paramsENDOR("start_EN");
                step_EN=paramsENDOR("step_EN");

                oneN=parameters.nuclear_nutation;
                Hfree_p=kehl_free_ham(parameters,paramsENDOR,spin_system,...
                    v_off_S,euler_angles);
                % Integration step for the Signal to account for oscillation
                t2=kehl_offset_step(v_off_S,off_1,Nint);

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
                    Hfree=full(hilb2liouv(sparse(Hfree),'comm'));

                    if parameters.Bterm==false

                        U1=full(propagator(spin_system,1i*sparse(R-1i*full(hilb2liouv(sparse(HRF),'comm'))),t(1)));
                        U2=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t2));
                    else
                        U1=kehl_rf_bterm_rlx(parameters,v_RF,Hfree_p,Iy,t(1),n_endor,R,spin_system);
                        U2=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t2));
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

function endor_amp=kehl_tensor_calc(spin_system,parameters)

    % Check consistency
    grumble(spin_system,parameters,[],[],[]);

    % Unpack context data
    constants=parameters.constants;
    paramsENDOR=parameters.paramsENDOR;
    EPR=parameters.epr;
    n_endor=parameters.n_endor;
    I=parameters.endor_spin_numbers;

    % Get cached operators and states
    ops=kehl_operator_basis(spin_system,parameters);
    Sx=ops.Sx;
    Sy=ops.Sy;
    Sz=ops.Sz;
    Ix=ops.Ix;
    Iy=ops.Iy;
    Iz=ops.Iz;

    % Unpack context maps
    t=parameters.pulse_times_s;
    Nint=8;

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
        B=B_sel(j);
        euler_angles=euler_sel(j,:);
        const_R=1/size(Sz,2)*constants('GE')*B/(2*pi*constants('K_B')*parameters.T);

        HF_zz=HF_zz_sel(j,:);

        NQI=zeros(n_endor,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);

        S=S_sel(j);

        endor_amp_tmp=zeros(1,Npts_EN);

        % Loop over nuclei
        for spin_idx=1:n_endor
            mI=-I(spin_idx):1:I(spin_idx);

            for jj=mI

                v_off_S=jj*HF_zz(spin_idx);
                off_1=-I(spin_idx)*HF_zz(spin_idx);

                [rho0]=2*Sz*Iz{spin_idx};

                start_EN=paramsENDOR("start_EN");
                step_EN=paramsENDOR("step_EN");

                oneE=parameters.electron_nutation;
                oneN=parameters.nuclear_nutation;

                nuc=spin_idx;
                Hfree_p=kehl_free_ham(parameters,paramsENDOR,spin_system,...
                    v_off_S,euler_angles);
                % Apply microwave pulses

                % Integration step for the Signal to account for oscillation
                t2=abs(kehl_offset_step(v_off_S,off_1,Nint));

                U2_p=full(propagator(spin_system,sparse(Hfree_p),t2));

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
                        U1=full(propagator(spin_system,sparse(HRF),t(1)));
                        U2=U2_p*full(propagator(spin_system,sparse(Hcorr),t2));
                    else
                        Hfree=Hfree_p;
                        U1=kehl_rf_bterm(parameters,v_RF,Hfree,Iy,t(1),n_endor,spin_system);
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

