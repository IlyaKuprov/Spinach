% Cross-polarisation ENDOR pulse sequence for the Kehl ENDOR context. Syntax:
%
%      endor_amp=endor_kehl_cp(spin_system,parameters,H,R,K)
%
% Parameters:
%
%   spin_system      - Spinach spin system structure.
%   parameters       - Kehl ENDOR context parameter structure.
%   H,R,K            - Spinach experiment signature matrices.
%
% Outputs:
%
%   endor_amp        - simulated cross-polarisation ENDOR amplitude array.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=endor_kehl_cp.m>

function endor_amp=endor_kehl_cp(spin_system,parameters,H,R,K)

    % Append sequence-specific parameters when requested by the context
    if nargin>=3&&ischar(H)&&strcmp(H,'parameters')
        endor_amp=kehl_cp_parameters(spin_system,parameters);
        return
    end
    % Check consistency
    grumble(spin_system,parameters,H,R,K);
    endor_amp=kehl_cp_calc(spin_system,parameters,R);
end

function parameters=kehl_cp_parameters(spin_system,parameters)

    % Check CP sweep inputs
    if ~isfield(parameters,'cp_start_hz')
        error('parameters.cp_start_mhz must be specified for CP ENDOR.');
    end
    if ~isfield(parameters,'cp_npoints')
        error('parameters.cp_npoints must be specified for CP ENDOR.');
    end
    if ~isfield(parameters,'cp_range_hz')
        error('parameters.cp_range_mhz must be specified for CP ENDOR.');
    end

    % Get pulse sequence timing and RF-field policy
    constants=parameters.constants;
    t=parameters.pulse_times_s;
    [rf_nutations,rf_auto]=kehl_rf_policy(parameters);

    % Set CP RF nutation fields
    if rf_auto==false
        if size(rf_nutations,2)==5
            parameters.prep_nutation=rf_nutations(1)*2*pi*1e6;
            parameters.spinlock_nutation=rf_nutations(2)*2*pi*1e6;
            parameters.cp_nutation=rf_nutations(3)*2*pi*1e3;
            parameters.electron_nutation=rf_nutations(4)*2*pi*1e6;
            parameters.nuclear_nutation=rf_nutations(5)*2*pi*1e3;
            parameters.pulse_width=parameters.prep_nutation/...
                (2*pi*constants('CONST1')*1e10);
        elseif size(rf_nutations,2)==3
            parameters.prep_nutation=rf_nutations(1)*2*pi*1e6;
            parameters.spinlock_nutation=rf_nutations(2)*2*pi*1e6;
            parameters.cp_nutation=rf_nutations(3)*2*pi*1e3;
            parameters.electron_nutation=2*pi/(4*t(7));
            parameters.nuclear_nutation=2*pi/(2*t(5));
            parameters.pulse_width=parameters.prep_nutation/...
                (2*pi*constants('CONST1')*1e10);
        elseif size(rf_nutations,2)==2
            if t(1)==0
                parameters.prep_nutation=2*pi/(4*t(7));
                parameters.spinlock_nutation=rf_nutations(1)*2*pi*1e6;
                parameters.cp_nutation=rf_nutations(2)*2*pi*1e3;
                parameters.pulse_width=parameters.spinlock_nutation/...
                    (2*pi*constants('CONST1')*1e10);
            else
                parameters.prep_nutation=2*pi/(4*t(1));
                parameters.spinlock_nutation=rf_nutations(1)*2*pi*1e6;
                parameters.cp_nutation=rf_nutations(2)*2*pi*1e3;
                parameters.pulse_width=parameters.prep_nutation/...
                    (2*pi*constants('CONST1')*1e10);
            end
            parameters.electron_nutation=2*pi/(4*t(7));
            parameters.nuclear_nutation=2*pi/(2*t(5));
        else
            error('parameters.rf_nutation_freqs has incompatible dimensions.');
        end
    else
        parameters.prep_nutation=2*pi/(4*t(7));
        parameters.spinlock_nutation=2*pi/(4*t(7));
        parameters.cp_nutation=2*pi/(2*t(5));
        parameters.electron_nutation=2*pi/(4*t(7));
        parameters.nuclear_nutation=2*pi/(2*t(5));
        parameters.pulse_width=parameters.electron_nutation/...
            (2*pi*constants('CONST1')*1e10);
    end

    % Append standard ENDOR sweep-axis data
    parameters=kehl_endor_axis(spin_system,parameters,'endor');

    % Append CP sweep-axis data
    if parameters.powder==false
        parameters.cp_npoints=1;
        parameters.cp_range_hz=0;
        parameters.y_coords=1;
    else
        if parameters.cp_npoints>1
            step_CP=parameters.cp_range_hz/(parameters.cp_npoints-1);
            y_coords=zeros(parameters.cp_npoints);
            for n=1:parameters.cp_npoints
                y_coords(n)=parameters.cp_start_hz+(n-1)*step_CP;
            end
            parameters.y_coords=y_coords(:,1)';
        else
            y_coords=parameters.cp_start_hz;
            parameters.y_coords=y_coords(:,1)';
        end
        parameters.paramsENDOR('y_coords')=parameters.y_coords;
    end

end

function endor_amp=kehl_cp_calc(spin_system,parameters,R)

    % Check consistency
    grumble(spin_system,parameters,[],[],[]);

    % Unpack context data
    constants=parameters.constants;
    paramsENDOR=parameters.paramsENDOR;
    EPR=parameters.epr;
    n_endor=parameters.n_endor;

    % Get cached operators and states
    ops=kehl_operator_basis(spin_system,parameters);
    Sx=ops.Sx;
    Sy=ops.Sy;
    Sy_state=ops.Sy_state;
    Ix=ops.Ix;
    Iy=ops.Iy;
    Iz=ops.Iz;
    t=parameters.pulse_times_s;
    Nint=8;
    Npts_CP=parameters.cp_npoints;

    if parameters.cp_npoints>1
        step_CP=parameters.cp_range_hz/(parameters.cp_npoints-1);
    else
        step_CP=0;
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
    for j=1:length(B_sel)
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

        % Loop over spin-manifold offsets
        for offset_idx=1:size(offsets,2)

            v_off_S=offsets(offset_idx);
            off_1=offsets(1);

            rho0=kehl_rho0(constants,paramsENDOR,B,geff,spin_system,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            % for tilted frame Rho relaxation times
            RRho=R;

            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            prep=parameters.prep_nutation;
            sl=parameters.spinlock_nutation;
            cp=parameters.cp_nutation;
            oneE=parameters.electron_nutation;
            oneN=parameters.nuclear_nutation;
            Hfree_p=kehl_free_ham(parameters,paramsENDOR,spin_system,...
                v_off_S,euler_angles);
            for mm=1:n_endor
                if parameters.powder==false
                    [v_CP,paramsENDOR]=kehl_cp_offset(parameters,paramsENDOR,...
                        v_off_S,HF_zz,NQI_zz,...
                        mm,offset_idx);
                end
            end

            % Integration step for the Signal to account for oscillation
            t11=kehl_offset_step(v_off_S,off_1,Nint);

            % Loop over CP frequencies
            for c=1:Npts_CP

                if parameters.powder==true
                    v_CP=parameters.cp_start_hz+step_CP*(c-1);
                end

                % Loop over RF frequencies
                parfor a=1:Npts_EN

                    if abs(HF_zz(1))>1/parameters.T2e*0.1

                        % Radiofrequency
                        v_RF=(start_EN+step_EN*(a-1));

                        Hcorr=zeros(size(Hfree_p));
                        % Hamiltonian for RF pulse (no HF enhancement)
                        HRF=Hfree_p;
                        HSL=Hfree_p+sl*Sy;

                        if parameters.Bterm==false
                            for mm=1:n_endor
                                HRF=HRF+2*pi*v_RF*Iz{mm}+oneN*Ix{mm};
                                HSL=HSL+2*pi*v_CP*Iz{mm}+cp*Ix{mm};
                            end
                        end

                        HSL_t=HSL;

                        for mm=1:n_endor
                            Hcorr=Hcorr+2*pi*v_CP*Iz{mm};
                        end

                        Hfree=Hfree_p+Hcorr;

                        Hnonsel=(Hfree+oneE*Sx);
                        Hprep=(Hfree+prep*Sx);

                        Hfree=Hfree;

                        if parameters.Bterm==false
                            U3=full(propagator(spin_system,1i*sparse(RRho-1i*HSL_t),t(3)));

                            U5=full(propagator(spin_system,1i*sparse(R-1i*HRF),t(5)));

                            U1=full(propagator(spin_system,1i*sparse(R-1i*Hprep),t(1)));
                            U2=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t(2)));

                            U4=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t(4)));

                            U6=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t(6)));
                            U7=full(propagator(spin_system,1i*sparse(R-1i*Hnonsel),t(7)));
                            U8=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t(8)));
                            U9=full(propagator(spin_system,1i*sparse(R-1i*Hnonsel),t(9)));
                            U10=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t(8)+t(7)/2));
                            U11=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t11));
                        else
                            t_stepCP=1/(v_CP*parameters.N_stepRF);
                            H_SL=HSL_t;
                            U_SL=eye(size(HSL_t));
                            for ll=1:parameters.N_stepRF
                                for mm=1:n_endor
                                    H_SL=H_SL+2*parameters.cp_nutation*Ix{mm}*cos(2*pi*v_CP*t_stepCP*(ll-1));
                                end
                                G=RRho/t_stepCP-1i*H_SL;
                                U_step=full(propagator(spin_system,1i*sparse(G),t_stepCP));
                                U_SL=U_step*U_SL;
                            end
                            U3=(U_SL^(t(3)*v_CP));
                            U5=kehl_rf_bterm(parameters,v_RF,Hfree_p,Iy,t(5),n_endor,spin_system,R);

                            U1=full(propagator(spin_system,1i*sparse(R-1i*Hprep),t(1)));
                            U2=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t(2)));

                            U4=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t(4)));

                            U6=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t(6)));
                            U7=full(propagator(spin_system,1i*sparse(R-1i*Hnonsel),t(7)));
                            U8=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t(8)));
                            U9=full(propagator(spin_system,1i*sparse(R-1i*Hnonsel),t(9)));
                            U10=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t(8)+t(7)/2));
                            U11=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t11));

                        end

                        % Evolve the densitymatrix
                        rho=rho0;
                        rho=U1*rho;
                        rho=U2*rho;

                        rho=U3*rho;
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
                            value_Sy=value_Sy+real(Sy_state'*rho);
                        end
                        endor_amp(a)=endor_amp(a)+(value_Sy*S/(Nint*size(offsets,2)));
                    end
                end
            end

        end
    end
end

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

