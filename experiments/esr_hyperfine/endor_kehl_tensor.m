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
    endor_amp=kehl_tensor_calc(spin_system,parameters,R);
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

function endor_amp=kehl_tensor_calc(spin_system,parameters,R)

    % Check consistency
    grumble(spin_system,parameters,[],[],[]);

    % Unpack context data
    constants=parameters.constants;
    paramsENDOR=parameters.paramsENDOR;
    EPR=parameters.epr;
    n_endor=parameters.n_endor;
    I=parameters.endor_spin_numbers;

    % Request Spinach operators and states directly
    electron_idx=parameters.electron_spin_idx;
    Sz=operator(spin_system,'Lz',electron_idx);
    Iy=cell(1,n_endor);
    SzIz_state=cell(1,n_endor);
    Iy_rf=sparse(size(Sz,1),size(Sz,2));
    Iz_rf=sparse(size(Sz,1),size(Sz,2));
    for n=1:n_endor
        spin_idx=parameters.endor_spins(n);
        Iy{n}=operator(spin_system,'Ly',spin_idx);
        SzIz_state{n}=state(spin_system,{'Lz','Lz'},...
            {electron_idx,spin_idx});
        Iy_rf=Iy_rf+Iy{n};
        Iz_rf=Iz_rf+operator(spin_system,'Lz',spin_idx);
    end

    % Unpack context maps

    t=parameters.pulse_times_s;
    Nint=8;

    B_sel=EPR("B_sel");
    euler_sel=EPR("euler_sel");
    HF_zz_sel=EPR("HF_zz_sel");

    S_sel=EPR("S_sel");

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

        S=S_sel(j);

        endor_amp_tmp=zeros(1,Npts_EN);

        % Loop over nuclei
        for spin_idx=1:n_endor
            mI=-I(spin_idx):1:I(spin_idx);

            for jj=mI

                v_off_S=jj*HF_zz(spin_idx);
                off_1=-I(spin_idx)*HF_zz(spin_idx);

                rho0=2*SzIz_state{spin_idx};
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
                        Hcorr=2*pi*v_RF*Iz_rf;
                        HRF=Hfree_p+Hcorr+oneN*Iy_rf;
                    end

                    Hfree=Hfree_p+Hcorr;

                    if parameters.Bterm==false

                        U1=full(propagator(spin_system,1i*sparse(R-1i*HRF),t(1)));
                        U2=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t2));
                    else
                        U1=kehl_rf_bterm(parameters,v_RF,Hfree_p,Iy,t(1),n_endor,spin_system,R);
                        U2=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t2));
                    end

                    % Evolve the densitymatrix

                    rho=rho0;
                    rho=U1*rho;

                    value_Sy=0;
                    for b=1
                        rho=U2*rho;
                        value_Sy=value_Sy+real(SzIz_state{spin_idx}'*rho);
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

