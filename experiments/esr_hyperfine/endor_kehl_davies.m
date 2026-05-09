% Davies ENDOR pulse sequence for the Kehl ENDOR context. Syntax:
%
%      [endor_amp,parameters]=endor_kehl_davies(spin_system,parameters,H,R)
%
% Parameters:
%
%   spin_system      - Spinach spin system structure.
%   parameters       - Kehl ENDOR context parameter structure.
%   H,R              - Spinach experiment signature matrices.
%
% Outputs:
%
%   endor_amp        - simulated Davies ENDOR amplitude array.
%   parameters       - parameter structure with derived sequence data.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=endor_kehl_davies.m>

function [endor_amp,parameters]=endor_kehl_davies(spin_system,parameters,H,R)

    % Check consistency
    grumble(spin_system,parameters,H,R);

    % Append sequence-specific parameters
    parameters=kehl_davies_parameters(parameters);

    % Select sequence-dependent EPR orientations
    parameters=kehl_sequence_context(spin_system,parameters);

    % Run the sequence calculation
    endor_amp=kehl_davies_calc(spin_system,parameters,R);
end

function parameters=kehl_davies_parameters(parameters)

    % Get pulse sequence timing and RF-field policy
    t=parameters.pulse_times;
    [rf_nutations,rf_auto]=kehl_rf_policy(parameters);

    % Set Davies RF nutation fields
    if rf_auto==false
        if size(rf_nutations,1)==3
            parameters.prep_nutation=rf_nutations(1);
            parameters.electron_nutation=rf_nutations(2);
            parameters.nuclear_nutation=rf_nutations(3);
        else
            error('parameters.rf_nutations has incompatible dimensions.');
        end
    else
        parameters.prep_nutation=pi/t(1);
        parameters.electron_nutation=pi/(2*t(5));
        parameters.nuclear_nutation=pi/t(3);
    end
    if ~isfield(parameters,'excite_width')
        parameters.excite_width=parameters.electron_nutation;
    end

    % Append standard ENDOR sweep-axis data
    parameters=kehl_endor_axis(parameters,'endor');

end

function endor_amp=kehl_davies_calc(spin_system,parameters,R)

    % Check consistency
    grumble(spin_system,parameters,[],R);

    % Unpack context data
    paramsENDOR=parameters.paramsENDOR;
    EPR=parameters.epr;
    n_endor=parameters.n_endor;

    % Request Spinach operators and states directly
    electron_idx=parameters.electron_spin_idx;
    Sx=operator(spin_system,'Lx',electron_idx);
    Sy_state=state(spin_system,'Ly',electron_idx);
    Iy=cell(1,n_endor);
    Iy_rf=sparse(size(Sx,1),size(Sx,2));
    Iz_rf=sparse(size(Sx,1),size(Sx,2));
    for n=1:n_endor
        spin_idx=parameters.endor_spins(n);
        Iy{n}=operator(spin_system,'Ly',spin_idx);
        Iy_rf=Iy_rf+Iy{n};
        Iz_rf=Iz_rf+operator(spin_system,'Lz',spin_idx);
    end

    % Unpack context maps

    t=parameters.pulse_times;
    Nint=8;

    B_sel=EPR("B_sel");
    euler_sel=EPR("euler_sel");
    HF_zz_sel=EPR("HF_zz_sel");

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
        euler_angles=euler_sel(j,:);

        HF_zz=HF_zz_sel(j,:);

        S=S_sel(j);
        offsets=offsets_sel(j,:);

        endor_amp_tmp=zeros(1,Npts_EN);

        if abs(HF_zz(1))>parameters.electron_r2_rate*0.1

            % Loop over spin-manifold offsets
            for offset_idx=1:size(offsets,2)

                v_off_S=offsets(offset_idx);
                off_1=offsets(1);

                rho0=-state(spin_system,'Lz',parameters.electron_spin_idx);

                start_EN=paramsENDOR("start_EN");
                step_EN=paramsENDOR("step_EN");

                prep=parameters.prep_nutation;
                oneE=parameters.electron_nutation;
                oneN=parameters.nuclear_nutation;
                Hfree_p=kehl_free_ham(parameters,paramsENDOR,spin_system,...
                    v_off_S,euler_angles);
                % Integration step for the Signal to account for oscillation
                t9=kehl_offset_step(v_off_S,off_1,Nint);

                % Loop over RF frequencies
                for a=1:Npts_EN

                    % Radiofrequency
                    v_RF=(start_EN+step_EN*(a-1));

                    Hcorr=zeros(size(Hfree_p));
                    HRF=Hfree_p;

                    if parameters.Bterm==false
                        Hcorr=v_RF*Iz_rf;
                        HRF=Hfree_p+Hcorr+oneN*Iy_rf;
                    end

                    Hfree=Hfree_p+Hcorr;

                    Hprep=(Hfree+prep*Sx);
                    Hnonsel=(Hfree+oneE*Sx);

                    if parameters.Bterm==false

                        U3=full(propagator(spin_system,1i*sparse(R-1i*HRF),t(3)));

                        U1=full(propagator(spin_system,1i*sparse(R-1i*Hprep),t(1)));
                        U2=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t(2)));

                        U4=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t(4)));
                        U5=full(propagator(spin_system,1i*sparse(R-1i*Hnonsel),t(5)));
                        U6=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t(6)));
                        U7=full(propagator(spin_system,1i*sparse(R-1i*Hnonsel),t(7)));

                        U8=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t(6)+t(5)/2));
                        U9=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t9));
                    else
                        U3=kehl_rf_bterm(parameters,v_RF,Hfree_p,Iy,t(3),n_endor,spin_system,R);

                        U1=full(propagator(spin_system,1i*sparse(R-1i*Hprep),t(1)));
                        U2=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t(2)));

                        U4=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t(4)));
                        U5=full(propagator(spin_system,1i*sparse(R-1i*Hnonsel),t(5)));
                        U6=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t(6)));
                        U7=full(propagator(spin_system,1i*sparse(R-1i*Hnonsel),t(7)));

                        U8=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t(6)+t(5)/2));
                        U9=full(propagator(spin_system,1i*sparse(R-1i*Hfree),t9));

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

                    value_Sy=0;
                    for b=1:Nint
                        rho=U9*rho;
                        value_Sy=value_Sy+real(Sy_state'*rho);
                    end

                    endor_amp_tmp(a)=endor_amp_tmp(a)+(value_Sy*S/(Nint*size(offsets,2)));
                end
            end
        end
        endor_amp=endor_amp+endor_amp_tmp;
    end
end

function grumble(spin_system,parameters,H,R)
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
end

