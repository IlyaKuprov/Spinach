% Spin-system metadata for the Kehl ENDOR context. Syntax:
%
%      parameters=kehl_context_spin_data(spin_system,parameters)
%
% Parameters:
%
%   spin_system      - Spinach spin system structure.
%   parameters       - Kehl ENDOR context parameter structure.
%
% Outputs:
%
%   parameters       - parameter structure with Kehl spin metadata.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_context_spin_data.m>

function parameters=kehl_context_spin_data(spin_system,parameters)

    % Check consistency
    grumble(spin_system,parameters);

    isotopes=spin_system.comp.isotopes;
    electron_idx=find(cellfun(@kehl_is_electron,isotopes),1);
    if isempty(electron_idx)
        error('spin_system must contain an electron spin.');
    end

    endor_spins=parameters.endor_spins;
    epr_spins=setdiff(1:numel(isotopes),electron_idx);

    [~,electron_mult]=spin(isotopes{electron_idx});
    n_endor=numel(endor_spins);

    parameters.electron_spin=(electron_mult-1)/2;
    parameters.electron_spin_idx=electron_idx;
    parameters.electron_isotope=isotopes{electron_idx};
    parameters.electron_r2_rate=kehl_spin_r2_rate(spin_system,electron_idx);
    parameters.g_matrix=kehl_electron_g_matrix(spin_system,electron_idx);
    parameters.g_iso=trace(parameters.g_matrix)/3;
    parameters.n_endor=n_endor;
    parameters.endor_isotopes=isotopes(endor_spins);
    parameters.endor_spin_numbers=kehl_spin_numbers(isotopes(endor_spins));

    A=zeros(3*n_endor,3);
    Q=zeros(3*n_endor,3);
    Q_used=false;
    for n=1:n_endor
        spin_idx=endor_spins(n);
        A(3*n-2:3*n,:)=kehl_coupling_matrix(spin_system,electron_idx,spin_idx);
        Q_block=kehl_coupling_matrix(spin_system,spin_idx,spin_idx);
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
        CS_block=kehl_nuclear_cs_matrix(spin_system,spin_idx);
        CS(3*n-2:3*n,:)=CS_block;
        CS_used=CS_used||any(CS_block(:));
    end
    if CS_used
        parameters.cs_matrix=CS;
    end
    parameters.cs_active=CS_used;

    pairs=kehl_infer_dipolar_pairs(spin_system,endor_spins);
    if ~isempty(pairs)
        D=zeros(3*size(pairs,1),3);
        for n=1:size(pairs,1)
            D(3*n-2:3*n,:)=kehl_coupling_matrix(spin_system,pairs(n,1),pairs(n,2));
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
        parameters.epr_spins=epr_spins;
        parameters.epr_isotopes=isotopes(epr_spins);
        parameters.n_epr=n_epr;
        parameters.epr_spin_numbers=kehl_spin_numbers(isotopes(epr_spins));
        g_N_EPR=zeros(n_epr,1);
        A_EPR=zeros(3*n_epr,3);
        Q_EPR=zeros(3*n_epr,3);
        EPR_Q_used=false;
        for n=1:n_epr
            spin_idx=epr_spins(n);

            % Store Spinach magnetogyric ratio in rad/s/T
            [gamma,~]=spin(isotopes{spin_idx});
            g_N_EPR(n)=gamma;
            A_EPR(3*n-2:3*n,:)=kehl_coupling_matrix(spin_system,electron_idx,spin_idx);
            Q_block=kehl_coupling_matrix(spin_system,spin_idx,spin_idx);
            Q_EPR(3*n-2:3*n,:)=Q_block;
            EPR_Q_used=EPR_Q_used||any(Q_block(:));
        end
        parameters.epr_gamma=g_N_EPR;
        parameters.epr_hfc_matrix=A_EPR;
        parameters.epr_nqi_matrix=Q_EPR;
        parameters.epr_nqi_active=EPR_Q_used;
    else
        parameters.epr_spins=[];
        parameters.epr_isotopes={};
        parameters.n_epr=0;
        parameters.epr_spin_numbers=zeros(0,1);
        parameters.epr_gamma=zeros(0,1);
        parameters.epr_hfc_matrix=zeros(0,3);
        parameters.epr_nqi_matrix=zeros(0,3);
        parameters.epr_nqi_active=false;
    end
end

% Read scalar R2 rate from Spinach relaxation data
function rate=kehl_spin_r2_rate(spin_system,spin_idx)

    % Default to no line-width threshold without scalar R2 data
    rate=0;

    % Use canonical Spinach R2 rates when present
    if isfield(spin_system,'rlx')&&isfield(spin_system.rlx,'r2_rates')&&...
            (numel(spin_system.rlx.r2_rates)>=spin_idx)&&...
            isnumeric(spin_system.rlx.r2_rates{spin_idx})&&...
            isscalar(spin_system.rlx.r2_rates{spin_idx})
        rate=spin_system.rlx.r2_rates{spin_idx};
    end

end

% Consistency enforcement
function grumble(spin_system,parameters)
    if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||...
            (~isfield(spin_system,'comp'))
        error('spin_system must be a Spinach spin system structure.');
    end
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
end

