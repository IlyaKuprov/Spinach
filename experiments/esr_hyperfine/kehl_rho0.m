% Initial density matrix for Kehl ENDOR kernels. Syntax:
%
%      rho0=kehl_rho0(constants,paramsENDOR,B,geff,spin_system,...
%                     parameters,HF_zz,HF_zy,HF_zx,NQI_zz)
%
% Parameters:
%
%   constants        - map containing physical constants.
%   paramsENDOR      - map containing ENDOR parameters.
%   B,geff           - magnetic field and effective g value.
%   spin_system      - full Spinach spin system.
%   parameters       - Kehl ENDOR context parameter structure.
%   HF_zz,HF_zy,HF_zx - effective hyperfine couplings.
%   NQI_zz           - effective quadrupole couplings.
%
% Outputs:
%
%   rho0             - starting density matrix.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_rho0.m>

function rho0=kehl_rho0(constants,paramsENDOR,B,geff,spin_system,...
        parameters,HF_zz,HF_zy,HF_zx,NQI_zz)
    if nargin<10
        NQI_zz=0;
    end

    % Check consistency
    grumble(constants,paramsENDOR,B,geff,spin_system,parameters,...
        HF_zz,HF_zy,HF_zx,NQI_zz);

    % Get explicit spin-system dimensions
    n_endor=parameters.n_endor;
    spin_numbers=parameters.endor_spin_numbers;

    % Get cached operators and states
    ops=kehl_operator_basis(spin_system,parameters);
    Sz=ops.Sz;
    Ix=ops.Ix;
    Iy=ops.Iy;
    Iz=ops.Iz;

    % Default quadrupole offsets to zero
    v_L=paramsENDOR("v_L");
    if NQI_zz==0
        NQI_zz=zeros(n_endor);
    end

    % Set the electron Zeeman Hamiltonian
    H_EZ=B*geff*constants("MU_B")/constants("H")*Sz;
    H_NZ=zeros(size(H_EZ));
    H_HF=zeros(size(H_EZ));
    H_NQI=zeros(size(H_EZ));

    % Add nuclear terms for the complete ENDOR spin system
    for n=1:n_endor
        H_NZ=H_NZ+v_L(n)*Iz{n};
        H_HF=H_HF+HF_zz(n)*Sz*Iz{n}+HF_zy(n)*(Sz*Iy{n})+...
            HF_zx(n)*(Sz*Ix{n});
        H_NQI=H_NQI+1/2*NQI_zz(n)*(3*Iz{n}*Iz{n}-...
            spin_numbers(n)*(spin_numbers(n)+1)*ops.eye);
    end

    % Assemble the static Hamiltonian
    H_S=H_EZ+H_NZ+H_HF+H_NQI;

    % Calculate the density matrix
    if parameters.temp_eff==true

        % Calculate the Boltzmann factor
        Boltz=expm(-H_S/(constants("K_B")*parameters.T));

        % Normalise the temperature-dependent density matrix
        rho0=Boltz/trace(Boltz);
    else

        % Use the cached electron Lz state without temperature correction
        rho0=ops.rho_z;
    end

    % Preserve only populations
    rho0=diag(diag(rho0));

end

% Consistency enforcement
function grumble(constants,paramsENDOR,B,geff,spin_system,...
        parameters,HF_zz,HF_zy,HF_zx,NQI_zz)
    if ~isa(constants,'containers.Map')
        error('constants must be a containers.Map object.');
    end
    if ~isa(paramsENDOR,'containers.Map')
        error('paramsENDOR must be a containers.Map object.');
    end
    if ~isnumeric(B)
        error('B must be numeric.');
    end
    if ~isnumeric(geff)
        error('geff must be numeric.');
    end
    if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||...
            (~isfield(spin_system,'comp'))
        error('spin_system must be a Spinach spin system structure.');
    end
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
    if ~isnumeric(HF_zz)
        error('HF_zz must be numeric.');
    end
    if ~isnumeric(HF_zy)
        error('HF_zy must be numeric.');
    end
    if ~isnumeric(HF_zx)
        error('HF_zx must be numeric.');
    end
    if ~isnumeric(NQI_zz)
        error('NQI_zz must be numeric.');
    end
end

