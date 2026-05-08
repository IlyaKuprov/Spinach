%KEHL_RHO0 Initial density matrix for Kehl ENDOR kernels.
%
%   RHO0=KEHL_RHO0(CONSTANTS,PARAMSENDOR,B,GEFF,OPERATOR_SPIN_SYSTEM,
%   PARAMETERS,HF_ZZ,HF_ZY,HF_ZX,NQI_ZZ) builds the initial density matrix
%   either from a Boltzmann factor or from the electron Lz state.
%
%   Inputs:
%
%      CONSTANTS            - map containing physical constants.
%      PARAMSENDOR          - map containing ENDOR parameters.
%      B,GEFF               - magnetic field and effective g value.
%      OPERATOR_SPIN_SYSTEM - reduced Spinach spin system for Kehl kernels.
%      PARAMETERS           - Kehl ENDOR context parameters.
%      HF_ZZ,HF_ZY,HF_ZX    - effective hyperfine couplings.
%      NQI_ZZ               - effective quadrupole couplings.
%
%   Output:
%
%      RHO0 - starting density matrix.
%
%   February 2024 A. Kehl (akehl@gwdg.de)
%   Spinach architecture migration May 2026 Talos

function rho0=kehl_rho0(constants,paramsENDOR,B,geff,operator_spin_system,...
                       parameters,HF_zz,HF_zy,HF_zx,NQI_zz)
if nargin<10
    NQI_zz=0;
end

% Check consistency
grumble(constants,paramsENDOR,B,geff,operator_spin_system,parameters,...
        HF_zz,HF_zy,HF_zx,NQI_zz);

% Get explicit spin-system dimensions
n_endor=parameters.n_endor;
n_spin_systems=parameters.n_spin_systems;
spin_numbers=parameters.endor_spin_numbers;

% Get cached operators and states
ops=kehl_operator_basis(operator_spin_system,n_endor,n_spin_systems);
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

% Add nuclear terms for separated-spin-system mode
if n_spin_systems>1
    m=1;
    H_NZ=H_NZ+v_L(m)*Iz{m};
    H_HF=H_HF+HF_zz(m)*Sz*Iz{m}+HF_zy(m)*(Sz*Iy{m})+...
         HF_zx(m)*(Sz*Ix{m});
    H_NQI=H_NQI+1/2*NQI_zz(m)*(3*Iz{m}*Iz{m}-...
          spin_numbers(m)*(spin_numbers(m)+1)*ops.eye);
else

    % Add nuclear terms for the complete ENDOR spin system
    for m=1:n_endor
        H_NZ=H_NZ+v_L(m)*Iz{m};
        H_HF=H_HF+HF_zz(m)*Sz*Iz{m}+HF_zy(m)*(Sz*Iy{m})+...
             HF_zx(m)*(Sz*Ix{m});
        H_NQI=H_NQI+1/2*NQI_zz(m)*(3*Iz{m}*Iz{m}-...
              spin_numbers(m)*(spin_numbers(m)+1)*ops.eye);
    end
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
function grumble(constants,paramsENDOR,B,geff,operator_spin_system,...
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
if (~isstruct(operator_spin_system))||(~isfield(operator_spin_system,'bas'))||...
   (~isfield(operator_spin_system,'comp'))
    error('operator_spin_system must be a Spinach spin system structure.');
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
