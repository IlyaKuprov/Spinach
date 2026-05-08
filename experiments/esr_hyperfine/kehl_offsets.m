%KEHL_OFFSETS Electron spin offsets for Kehl ENDOR kernels.
%
%   OFFSETS=KEHL_OFFSETS(CONSTANTS,PARAMETERS,OPERATOR_SPIN_SYSTEM,
%   PARAMSENDOR,B,GEFF,HF_ZZ,NQI_ZZ) calculates the electron spin offsets
%   for all spin manifolds from hyperfine and quadrupole couplings.
%
%   Inputs:
%
%      CONSTANTS            - map containing physical constants.
%      PARAMETERS           - Kehl ENDOR context parameters.
%      OPERATOR_SPIN_SYSTEM - reduced Spinach spin system for Kehl kernels.
%      PARAMSENDOR          - map containing ENDOR parameters.
%      B,GEFF               - magnetic field and effective g value.
%      HF_ZZ,NQI_ZZ         - effective hyperfine and quadrupole couplings.
%
%   Output:
%
%      OFFSETS - offsets for all spin manifolds.
%
%   February 2024 A. Kehl (akehl@gwdg.de)
%   Spinach architecture migration May 2026 Talos

function offsets=kehl_offsets(constants,parameters,operator_spin_system,...
                              paramsENDOR,B,geff,HF_zz,NQI_zz)

% Check consistency
grumble(constants,parameters,operator_spin_system,paramsENDOR,B,geff,...
        HF_zz,NQI_zz);

% Get cached operators
ops=kehl_operator_basis(operator_spin_system,parameters.n_endor,...
                        parameters.n_spin_systems);
Sz=ops.Sz;
Sx=ops.Sx;
Ix=ops.Ix;
Iy=ops.Iy;
Iz=ops.Iz;

% Get explicit context data
n_endor=parameters.n_endor;
n_spin_systems=parameters.n_spin_systems;
spin_numbers=parameters.endor_spin_numbers;
v_L=paramsENDOR("v_L");
offsets=[];

% Loop over nuclei if they are simulated as separate spin systems
for n=1:n_spin_systems

    % Build the spin Hamiltonian
    H_EZ=B*geff*constants("MU_B")/constants("H")*Sz;
    H_NZ=zeros(size(Sz));
    H_HF=zeros(size(Sz));
    H_NQI=zeros(size(Sz));

    % Use one nuclear operator set in separated-spin-system mode
    if n_spin_systems>1
        H_NZ=H_NZ-v_L(n)*Iz{1};
        H_HF=H_HF+2*pi*HF_zz(n)*(Sz*Iz{1});
        H_NQI=H_NQI+pi*NQI_zz(n)*(3*Iz{1}*Iz{1}-...
              spin_numbers(n)*(spin_numbers(n)+1)*ops.eye);
    else

        % Use all ENDOR nuclear operators in compact mode
        for mm=1:n_endor
            H_NZ=H_NZ-v_L(mm)*Iz{mm};
            H_HF=H_HF+2*pi*HF_zz(mm)*(Sz*Iz{mm});
            H_NQI=H_NQI+pi*NQI_zz(mm)*(3*Iz{mm}*Iz{mm}-...
                  spin_numbers(mm)*(spin_numbers(mm)+1)*ops.eye);
        end
    end

    % Diagonalise the electron transition Hamiltonian
    H_S=H_EZ+H_NZ+H_HF+H_NQI;
    [V_EPR,E_EPR]=eig(H_S);
    E_EPR=real(diag(E_EPR));

    % Preallocate transition arrays
    freq_EPR=zeros(1,length(V_EPR)*length(V_EPR));
    trans_prob_EPR=zeros(1,length(V_EPR)*length(V_EPR));

    % Calculate transitions between all elements
    q=1;
    for x=1:length(V_EPR)
        for y=x+1:length(V_EPR)

            % Electron transition probability in the eigenbasis
            trans_prob_EPR(q)=abs(round((V_EPR(:,x))'*Sx*(V_EPR(:,y)),9))^2;

            % Transition frequency in Hz
            freq_EPR(q)=abs(E_EPR(x)-E_EPR(y));
            q=q+1;
        end
    end

    % Keep only allowed transitions
    freq_EPR=freq_EPR(logical(trans_prob_EPR));
    offsets_tmp=zeros(1,size(freq_EPR,2));

    % Calculate offsets around the transition centre
    if mod(size(freq_EPR,2),2)==0
        for mm=1:size(freq_EPR,2)/2
            offsets_tmp(mm)=-(freq_EPR(size(freq_EPR,2)+1-mm)-freq_EPR(mm))/2;
        end
        for mm=size(freq_EPR,2)/2+1:size(freq_EPR,2)
            offsets_tmp(mm)=-offsets_tmp(size(freq_EPR,2)+1-mm);
        end
    else
        center=size(freq_EPR,2)/2+0.5;
        for mm=1:size(freq_EPR,2)
            offsets_tmp(mm)=freq_EPR(mm)-freq_EPR(center);
        end
    end

    % Store offsets for this separated spin system
    r=size(freq_EPR,2);
    offsets((1+(n-1)*r):n*r)=offsets_tmp/(2*pi);
end

end

% Consistency enforcement
function grumble(constants,parameters,operator_spin_system,paramsENDOR,B,...
                 geff,HF_zz,NQI_zz)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if (~isstruct(operator_spin_system))||(~isfield(operator_spin_system,'bas'))||...
   (~isfield(operator_spin_system,'comp'))
    error('operator_spin_system must be a Spinach spin system structure.');
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
if ~isnumeric(HF_zz)
    error('HF_zz must be numeric.');
end
if ~isnumeric(NQI_zz)
    error('NQI_zz must be numeric.');
end
end
