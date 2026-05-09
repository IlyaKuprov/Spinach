% Electron spin offsets for Kehl ENDOR kernels. Syntax:
%
%      offsets=kehl_offsets(constants,parameters,spin_system,...
%                           paramsENDOR,B,geff,HF_zz,NQI_zz)
%
% Parameters:
%
%   constants        - map containing physical constants.
%   parameters       - Kehl ENDOR context parameter structure.
%   spin_system      - full Spinach spin system.
%   paramsENDOR      - map containing ENDOR parameters.
%   B,geff           - magnetic field and effective g value.
%   HF_zz,NQI_zz     - effective hyperfine and quadrupole couplings.
%
% Outputs:
%
%   offsets          - electron spin offsets for all spin manifolds.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_offsets.m>

function offsets=kehl_offsets(constants,parameters,spin_system,...
        paramsENDOR,B,geff,HF_zz,NQI_zz)

    % Check consistency
    grumble(constants,parameters,spin_system,paramsENDOR,B,geff,...
        HF_zz,NQI_zz);

    % Get cached operators
    ops=kehl_operator_basis(spin_system,parameters);
    Sz=ops.Sz;
    Sx=ops.Sx;
    Iz=ops.Iz;

    % Get explicit context data
    n_endor=parameters.n_endor;
    spin_numbers=parameters.endor_spin_numbers;
    v_L=paramsENDOR("v_L");

    % Build the spin Hamiltonian
    H_EZ=B*geff*constants("MU_B")/constants("H")*Sz;
    H_NZ=zeros(size(Sz));
    H_HF=zeros(size(Sz));
    H_NQI=zeros(size(Sz));

    % Use all ENDOR nuclear operators in compact mode
    for n=1:n_endor
        H_NZ=H_NZ-v_L(n)*Iz{n};
        H_HF=H_HF+2*pi*HF_zz(n)*(Sz*Iz{n});
        H_NQI=H_NQI+pi*NQI_zz(n)*(3*Iz{n}*Iz{n}-...
            spin_numbers(n)*(spin_numbers(n)+1)*ops.eye);
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
    for n=1:length(V_EPR)
        for k=n+1:length(V_EPR)

            % Electron transition probability in the eigenbasis
            trans_prob_EPR(q)=abs(round((V_EPR(:,n))'*Sx*(V_EPR(:,k)),9))^2;

            % Transition frequency in Hz
            freq_EPR(q)=abs(E_EPR(n)-E_EPR(k));
            q=q+1;
        end
    end

    % Keep only allowed transitions
    freq_EPR=freq_EPR(logical(trans_prob_EPR));
    offsets=zeros(1,size(freq_EPR,2));

    % Calculate offsets around the transition centre
    if mod(size(freq_EPR,2),2)==0
        for n=1:size(freq_EPR,2)/2
            offsets(n)=-(freq_EPR(size(freq_EPR,2)+1-n)-freq_EPR(n))/2;
        end
        for n=size(freq_EPR,2)/2+1:size(freq_EPR,2)
            offsets(n)=-offsets(size(freq_EPR,2)+1-n);
        end
    else
        center=size(freq_EPR,2)/2+0.5;
        for n=1:size(freq_EPR,2)
            offsets(n)=freq_EPR(n)-freq_EPR(center);
        end
    end

    % Convert angular offsets to Hz
    offsets=offsets/(2*pi);

end

% Consistency enforcement
function grumble(constants,parameters,spin_system,paramsENDOR,B,geff,...
        HF_zz,NQI_zz)
    if ~isa(constants,'containers.Map')
        error('constants must be a containers.Map object.');
    end
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
    if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||...
            (~isfield(spin_system,'comp'))
        error('spin_system must be a Spinach spin system structure.');
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

