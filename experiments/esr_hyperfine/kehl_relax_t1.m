% T1 relaxation superoperator for Kehl ENDOR kernels. Syntax:
%
%      R1=kehl_relax_t1(rho0,O,T1)
%
% Parameters:
%
%   rho0             - equilibrium density matrix.
%   O                - relaxation-inducing operator.
%   T1               - longitudinal relaxation time.
%
% Outputs:
%
%   R1               - T1 relaxation superoperator.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_relax_t1.m>

function R1=kehl_relax_t1(rho0,O,T1)

    % Check consistency
    grumble(rho0,O,T1);
    R1=zeros(size(O).^2);

    % equilibrium population of states are diagonal values of rho0
    population=diag(rho0);
    % Calculate the equilibrium population ratio
    epsilon=(population*population'.^-1)';

    % Use the epsilon factor from
    % Feintuch, Vega in EPR spectroscopy (2018)
    % (only upper right triangle used later)
    Thermal=epsilon./(1+epsilon)./T1;
    Thermal2=1./(1+epsilon)./T1;

    S=length(O);
    [j k]=find(triu(O));

    nj=(j+(j-1)*S);
    nk=(k+(k-1)*S);

    M=numel(j);
    % add relaxation values for all spins
    for m=1:M
        R1(nj(m),nj(m))=R1(nj(m),nj(m))-Thermal(j(m),k(m));
        R1(nk(m),nj(m))=R1(nk(m),nj(m))+Thermal(j(m),k(m));
        R1(nj(m),nk(m))=R1(nj(m),nk(m))+Thermal2(j(m),k(m));
        R1(nk(m),nk(m))=R1(nk(m),nk(m))-Thermal2(j(m),k(m));
    end
end

function grumble(rho0,O,T1)
    if ~isnumeric(rho0)
        error('rho0 must be numeric.');
    end
    if ~isnumeric(O)
        error('O must be numeric.');
    end
    if ~isnumeric(T1)
        error('T1 must be numeric.');
    end
end

