% generates the relaxation superoperator for T1 relaxation in
% Liouville space
%
% input parameters:
% rho0: the equilibrium density matrix
% O: the Relaxation inducing operator (usually Sx or Ix)
% T1: the relaxation time in s
%
% output parameters:
% R1: the R1 relaxation superoperator
%
% February 2024 A. Kehl (akehl@gwdg.de)

% preallocate matrix

function R1=kehl_relax_t1(rho0,O,T1)

    % Check consistency
    grumble(rho0,O,T1);
    R1=zeros(size(O).^2);

    % equilibrium population of states are diagonal values of rho0
    population=diag(rho0);
    % epsilon gives the ratio between the equl. densities
    epsilon=(population*population'.^-1)';

    % equivalent to the epsilon factor in
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

