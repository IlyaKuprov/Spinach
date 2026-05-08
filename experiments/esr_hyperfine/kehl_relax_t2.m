% T2 relaxation superoperator for Kehl ENDOR kernels. Syntax:
%
%      RT2=kehl_relax_t2(O,T2)
%
% Parameters:
%
%   O                - relaxation-inducing operator.
%   T2               - transverse relaxation time.
%
% Outputs:
%
%   RT2              - T2 relaxation superoperator.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_relax_t2.m>

function RT2=kehl_relax_t2(O,T2)

    % Check consistency
    grumble(O,T2);
    RT2=zeros(size(O).^2);

    ind=find(hilb2liouv(O,'statevec'));

    RT2(ind,ind)=-1/T2;
    RT2=diag(diag(RT2));
end

function grumble(O,T2)
    if ~isnumeric(O)
        error('O must be numeric.');
    end
    if ~isnumeric(T2)
        error('T2 must be numeric.');
    end
end

