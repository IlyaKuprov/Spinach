% generates the relaxation superoperator for T1 relaxation in
%
% input parameters:
% O: the Relaxation inducing operator (usually Sx or Ix)
% T2: the relaxation time in s
%
% output parameters:
% RT2: the R2 relaxation superoperator
%
% February 2024 A. Kehl (akehl@gwdg.de)


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

