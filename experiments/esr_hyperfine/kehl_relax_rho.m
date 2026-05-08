% adds the relaxation parameters to the spinsystem
% if a tilted frame (rho) needs t be considered, i.e. CP/SL ENDOR
% is directly called for the set up
%
% input parameters:
% spinSys_in: the Map describing the spin system
% T1e: T1e relaxation time in s
% T1n: T1n relaxation time in s
% T2e: T2e relaxation time in s
% T2n: T2n relaxation time in s
% T2dq: T2 double and zero quantum time in s
%
% output parameters:
% spinSys_out: updated Map describing the spin system
%
% April 2024 A. Kehl (akehl@gwdg.de)


function spinSys_out=kehl_relax_rho(spinSys_in,T1e,T1n,T2e,T2n,T2dq,T1eRho,T1nRho,T2eRho,T2nRho,T2dqRho)

    % Check consistency
    grumble(spinSys_in,T1e,T1n,T2e,T2n,T2dq,T1eRho,T1nRho,T2eRho,T2nRho,T2dqRho);
    spinSys_out=spinSys_in;

    spinSys_out('T1e')=T1e;
    spinSys_out('T1n')=T1n;
    spinSys_out('T2e')=T2e;
    spinSys_out('T2n')=T2n;
    spinSys_out('T2dq')=T2dq;
    spinSys_out('T1eR')=T1eRho;
    spinSys_out('T1nR')=T1nRho;
    spinSys_out('T2eR')=T2eRho;
    spinSys_out('T2nR')=T2nRho;
    spinSys_out('T2dqR')=T2dqRho;

end

function grumble(spinSys_in,T1e,T1n,T2e,T2n,T2dq,T1eRho,T1nRho,T2eRho,T2nRho,T2dqRho)
if ~isa(spinSys_in,'containers.Map')
    error('spinSys_in must be a containers.Map object.');
end
if ~isnumeric(T1e)
    error('T1e must be numeric.');
end
if ~isnumeric(T1n)
    error('T1n must be numeric.');
end
if ~isnumeric(T2e)
    error('T2e must be numeric.');
end
if ~isnumeric(T2n)
    error('T2n must be numeric.');
end
if ~isnumeric(T2dq)
    error('T2dq must be numeric.');
end
if ~isnumeric(T1eRho)
    error('T1eRho must be numeric.');
end
if ~isnumeric(T1nRho)
    error('T1nRho must be numeric.');
end
if ~isnumeric(T2eRho)
    error('T2eRho must be numeric.');
end
if ~isnumeric(T2nRho)
    error('T2nRho must be numeric.');
end
if ~isnumeric(T2dqRho)
    error('T2dqRho must be numeric.');
end
end

