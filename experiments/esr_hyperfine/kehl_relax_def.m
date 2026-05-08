% adds the relaxation parameters to the spinsystem
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
% February 2024 A. Kehl (akehl@gwdg.de)


function spinSys_out=kehl_relax_def(spinSys_in,T1e,T1n,T2e,T2n,T2dq)

    % Check consistency
    grumble(spinSys_in,T1e,T1n,T2e,T2n,T2dq);
    spinSys_out=spinSys_in;

    spinSys_out('T1e')=T1e;
    spinSys_out('T1n')=T1n;
    spinSys_out('T2e')=T2e;
    spinSys_out('T2n')=T2n;
    spinSys_out('T2dq')=T2dq;

end

function grumble(spinSys_in,T1e,T1n,T2e,T2n,T2dq)
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
end

