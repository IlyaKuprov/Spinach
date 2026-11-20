% Compares periodic acquisition with direct chronological evolution.
% Syntax: result=test_periodic_acquire()
%
% Outputs:
%
%    result - complex signal and input boundary checks
%
% talos@spindynamics.org

function result=test_periodic_acquire()

% Define a non-Hermitian receiver and two noncommuting ticks
result=new_test_result('kernel/periodic_acquire','Periodic complex acquisition',...
                       'Periodic propagation must match chronological products.');
spin_system.bas.formalism='zeeman-hilb';
spin_system.tols.liouv_zero=1e-10;
parameters.rate=2500;
parameters.npoints=21;
parameters.dwell=1/(2*parameters.rate);
parameters.start_tick=2;
parameters.rho0=[0 1; 1 0]/2;
parameters.coil=[0 1; 0 0];
Z=diag([1 -1])/2;
X=parameters.rho0;
Y=[0 -1i; 1i 0]/2;
H={2*pi*720*Z,2*pi*560*X};
parameters.tick_rf={2*pi*170*X,2*pi*440*Y};
zero=zeros(2);

% Walk the shifted drift and RF in chronological order
steps={expm(-1i*(H{2}+parameters.tick_rf{2})*parameters.dwell),...
       expm(-1i*(H{1}+parameters.tick_rf{1})*parameters.dwell)};
rho=parameters.rho0;
ref=complex(zeros(1,parameters.npoints));
for n=1:numel(ref)
    ref(n)=trace(parameters.coil'*rho);
    U=steps{mod(n-1,2)+1};
    rho=U*rho*U';
end
fid=periodic_acquire(spin_system,parameters,H,zero,zero);
result=test_close(result,'shifted combined period',fid,ref,1e-12,0,...
                  'Keep RF and drift aligned across noncommuting ticks.');

% Check signed static precession against its complex exponential
parameters.start_tick=1;
parameters.tick_rf={};
parameters.dwell=1/parameters.rate;
fid=periodic_acquire(spin_system,parameters,H(1),zero,zero);
ref=0.5*exp(-2i*pi*720*(0:parameters.npoints-1)*parameters.dwell);
result=test_close(result,'signed frequency',fid,ref,1e-12,0,...
                  'The non-Hermitian receiver resolves the sign of phase evolution.');

% Accept the one-tick boundary when the inferred count rounds down slightly
parameters.rate=6000;
parameters.dwell=1/(6000*41);
fid=periodic_acquire(spin_system,parameters,repmat(H(1),1,41),zero,zero);
ref=0.5*exp(-2i*pi*720*(0:parameters.npoints-1)*parameters.dwell);
result=test_close(result,'one-tick roundoff',fid,ref,1e-12,0,...
                  'A sub-ulp deficit in one tick must not reject a synchronous dwell.');

% Resolve the two-ulp deficit from a different rotor grid
parameters.rate=1000;
parameters.dwell=1/(1000*115);
fid=periodic_acquire(spin_system,parameters,repmat(H(1),1,115),zero,zero);
ref=0.5*exp(-2i*pi*720*(0:parameters.npoints-1)*parameters.dwell);
result=test_close(result,'accumulated one-tick roundoff',fid,ref,1e-12,0,...
                  'The count recovery must allow the combined arithmetic roundoff.');

% Restore the single-tick rotor period for the rejection checks
parameters.dwell=1/parameters.rate;

% Reject asynchronous acquisition, a nonperiodic waveform, and dissipation
for n=1:3

    % Introduce one invalid input at a time
    malformed=parameters;
    R=zero;
    if n==1
        malformed.dwell=malformed.dwell/3;
        reason='integer number of rotor ticks';
    elseif n==2
        malformed.tick_rf={zero,zero};
        reason='exactly one period';
    else
        R=eye(2);
        reason='zero relaxation and kinetics';
    end

    % Check the explicit rejection rather than accepting any error
    rejected=false;
    try
        periodic_acquire(spin_system,malformed,H(1),R,zero);
    catch err
        rejected=contains(err.message,reason);
    end
    result=test_true(result,reason,rejected,'Unsupported input must not be rounded or ignored.');

end

end


