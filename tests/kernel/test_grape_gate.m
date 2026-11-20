% Checks gate fidelity and derivatives through the normal objective path.
% Syntax: result=test_grape_gate()
%
% Talos <talos@spindynamics.org>

function result=test_grape_gate()

% Prepare two noncommuting controls with unequal dwells
sys.isotopes={'1H'};
sys.magnet=1;
bas.formalism='zeeman-hilb';
bas.approximation={'none'};
spin_system=test_spin_system(sys,struct(),bas);
S=pauli(2);
control.isotopes={'1H'};
control.channels=[1;1];
control.operators={S.x,S.y};
control.drifts={{0.7*S.z}};
control.pulse_dt=[0.04 0.08 0.06];
control.pwr_levels=1;
control.penalties={'none'};
control.p_weights=0;
control.gate_targets={exp(0.7i)*eye(2)};
control.method='lbfgs';
result=new_test_result('kernel/grape_gate','Gate fidelity','Checks unitary objective and derivatives.');

% A global phase must not change an identity target
control.drifts={{zeros(2)}};
local=optimcon(spin_system,control);
[~,fid]=grape_xy(zeros(2,3),local);
result=test_close(result,'global phase',fid(1),1,1e-12,0,'The trace magnitude discards global phase.');

% Compare the full gate against independently ordered dense exponentials
control.drifts={{0.7*S.z}};
target=expm(-1i*0.3*S.y);
control.gate_targets={target};
local=optimcon(spin_system,control);
wave=[3 1 -2;1 -4 2];
[~,fid,grad]=grape_xy(wave,local);
[~,value_only]=grape_xy(wave,local);
result=test_close(result,'value-only gate',value_only,fid,0,0,...
                  'Streaming and gradient paths must return the same objective.');
U=eye(2);
for n=1:3

    % Chronologically concatenate the slice gates
    H=0.7*S.z+wave(1,n)*S.x+wave(2,n)*S.y;
    U=expm(-1i*H*control.pulse_dt(n))*U;

end
result=test_close(result,'trace overlap',fid(1),abs(trace(target'*U)/2)^2,1e-12,0,...
                  'The full gate must retain pulse order.');
direction=sin(wave);
delta=1e-5;
[~,plus]=grape_xy(wave+delta*direction,local);
[~,minus]=grape_xy(wave-delta*direction,local);
result=test_close(result,'directional derivative',sum(grad(:,:,1).*direction,'all'),...
                  (plus(1)-minus(1))/(2*delta),1e-9,0,'The suffix and prefix give the trace derivative.');

% Compare neutral settings with the default objective
neutral=control;
neutral.dead_time=0;
neutral.bsiegert=false;
[~,neutral_fid,neutral_grad]=grape_xy(wave,optimcon(spin_system,neutral));
result=test_close(result,'neutral gate options',neutral_fid,fid,0,0,'Disabled options preserve fidelity.');
result=test_close(result,'neutral gate derivative',neutral_grad,grad,0,0,'Disabled options preserve derivatives.');

% Reject physical settings the gate worker does not propagate
options={'bsiegert','dead_time'};
settings={true,0.2};
for n=1:2

    % Exercise both setup and the public objective on stock code
    active=control;
    active.(options{n})=settings{n};
    refused=false;
    try
        [~,actual]=grape_xy(wave,optimcon(spin_system,active));
        fprintf('Accepted %s; objective change %.12g\n',options{n},actual(1)-fid(1));
    catch err
        refused=strcmp(err.identifier,['Spinach:optimcon:gate_' options{n}]);
    end
    result=test_true(result,['reject ' options{n}],refused,'Unsupported physical terms must not be ignored.');

end

% Reject an explicit state fidelity that would otherwise be ignored
control.fidelity='real';
refused=false;
try
    optimcon(spin_system,control);
catch err
    refused=strcmp(err.message,'gate targets require square fidelity.');
end
result=test_true(result,'reject real fidelity',refused,'Gate mode only supports squared trace overlap.');

% Require an actual unitary with the Hilbert dimension
control.fidelity='square';
targets={eye(3),2*eye(2),ones(2,2,2)};
for n=1:numel(targets)

    % Check wrong dimensions, scaling, and rank of the array
    control.gate_targets=targets(n);
    refused=false;
    try
        optimcon(spin_system,control);
    catch err
        refused=strcmp(err.message,'each gate target must be unitary and match the Hilbert dimension.');
    end
    result=test_true(result,sprintf('reject target %d',n),refused,...
                     'Gate targets must be square unitaries at the Hilbert dimension.');

end

end


