% Finite central spin-lock drive including the coupled nuclear drift.
% Piecewise constant drive in rad/s; independent full-Hilbert expm check.
% This small illustrative model is not a published material reproduction.
% Syntax: result=gcce_finite_drive()
%
% ilya.kuprov@weizmann.ac.il

function result=gcce_finite_drive()

% Magnet and isotopes; zero electron g isolates the model energy scales
sys.magnet=0; sys.isotopes={'E','13C'};
sys.output='hush'; sys.parallel={'local',1};
inter.zeeman.matrix=repmat({zeros(3)},1,2);
inter.coupling.matrix=cell(2);

% Standard Spinach coupling and zero-field-splitting tensors in Hz
inter.coupling.matrix{1,2}=diag([1100 900 2200]);

% Compile metadata, without a full-bath basis
spin_system=create(sys,inter);

% Generalised CCE propagates the entire central density matrix
parameters.central_spins=1;
parameters.bath_spins=2;
parameters.method='gcce'; parameters.expansion='standard';
parameters.quantity='density'; parameters.order=1;
v=[1;1]/sqrt(2); parameters.central.rho=v*v';
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf; parameters.meanfield='none';
parameters.graph=struct('mode','complete');
parameters.times=[0 linspace(20e-6,200e-6,30)]; parameters.sequence=struct([]);
parameters.bath_state=struct('mode','mixed_product');
parameters.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);

% A finite event is driven while the complete cluster Hamiltonian remains active
parameters.sequence=struct('type','finite','target',1,'axis','x',...
    'amplitude',4e4,'detuning',0,'placement','fraction','at',0,'duration',20e-6);

% Run full order for the small bath, retaining all central coherences
result=cce(spin_system,parameters);
fprintf('gcce_finite_drive: trace error %.3g; minimum eigenvalue %.3g.\n',...
        max(result.diagnostics.trace_error),min(result.diagnostics.minimum_eigenvalue));

% Independent dense two-spin propagator during drive plus free evolution
x=[0 1;1 0]/2; y=[0 -1i;1i 0]/2; z=diag([1 -1])/2;
h=2*pi*(1100*kron(x,x)+900*kron(y,y)+2200*kron(z,z));
rho0=kron(parameters.central.rho,eye(2)/2); exact=zeros(size(result.value));
for n=1:numel(parameters.times)

    % A fixed 20-microsecond drive is followed by variable free evolution
    t=parameters.times(n);
    duration=min(t,20e-6);
    u=expm(-1i*h*(t-duration))*expm(-1i*(h+4e4*kron(x,eye(2)))*duration);
    fullrho=u*rho0*u'; small=zeros(2);
    for k=1:2

        % Trace over the explicit nuclear tensor factor
        ids=k:2:4; small=small+fullrho(ids,ids);
    end
    exact(n,:)=small(:).';
end
result.reference_error=max(abs(result.value-exact),[],'all');
fprintf('Finite-drive dense expm maximum error: %.3g.\n',result.reference_error);

% Plot central populations without silently repairing the density matrix
populations=zeros(numel(result.times),2);
for n=1:numel(result.times)

    % Read the diagonal in the documented central tensor-product order
    populations(n,:)=real(diag(result.density(:,:,n))).';
end
kfigure(); plot(1e6*result.times,populations);
xlabel('time / microseconds'); ylabel('central populations'); kgrid();

end


