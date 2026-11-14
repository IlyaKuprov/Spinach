% Spin-3/2 vacancy-like centre with quadrupolar nuclear environment.
% doi:10.1103/PhysRevA.109.022603; lattice and transition geometry gate.
% This small illustrative model is not a published material reproduction.
% Syntax: result=gcce_highspin_siv()
%
% ilya.kuprov@weizmann.ac.il

function result=gcce_highspin_siv()

% Magnet and isotopes; zero electron g isolates the model energy scales
sys.magnet=0; sys.isotopes={'E4','13C','29Si'};
sys.output='hush'; sys.parallel={'local',1};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Standard Spinach coupling and zero-field-splitting tensors in Hz
inter.coupling.matrix{1,1}=diag([-1 -1 2])*1e4;
inter.coupling.matrix{1,2}=diag([1700 1300 2400]);
inter.coupling.matrix{1,3}=diag([800 900 1500]);
inter.coupling.matrix{2,3}=diag([100 100 -200]);

% Compile metadata, without a full-bath basis
spin_system=create(sys,inter);

% Generalised CCE propagates the entire central density matrix
parameters.central_spins=[1];
parameters.bath_spins=[2 3];
parameters.method='gcce'; parameters.expansion='standard';
parameters.quantity='density'; parameters.order=2;
v=[0;1;1;0]/sqrt(2); parameters.central.rho=v*v';
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf; parameters.meanfield='none';
parameters.graph=struct('mode','complete');
parameters.times=linspace(0,200e-6,31); parameters.sequence=struct([]);
parameters.bath_state=struct('mode','mixed_product');
parameters.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);

% Run full order for the small bath, retaining all central coherences
result=cce(spin_system,parameters);
fprintf('gcce_highspin_siv: trace error %.3g; minimum eigenvalue %.3g.\n',...
        max(result.diagnostics.trace_error),min(result.diagnostics.minimum_eigenvalue));
fprintf('Published T2: 1.13/1.16 ms, n: 2.38/2.05; toy not fitted.\n');

% Plot central populations without silently repairing the density matrix
populations=zeros(numel(result.times),4);
for n=1:numel(result.times)

    % Read the diagonal in the documented central tensor-product order
    populations(n,:)=real(diag(result.density(:,:,n))).';
end
kfigure(); plot(1e6*result.times,populations);
xlabel('time / microseconds'); ylabel('central populations'); kgrid();

end

