% Low-field clock-like spin-one centre with transverse ZFS and a nuclear bath.
% doi:10.1103/PRXQuantum.2.010311; source D/E and atomistic bath gate.
% This small illustrative model is not a published material reproduction.
% Syntax: result=gcce_sic_clock()
%
% ilya.kuprov@weizmann.ac.il

function result=gcce_sic_clock()

% Magnet and isotopes; zero electron g isolates the model energy scales
sys.magnet=0; sys.isotopes={'E3','13C','29Si'};
sys.output='hush'; sys.parallel={'local',1};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Standard Spinach coupling and zero-field-splitting tensors in Hz
inter.coupling.matrix{1,1}=diag([-1.0e4+2e3 -1.0e4-2e3 2e4]);
inter.coupling.matrix{1,2}=diag([2200 1800 3000]);
inter.coupling.matrix{1,3}=diag([900 1200 1700]);
inter.coupling.matrix{2,3}=diag([120 120 -240]);

% Compile metadata, without a full-bath basis
spin_system=create(sys,inter);

% Generalised CCE propagates the entire central density matrix
parameters.central_spins=1;
parameters.bath_spins=[2 3];
parameters.method='gcce'; parameters.expansion='standard';
parameters.quantity='density'; parameters.order=2;
v=[1;0;1]/sqrt(2); parameters.central.rho=v*v';
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf; parameters.meanfield='none';
parameters.graph=struct('mode','complete');
parameters.times=linspace(0,200e-6,31); parameters.sequence=struct([]);
parameters.bath_state=struct('mode','mixed_product');
parameters.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);

% Run full order for the small bath, retaining all central coherences
result=cce(spin_system,parameters);
fprintf('gcce_sic_clock: trace error %.3g; minimum eigenvalue %.3g.\n',...
        max(result.diagnostics.trace_error),min(result.diagnostics.minimum_eigenvalue));
fprintf('Clock-like toy, not the published ~1.15 ms SiC ensemble echo.\n');

% Independently construct the full dense Hamiltonian from the input Hz tensors
dimensions=spin_system.comp.mults; dimension=prod(dimensions);
operators=cell(numel(dimensions),3); h=zeros(dimension);
for j=1:numel(dimensions)

    % Local angular momentum matrices and their tensor-product embeddings
    q=(dimensions(j)-1)/2; m=q:-1:-q;
    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
    for a=1:3
        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
            kron(local{a},eye(prod(dimensions(j+1:end)))));
    end
end
for j=1:numel(dimensions)
    for k=j:numel(dimensions)

        % Spinach ordered pair matrices transpose into Cartesian operators
        tensor=inter.coupling.matrix{j,k};
        if isempty(tensor), continue; end
        for a=1:3
            for b=1:3
                h=h+2*pi*tensor(b,a)*operators{j,a}*operators{k,b};
            end
        end
    end
end

% Propagate an explicit product density and trace over the whole finite bath
d0=size(parameters.central.rho,1); db=dimension/d0;
rho0=kron(parameters.central.rho,eye(db)/db); reference=zeros(size(result.value));
for n=1:numel(result.times)

    % This independent calculation uses expm, not the CCE eigensolver
    u=expm(-1i*h*result.times(n)); rho=u*rho0*u'; reduced=zeros(d0);
    for k=1:db
        ids=k:db:dimension; reduced=reduced+rho(ids,ids);
    end
    reference(n,:)=reduced(:).';
end
result.reference_error=max(abs(result.value-reference),[],'all');
fprintf('Independent dense density maximum error: %.3g.\n',result.reference_error);

% Plot central populations without silently repairing the density matrix
populations=zeros(numel(result.times),3);
for n=1:numel(result.times)

    % Read the diagonal in the documented central tensor-product order
    populations(n,:)=real(diag(result.density(:,:,n))).';
end
kfigure(); plot(1e6*result.times,populations);
xlabel('time / microseconds'); ylabel('central populations'); kgrid();

end


