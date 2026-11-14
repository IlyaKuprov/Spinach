% Partitioned CCE for a strongly coupled electron pair in a planar layer.
% Schuetz et al., doi:10.1103/PhysRevResearch.6.043303; finite exact patch.
% The deliberately small model is not a published material reproduction.
% Syntax: result=pcce_electron_layer()
%
% ilya.kuprov@weizmann.ac.il

function result=pcce_electron_layer()

% Magnet and isotopes
sys.magnet=0.01; sys.isotopes={'E','E','E','E'};
sys.output='hush'; sys.parallel={'local',1};

% Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
inter.coordinates={[0 0 0];[20 0 0];[21 0 0];[0 30 0]};
inter.zeeman.matrix=repmat({zeros(3)},1,4);
inter.coupling.matrix=cell(4);

% Explicit illustrative hyperfine and bath coupling tensors in Hz
inter.coupling.matrix{1,2}=[0 0 0;0 0 0;15000 0 8000].';
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;-12000 0 6000].';
inter.coupling.matrix{1,4}=[0 0 0;0 0 0;4000 0 3000].';
inter.coupling.matrix{2,3}=diag([2e5 2e5 -4e5]);
inter.coupling.matrix{3,4}=diag([300 300 -600]);

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

% Compile physical metadata without a full-bath basis
spin_system=create(sys,inter);

% Conditional central-spin Hahn coherence and a complete small-bath graph
parameters.central_spins=1; parameters.bath_spins=2:4;
parameters.expansion='standard'; parameters.method='cce'; parameters.quantity='coherence'; parameters.order=3;
parameters.central.alpha=[1;0]; parameters.central.beta=[0;1];
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf; parameters.meanfield='none';
parameters.graph=struct('mode','complete');
parameters.times=linspace(0,200e-6,41);
parameters.sequence=struct('type','ideal','target',1,'axis','x',...
                           'angle',pi,'placement','fraction','at',0.5);
parameters.bath_state=struct('mode','mixed_product');
parameters.numerics=struct('atol',1e-5,'rtol',1e-5,'max_memory_bytes',2^28);

% Freeze the strongly coupled pair as one partition
parameters.expansion='partition'; parameters.groups={[2 3],4};
parameters.partition=struct('mode','user','max_dimension',4);
parameters.order=2;

% Run the cluster expansion and report its physical input boundary
result=cce(spin_system,parameters);
fprintf('pcce_electron_layer: %d/%d valid values; finite-cluster demonstration.\n',...
        nnz(result.valid),numel(result.valid));

% Compare with a separate unpartitioned full-order expansion
full_parameters=rmfield(parameters,{'partition','groups'});
full_parameters.expansion='standard'; full_parameters.order=3;
full_result=cce(spin_system,full_parameters);
result.reference_error=max(abs(result.value-full_result.value));
fprintf('Partition/full-order maximum difference: %.3g (same Hamiltonian).\n',...
        result.reference_error);

% Independently assemble conditional bath Hamiltonians from input Hz tensors
bath=parameters.bath_spins; dimensions=spin_system.comp.mults(bath);
dimension=prod(dimensions); operators=cell(numel(bath),3);
for j=1:numel(bath)

    % Angular momentum matrices in descending magnetic-quantum-number order
    q=(dimensions(j)-1)/2; m=q:-1:-q;
    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
    for a=1:3

        % Embed a physical bath operator with dense Kronecker products
        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
            kron(local{a},eye(prod(dimensions(j+1:end)))));
    end
end
bath_h=zeros(dimension); difference=zeros(dimension);
for j=1:numel(bath)

    % Nuclear Zeeman terms; the illustrative electron g tensors are zero
    id=bath(j);
    if sys.isotopes{id}(1)~='E'
        bath_h=bath_h-spin(sys.isotopes{id})*sys.magnet*operators{j,3};
    end
    for a=1:3

        % Conditional hyperfine difference and on-site quadratic terms
        if ~isempty(inter.coupling.matrix{1,id})
            difference=difference+2*pi*inter.coupling.matrix{1,id}(a,3)*operators{j,a};
        end
        for b=1:3
            if ~isempty(inter.coupling.matrix{id,id})
                bath_h=bath_h+2*pi*inter.coupling.matrix{id,id}(a,b)*operators{j,a}*operators{j,b};
            end
        end
    end
    for k=j+1:numel(bath)
        for a=1:3
            for b=1:3

                % Explicit pair tensors, with no call to the CCE Hamiltonian builder
                tensor=inter.coupling.matrix{id,bath(k)};
                if ~isempty(tensor)
                    bath_h=bath_h+2*pi*tensor(a,b)*operators{j,a}*operators{k,b};
                end
            end
        end
    end
end
ha=bath_h+difference/2; hb=bath_h-difference/2;
reference=complex(zeros(numel(result.times),1));
for n=1:numel(result.times)

    % Independent full-bath Hahn propagators for a maximally mixed state
    ua=expm(-1i*ha*result.times(n)/2); ub=expm(-1i*hb*result.times(n)/2);
    reference(n)=trace(ub*ua*(ua*ub)')/dimension;
end
result.reference_error=max(abs(result.value-reference));
fprintf('Independent dense Hahn maximum complex error: %.3g.\n',result.reference_error);

% Plot the conditional complex coherence
kfigure(); plot(1e6*result.times,real(result.value),'b-',...
                1e6*result.times,imag(result.value),'r--');
xlabel('total echo time / microseconds'); ylabel('normalised coherence');
legend('real','imaginary'); kgrid();

end

