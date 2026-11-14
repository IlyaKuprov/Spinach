% Mixed 13C/29Si nuclear bath around a divacancy-like centre.
% Seo et al., doi:10.1038/ncomms12935; full 4H-SiC geometry is gated.
% The deliberately small model is not a published material reproduction.
% Syntax: result=cce_sic_divacancy()
%
% ilya.kuprov@weizmann.ac.il

function result=cce_sic_divacancy()

% Magnet and isotopes
sys.magnet=0.01; sys.isotopes={'E','13C','29Si'};
sys.output='hush'; sys.parallel={'local',1};

% Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
inter.coordinates={[0 0 0];[4 0 1];[0 5 2]};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Explicit illustrative hyperfine and bath coupling tensors in Hz
inter.coupling.matrix{1,2}=[0 0 0;0 0 0;15000 0 8000].';
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;9000 0 -4000].';
inter.coupling.matrix{2,3}=diag([250 250 -500]);

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

% Compile physical metadata without a full-bath basis
spin_system=create(sys,inter);

% Conditional central-spin Hahn coherence and a complete small-bath graph
parameters.central_spins=1; parameters.bath_spins=2:3;
parameters.expansion='standard'; parameters.method='cce'; parameters.quantity='coherence'; parameters.order=2;
parameters.central.alpha=[1;0]; parameters.central.beta=[0;1];
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf; parameters.meanfield='none';
parameters.graph=struct('mode','complete');
parameters.times=linspace(0,200e-6,41);
parameters.sequence=struct('type','ideal','target',1,'axis','x',...
                           'angle',pi,'placement','fraction','at',0.5);
parameters.bath_state=struct('mode','mixed_product');
parameters.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);


% Run the cluster expansion and report its physical input boundary
result=cce(spin_system,parameters);
fprintf('cce_sic_divacancy: %d/%d valid values; finite-cluster demonstration.\n',...
        nnz(result.valid),numel(result.valid));
fprintf('Published 65 G: T2=1.18 ms,n=2.27; 115 G: 1.26 ms,n=2.45; no toy fit.\n');

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

