% Test: kernel/cce_exact | CCE exact spin-manifold validation
% Independent full-density receiver reference for four-spin projected CCE.
%
% Syntax:
%
%                         result=test_cce_exact()
%
% Outputs:
%
%     result - standard registered-test result structure
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_exact()

% Specify independent full-Hilbert and subset comparisons
result=new_test_result('kernel/cce_exact','CCE exact fixed-receiver validation',...
                       'Four bath spins, pure/mixed states, FID/Hahn and raw subsets.');

% Build a deterministic model without sharing the CCE Hamiltonian builder
nb=4;
sys.isotopes=repmat({'13C'},1,nb+1);
sys.magnet=0;
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'local',1};
sys.tols.cce_atol=1e-8;
sys.tols.cce_rtol=1e-8;
sys.tols.cce_mem=2^29;
inter.zeeman.matrix=repmat({zeros(3)},1,nb+1);
inter.coupling.matrix=cell(nb+1);
spin_system=create(sys,inter);
fields=[0 130 210 170 290];
tensor=cell(nb+1);
for j=1:nb

    % Set each bath field and observed-spin coupling
    spin_system.inter.zeeman.matrix{j+1}(3,3)=fields(j+1);
    tensor{1,j+1}=[0 0 0;0 0 0;31*j -17*j 19*j];
    spin_system=cce_test_write(spin_system,[1 j+1],tensor{1,j+1}.');
    for k=j+1:nb

        % Set the bath-pair interaction tensor
        tensor{j+1,k+1}=diag([11 -7 13])/(j+k);
        spin_system=cce_test_write(spin_system,[j+1 k+1],tensor{j+1,k+1}.');

    end

end
p.observed=1;
p.cce.bath_spins=2:5;
p.cce.approximation='projected';
p.cce.order=4;
p.transition.alpha=[1;0];
p.transition.beta=[0;1];
p.orientation=[0 0 0];
p.cce.interaction_cutoff=Inf;
p.cce.meanfield='none';
p.times=linspace(0,0.018,13);
maximum=0;
max_bound=0;
for pure=[false true]

    % Prepare normalised complex product states or maximally mixed factors
    states=repmat({eye(2)/2},1,nb);
    if pure
        for j=1:nb

            % Prepare each pure bath-spin density
            v=[cos(j/7);exp(1i*j/3)*sin(j/7)];
            states{j}=v*v';

        end
    end
    p.bath_state=states;
    for hahn=[false true]

        % Select free evolution or the Hahn experiment
        pulse_sequence=@free_evolution;
        if hahn
            pulse_sequence=@hahn_echo;
        end
        p.transition.spins=p.observed;
        [value,actual]=clexpn(spin_system,pulse_sequence,p,'labframe');

        % Validate every raw subset against independent Kronecker/expm operators
        for record=actual.raw{1}

            % Build an independent reference for each retained subset
            sites=[];
            if record.cluster_id>0
                sites=actual.clusters(record.cluster_id).sites;
            end
            n=numel(sites);
            dim=2^n;
            x=[0 1;1 0]/2;
            y=[0 -1i;1i 0]/2;
            z=diag([1 -1])/2;
            local={x,y,z};
            ops=cell(n,3);
            hb=zeros(dim);
            coupling=zeros(dim);
            rho=1;
            for j=1:n

                % Embed each bath spin in the retained tensor product
                rho=kron(rho,states{sites(j)-1});
                for a=1:3

                    % Embed each Cartesian spin component
                    ops{j,a}=kron(eye(2^(j-1)),kron(local{a},eye(2^(n-j))));

                end
                hb=hb+fields(sites(j))*ops{j,3};
                for a=1:3

                    % Accumulate the observed-spin interaction
                    coupling=coupling+tensor{1,sites(j)}(3,a)*ops{j,a};

                end

            end
            for j=1:n

                % Accumulate the bath-pair interactions
                for k=j+1:n

                    % Select the second bath spin
                    for a=1:3

                        % Select the first Cartesian component
                        for b=1:3

                            % Accumulate each Cartesian coupling term
                            hb=hb+tensor{sites(j),sites(k)}(a,b)*ops{j,a}*ops{k,b};

                        end

                    end

                end

            end
            ha=hb+coupling/2;
            hb=hb-coupling/2;
            initial=kron([1 1;1 1]/2,rho);
            receiver=kron([0 1;0 0],eye(dim));
            inversion=kron([0 -1i;-1i 0],eye(dim));
            reference=zeros(numel(p.times),1);
            for k=1:numel(p.times)

                % Propagate and detect the current observation time
                t=p.times(k);
                if hahn
                    ua=expm(-1i*hb*t/2)*expm(-1i*ha*t/2);
                    ub=expm(-1i*ha*t/2)*expm(-1i*hb*t/2);
                else
                    ua=expm(-1i*ha*t);
                    ub=expm(-1i*hb*t);
                end
                U=blkdiag(ua,ub);
                if hahn&&t>0
                    U=inversion*U;
                end
                reference(k)=2*trace(receiver'*U*initial*U');

            end
            discrepancy=max(abs(record.value-reference));
            maximum=max(maximum,discrepancy);
            bound=max(record.abs_error_bound);
            max_bound=max(max_bound,bound);
            result=test_true(result,'raw subset exact',discrepancy<=bound+dim^3*eps*100,'Independent analytic or dense reference.');

        end
        result=test_true(result,'full order bypass',isequal(value,actual.raw{1}(end).value),'Independent analytic or dense reference.');

    end

end
metrics.maximum_raw_error=maximum;
metrics.maximum_bound=max_bound;

% Report the measured metrics
disp(metrics);

end


