% Independent dense branch reference for four-spin conditional CCE. Syntax:
%
%                         result=test_cce_exact()
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_exact()

% Specify independent full-Hilbert and subset comparisons
result=new_test_result('kernel/cce_exact','CCE exact branch validation',...
                       'Four bath spins, pure/mixed states, FID/Hahn and raw subsets.');

% Build a deterministic model without sharing the CCE Hamiltonian builder
nb=4; sys.isotopes=repmat({'13C'},1,nb+1); sys.magnet=0;
sys.output='hush'; sys.disable={'hygiene'}; sys.parallel={'local',1};
inter.zeeman.matrix=repmat({zeros(3)},1,nb+1);
inter.coupling.matrix=cell(nb+1); spin_system=create(sys,inter);
fields=[0 130 210 170 290]; tensor=cell(nb+1);
for j=1:nb
    spin_system.inter.zeeman.matrix{j+1}(3,3)=fields(j+1);
    tensor{1,j+1}=[0 0 0;0 0 0;31*j -17*j 19*j];
    spin_system.inter.coupling.matrix{1,j+1}=tensor{1,j+1}.';
    for k=j+1:nb
        tensor{j+1,k+1}=diag([11 -7 13])/(j+k);
        spin_system.inter.coupling.matrix{j+1,k+1}=tensor{j+1,k+1}.';
    end
end
p.central_spins=1; p.bath_spins=2:5; p.method='cce'; p.expansion='standard';
p.quantity='coherence'; p.order=4; p.central.alpha=[1;0]; p.central.beta=[0;1];
p.assumptions='labframe'; p.orientation=[0 0 0]; p.interaction_cutoff=Inf;
p.graph=struct('mode','complete','score','coupling_norm'); p.meanfield='none';
p.times=linspace(0,0.018,13); p.sequence=[];
p.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^29);
maximum=0; max_bound=0;
for pure=[false true]

    % Prepare normalized complex product states or maximally mixed factors
    states=repmat({eye(2)/2},1,nb);
    if pure
        for j=1:nb
            v=[cos(j/7);exp(1i*j/3)*sin(j/7)]; states{j}=v*v';
        end
    end
    p.bath_state=struct('mode','explicit_product','rho',{states});
    for hahn=[false true]
        p.sequence=[];
        if hahn
            p.sequence=struct('type','ideal','target',1,'placement','fraction',...
                              'at',0.5,'axis','x','angle',pi);
        end
        actual=cce(spin_system,p);

        % Validate every raw subset against independent Kronecker/expm operators
        for record=actual.raw{1}
            sites=record.sites; n=numel(sites); dim=2^n;
            x=[0 1;1 0]/2; y=[0 -1i;1i 0]/2; z=diag([1 -1])/2;
            local={x,y,z}; ops=cell(n,3); hb=zeros(dim); coupling=zeros(dim); rho=1;
            for j=1:n
                rho=kron(rho,states{sites(j)-1});
                for a=1:3
                    ops{j,a}=kron(eye(2^(j-1)),kron(local{a},eye(2^(n-j))));
                end
                hb=hb+fields(sites(j))*ops{j,3};
                for a=1:3
                    coupling=coupling+tensor{1,sites(j)}(3,a)*ops{j,a};
                end
            end
            for j=1:n
                for k=j+1:n
                    for a=1:3
                        for b=1:3
                            hb=hb+tensor{sites(j),sites(k)}(a,b)*ops{j,a}*ops{k,b};
                        end
                    end
                end
            end
            ha=hb+coupling/2; hb=hb-coupling/2;
            reference=zeros(numel(p.times),1);
            for k=1:numel(p.times)
                t=p.times(k);
                if hahn
                    ua=expm(-1i*hb*t/2)*expm(-1i*ha*t/2);
                    ub=expm(-1i*ha*t/2)*expm(-1i*hb*t/2);
                else
                    ua=expm(-1i*ha*t); ub=expm(-1i*hb*t);
                end
                reference(k)=trace(ua*rho*ub');
            end
            discrepancy=max(abs(record.value-reference)); maximum=max(maximum,discrepancy);
            bound=max(record.abs_error_bound); max_bound=max(max_bound,bound);
            result=test_true(result,'raw subset exact',discrepancy<=bound+dim^3*eps*100,'Independent analytic or dense reference.');
        end
        result=test_true(result,'full order bypass',actual.exact_reference,'Independent analytic or dense reference.');
    end
end
metrics.maximum_raw_error=maximum; metrics.maximum_bound=max_bound;

% Report the measured metrics
disp(metrics);

end

