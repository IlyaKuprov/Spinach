% Exhaustive small-graph CCE family and physical grouping tests. Syntax:
%
%                    result=test_cce_clusters()
%
% Outputs:
%    result - Spinach regression result and explanatory checks
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_clusters()

% Initialise independent exhaustive family checks
result=new_test_result('kernel/cce_clusters','CCE graphs and subsets',...
                       'Check frontier enumeration, closure, powers and physical group costs.');
counts=0; max_error=0;
for n=0:7
    for topology=1:5

        % Construct path, cycle, complete, empty and deterministic dense graphs
        adjacency=false(n);
        for a=1:n
            for b=a+1:n
                edge=(topology==1&&b==a+1)||(topology==2&&(b==a+1||(a==1&&b==n)))||...
                     topology==3||(topology==5&&mod(a*b+a+b,3)~=0);
                adjacency(a,b)=edge; adjacency(b,a)=edge;
            end
        end
        model=struct('physical_ids',1:n+1,'central_ids',1,'bath_ids',2:n+1,...
                     'multiplicities',2*ones(1,n+1),'coordinates_A',[(0:n).' zeros(n+1,2)],...
                     'pairs',struct('ids',{},'tensor',{}));
        parameters=struct('method','cce','order',0,'expansion','standard',...
                          'graph',struct('mode','explicit','adjacency',sparse(adjacency)));
        graph=cce_graph(model,parameters);
        for order=0:n

            % Exhaust all subsets independently, then close connected seeds
            parameters.order=order; clusters=cce_clusters(model,graph,parameters);
            reference=false(1,2^n-1);
            for mask=1:2^n-1
                vertices=find(bitget(mask,1:n));
                if numel(vertices)>order, continue; end
                reachable=vertices(1); old=0;
                while numel(reachable)~=old
                    old=numel(reachable);
                    reachable=union(reachable,intersect(vertices,find(any(adjacency(reachable,:),1))));
                end
                if numel(reachable)==numel(vertices)
                    for submask=1:mask
                        if bitand(mask,submask)==submask, reference(submask)=true; end
                    end
                end
            end
            observed=false(size(reference));
            for k=1:numel(clusters)
                observed(sum(2.^(clusters(k).groups-1)))=true;
            end
            assert(isequal(observed,reference),'frontier family differs from brute-force closure');
            counts=counts+1;

            % Verify each Mobius incidence identity using independent bitmasks
            for k=1:numel(clusters)
                superset=false(1,numel(clusters));
                for j=1:numel(clusters)
                    superset(j)=all(ismember(clusters(k).groups,clusters(j).groups));
                end
                max_error=max(max_error,abs(sum(double([clusters(superset).mobius_power]))-1));
            end

        end
    end
end
result=test_true(result,'Brute-force families through seven vertices',counts==180,...
                 sprintf('%d topology/order cases exactly equal exhaustive subset closure',counts));
result=test_close(result,'Exact Mobius incidence identities',max_error,0,0,0,...
                  'signed coefficients sum to one over every included superset');

% Exhaust every labelled graph through four vertices at every expansion order
all_graph_cases=0;
for n=1:4
    edges=nchoosek(1:max(n,2),2);
    if n==1, edges=zeros(0,2); end
    for graph_mask=0:2^size(edges,1)-1
        adjacency=false(n);
        for edge=1:size(edges,1)
            if bitget(graph_mask,edge)
                adjacency(edges(edge,1),edges(edge,2))=true;
                adjacency(edges(edge,2),edges(edge,1))=true;
            end
        end
        model=struct('physical_ids',1:n+1,'central_ids',1,'bath_ids',2:n+1,...
                     'multiplicities',2*ones(1,n+1));
        graph=struct('bath_ids',2:n+1,'adjacency',sparse(adjacency),'scores',double(adjacency));
        parameters=struct('order',0,'expansion','standard');
        for order=0:n
            parameters.order=order; clusters=cce_clusters(model,graph,parameters);
            reference=false(1,2^n-1);
            for mask=1:2^n-1
                vertices=find(bitget(mask,1:n));
                if numel(vertices)>order, continue; end
                subgraph=adjacency(vertices,vertices)|eye(numel(vertices));
                for pivot=1:numel(vertices)
                    subgraph=subgraph|(subgraph(:,pivot)&subgraph(pivot,:));
                end
                if all(subgraph,'all')
                    for subset=1:mask
                        if bitand(mask,subset)==subset, reference(subset)=true; end
                    end
                end
            end
            observed=false(size(reference));
            for k=1:numel(clusters)
                observed(sum(2.^(clusters(k).groups-1)))=true;
            end
            assert(isequal(observed,reference),'exhaustive labelled-graph comparison failed');
            all_graph_cases=all_graph_cases+1;
        end
    end
end
result=test_true(result,'Every labelled graph through four vertices',all_graph_cases==360,...
                 '360 graph/order cases exactly equal independent Boolean-transitive-closure enumeration');

% Check named families against explicit signed coefficient references
model=struct('physical_ids',1:4,'central_ids',1,'bath_ids',2:4,...
             'multiplicities',[2 2 2 2]);
parameters=struct('order',2,'expansion','standard');
expected={int64([0 -1 0 1 1]),int64([-1 -1 -1 1 1 1]),int64([0 0 1 1])};
for topology=1:3
    adjacency=logical([0 1 0;1 0 1;0 1 0]);
    if topology==2, adjacency=~eye(3); end
    if topology==3, adjacency(2,3)=false; adjacency(3,2)=false; end
    graph=struct('bath_ids',2:4,'adjacency',sparse(adjacency),'scores',double(adjacency));
    named=cce_clusters(model,graph,parameters);
    result=test_true(result,['Named family powers ' num2str(topology)],...
                     isequal([named.mobius_power],expected{topology}),...
                     'chain, triangle and disconnected pair-plus-singleton coefficients are exact');
end

% Check cap selection, mandatory closure and physical composite dimensions
model=struct('physical_ids',[10 20 30 40 50],'central_ids',10,'bath_ids',[20 30 40 50],...
             'multiplicities',[3 2 3 2 2],'coordinates_A',zeros(5,3),...
             'pairs',struct('ids',{},'tensor',{}));
parameters=struct('method','cce','order',3,'expansion','partition',...
                  'groups',{{[20 30],40,50}},'clusters',{{[1 2 3]}},...
                  'cluster_cap',[0 0 1],'graph',struct('mode','complete'));
graph=cce_graph(model,parameters); clusters=cce_clusters(model,graph,parameters);
result=test_true(result,'Explicit closure and cap overrun',numel(clusters)==7&&...
                 isequal(clusters(1).cap_overrun,[3 3 0])&&...
                 isequal([clusters.mobius_power],int64([0 0 0 0 0 0 1])),...
                 'one capped triple adds every subset without silently violating the reported cap');
result=test_close(result,'Composite physical Hilbert dimension',clusters(end).dimension,72,0,0,...
                  'central spin one and bath multiplicities 2*3*2*2 are retained exactly');

% Measure central-mediated edges without bath-bath interactions
model.pairs=struct('ids',{[10 20],[10 40]},'tensor',{eye(3),2*eye(3)});
parameters=struct('method','gcce','order',2,'expansion','standard',...
                  'graph',struct('mode','coupling','threshold',Inf));
graph=cce_graph(model,parameters);
result=test_true(result,'Conservative central-mediated edge',graph.adjacency(1,3)&&...
                 nnz(graph.adjacency)==2&&nnz(graph.pair_norm)==0,...
                 'gCCE retains bath correlation paths through the central register at J12=0');
parameters.graph.approximate=true; approximate=cce_graph(model,parameters);
result=test_true(result,'Explicit approximate graph',nnz(approximate.adjacency)==0&&approximate.approximate,...
                 'dropping central-mediated correlations requires an explicit approximation');

% Independent full-Hilbert counterexample to bath-only gCCE connectivity
plus=[0 1;0 0]; minus=plus'; identity=eye(2); up=[1;0]; down=[0;1];
hamiltonian=kron(kron(plus,minus),identity)+kron(kron(minus,plus),identity)+...
            kron(kron(plus,identity),minus)+kron(kron(minus,identity),plus);
initial=kron(kron(up,down),down); time=0.7;
propagator=expm(-1i*time*hamiltonian); evolved=propagator*initial;
exact=sum(abs(evolved(1:4)).^2); pairwise=cos(time)^4;
result=test_close(result,'Central-mediated full-Hilbert reference',exact,cos(sqrt(2)*time)^2,...
                  32*eps,0,'two uncoupled bath spins exchange through one central spin');
result=test_true(result,'Central mediation is not disconnected factorisation',abs(exact-pairwise)>0.04,...
                 sprintf('exact %.12f versus disconnected-singleton %.12f',exact,pairwise));

% Independently check the actual two-spin interaction norm and distance cutoffs
model.pairs=struct('ids',[20 40],'tensor',4*diag([0 0 1]));
parameters.method='cce'; parameters.graph=struct('mode','coupling','threshold',1);
graph=cce_graph(model,parameters);
result=test_close(result,'Physical operator norm',graph.pair_norm(1,3),1,8*eps,0,...
                  '4*Iz*Iz has spectral norm one for two physical spin halves');
parameters.interaction_cutoff=0; other=cce_graph(model,parameters);
result=test_true(result,'Hamiltonian cutoff independence',isequal(graph,other),...
                 'graph selection never reads the physical Hamiltonian cutoff');
parameters.graph.score='resonance'; parameters.graph.mismatch=3*ones(4);
parameters.graph.mismatch_scale=4; graph=cce_graph(model,parameters);
result=test_close(result,'Explicit resonance score',graph.scores(1,3),0.8,8*eps,0,...
                  'operator norm divided by hypot(1,mismatch/scale)');
parameters.graph=struct('mode','distance','radius',1);
model.coordinates_A=[0 0 0;1 0 0;2 0 0;3 0 0;4 0 0];
parameters.bath_radius=2; graph=cce_graph(model,parameters);
result=test_true(result,'Union-of-centres bath radius',isequal(graph.bath_ids,[20 30])&&...
                 isequal(graph.omitted_bath_ids,[40 50])&&nnz(graph.adjacency)==2,...
                 'bath filtering and the separate physical graph radius preserve ID ordering');

% Supplied scores set deterministic cap selection before closure
model=struct('physical_ids',1:5,'central_ids',1,'bath_ids',2:5,...
             'multiplicities',2*ones(1,5),'coordinates_A',zeros(5,3),...
             'pairs',struct('ids',{},'tensor',{}));
parameters=struct('method','cce','order',2,'expansion','standard',...
                  'cluster_cap',[0 1],'graph',struct('mode','complete',...
                  'score','supplied','scores',[0 1 1 1;1 0 7 1;1 7 0 1;1 1 1 0]));
graph=cce_graph(model,parameters); capped=cce_clusters(model,graph,parameters);
result=test_true(result,'Deterministic score cap',isequal(capped(end).groups,[2 3])&&...
                 isequal(capped(1).cap_overrun,[2 0]),...
                 'highest supplied edge wins and its two singleton closure additions are reported');
parameters.graph.scores=ones(4)-eye(4); graph=cce_graph(model,parameters);
capped=cce_clusters(model,graph,parameters);
result=test_true(result,'Lexicographic score ties',isequal(capped(end).groups,[1 2]),...
                 'equal scores produce a unique group-ID tie order');

% A single composite group remains exact without replacing it by a pseudo-spin
parameters=struct('method','cce','order',1,'expansion','disjoint',...
                  'groups',{{2:5}},'graph',struct('mode','complete'));
graph=cce_graph(model,parameters); grouped=cce_clusters(model,graph,parameters);
result=test_true(result,'One-group exact physical cost',isscalar(grouped)&&...
                 grouped.dimension==32&&grouped.bath_dimension==16&&...
                 isequal(grouped.sites,2:5),...
                 'group order one keeps all four spin-half bath degrees of freedom');

% One hundred sparse bath vertices do not trigger a global powerset
large=100; model=struct('physical_ids',1:large+1,'central_ids',1,...
                       'bath_ids',2:large+1,'multiplicities',2*ones(1,large+1));
graph=struct('bath_ids',2:large+1,'adjacency',sparse(false(large)),...
             'scores',sparse(large,large));
parameters=struct('order',1,'expansion','standard','cluster_cap',1);
large_family=cce_clusters(model,graph,parameters);
result=test_true(result,'Sparse large-bath capped enumeration',isscalar(large_family)&&...
                 isequal(large_family.sites,2),...
                 'frontier generation selects one singleton without allocating 2^N subsets');

% Restore a small physical graph for input-refusal tests
model=struct('physical_ids',1:3,'central_ids',1,'bath_ids',2:3,...
             'multiplicities',[2 2 2]);
graph=struct('bath_ids',2:3,'adjacency',sparse(logical([0 1;1 0])),...
             'scores',[0 1;1 0]);

% Reject duplicate explicit cluster declarations rather than silently deduplicating
parameters=struct('method','cce','order',2,'expansion','standard',...
                  'clusters',{{[1 2],[1 2]}}); rejected=false;
try
    cce_clusters(model,graph,parameters);
catch exception
    rejected=strcmp(exception.identifier,'cce_clusters:explicit');
end
result=test_true(result,'Duplicate explicit-family refusal',rejected,...
                 'duplicate user clusters are invalid input, not graph-generation duplicates');

end
