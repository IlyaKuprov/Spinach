% Test: kernel/cce_reduction | CCE deterministic reduction
% Complex CCE algebra, gather ordering, and additive polynomial tests.
%
%                     result=test_cce_reduction()
%
% Outputs:
%    result - regression result with measured algebraic errors
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_reduction()

% Compile a nontrivial signed-power family with disconnected subsets
result=new_test_result('kernel/cce_reduction','CCE deterministic reduction',...
                       'Compare direct factors, signed powers and additive inclusion-exclusion.');
model.cce=struct();
model.tols=struct('cce_atol',1e-9,'cce_rtol',1e-9,'cce_mem',Inf);
model_p.observed=1;
model_p.cce.bath_spins=2:5;
model.comp.mults=2*ones(1,5);
model.comp.nspins=numel(model.comp.mults);
adjacency=logical([0 1 0 0;1 0 1 0;0 1 0 1;0 0 1 0]);
graph=struct('bath_ids',2:5,'adjacency',sparse(adjacency),'scores',double(adjacency));
parameters=struct('method','cce','order',3,'expansion','standard');
model_p.cce.order=parameters.order;
model_p.cce.expansion=parameters.expansion;
clusters=cce_clusters(model,graph,model_p);
count=numel(clusters);
time=linspace(0,2,41).';
baseline=exp(0.13i*time);
factors=complex(zeros(numel(time),count));
raw=repmat(struct('cluster_id',uint64(0),'value',baseline,...
                  'abs_error_bound',zeros(size(time)),'valid',true(size(time))),count+1,1);
for k=1:count

    % Construct synthetic irreducibles, then raw subset products independently
    factors(:,k)=exp((-0.003*k+1i*(0.2*k-0.7))*time);
    raw(k+1).cluster_id=clusters(k).id;
    raw(k+1).value=baseline.*factors(:,k);
    for j=1:k-1
        if all(ismember(clusters(j).groups,clusters(k).groups))
            raw(k+1).value=raw(k+1).value.*factors(:,j);
        end
    end

end

% Compare both algorithms to the known product of irreducibles
reference=baseline.*prod(factors,2);
powers=cce_reduce(model,raw,clusters,parameters);
parameters.reducer='direct';
direct=cce_reduce(model,raw,clusters,parameters);
result=test_true(result,'Resolved direct and powers masks',all(powers.valid,'all')&&all(direct.valid,'all'),...
                 'both independently propagated conditioning bounds fit the caller budget');
result=test_close(result,'Signed powers versus known factors',powers.value,reference,3e-13,3e-13,...
                  'complex phases and negative Mobius exponents are retained');
result=test_close(result,'Direct recursion versus signed powers',direct.value,powers.value,3e-13,3e-13,...
                  'the same closed family has two independent multiplicative formulations');
result=test_close(result,'Unit initial value',powers.by_order(1,:,:),ones(size(powers.by_order(1,:,:))),...
                  8*eps,0,'every prefix has t=0 normalisation exactly one');

% Prove the gather permutation cannot alter deterministic arithmetic
parameters.reducer='powers';
permutation=[count+1:-2:1,count:-2:1];
gathered=cce_reduce(model,raw(permutation),clusters,parameters);
result=test_true(result,'Permutation-gather bitwise identity',isequaln(gathered,powers),...
                 'cluster IDs canonicalise all gather schedules before compensated reduction');

% Check every order against the product of known lower-order irreducibles
prefix_error=0;
for order=0:parameters.order
    expected=baseline.*prod(factors(:,[clusters.order]<=order),2);
    prefix_error=max(prefix_error,max(abs(powers.by_order(:,:,order+1)-expected)));
end
result=test_close(result,'All order prefixes',prefix_error,0,3e-13,0,...
                  'powers are recomputed for each prefix rather than reused from maximum order');

% Periodic segments cancel all but the two adjacent length tiers
cycle=graph;
cycle.adjacency(1,4)=true;
cycle.adjacency(4,1)=true;
segment_par=parameters;
segment_par.expansion='periodic_segments';
segment_p=model_p;
segment_p.cce.expansion='periodic_segments';
segments=cce_clusters(model,cycle,segment_p);
segment_raw=repmat(raw(1),numel(segments)+1,1);
for k=1:numel(segments)

    % Assign unrelated complex traces so an extra baseline factor cannot cancel
    segment_raw(k+1).cluster_id=segments(k).id;
    segment_raw(k+1).value=exp((-0.01*k+0.07i*k)*time);

end
periodic=cce_reduce(model,segment_raw,segments,segment_par);
orders=[0 segments.order];
result=test_true(result,'Only periodic segments',isequal(segments(1).actual_counts,[4 4 4]),...
                 'the cycle includes the wrap-around bond without disconnected subsets');
for n=1:3

    % Compare signed powers with adjacent-order products of the unnormalised traces
    expected=prod([segment_raw(orders==n).value],2);
    if n>1
        expected=expected./prod([segment_raw(orders==n-1).value],2);
    end
    result=test_close(result,sprintf('Periodic segment order %d',n),...
                      periodic.by_order(:,:,n+1),expected,3e-13,3e-13,...
                      'each smaller segment cancels except the adjacent overlap');

end
segment_raw(find(orders==2,1)).value(2)=0;
periodic=cce_reduce(model,segment_raw,segments,segment_par);
result=test_true(result,'Unresolved periodic overlap',~periodic.valid(2)&&...
                 strcmp(periodic.reason{2},'denominator_unresolved'),...
                 'a zero length-two trace makes the length-three ratio unresolved');

% Additive quadratic set functions have known singleton and pair cumulants
parameters.method='noise';
parameters.order=4;
parameters.clusters={1:4};
model_p.cce.order=parameters.order;
model_p.cce.clusters=parameters.clusters;
clusters=cce_clusters(model,graph,model_p);
raw=repmat(raw(1),numel(clusters)+1,1);
raw(1).value=2+0.3*time;
linear=[0.2 0.4 0.6 0.8];
quadratic=reshape(1:16,4,4)/100;
for k=1:numel(clusters)
    sites=clusters(k).groups;
    value=2+0.3*time+sum(linear(sites))*time;
    for a=sites
        for b=sites(sites>a)
            value=value+quadratic(a,b)*time.^2;
        end
    end
    raw(k+1).cluster_id=clusters(k).id;
    raw(k+1).value=value;
end
additive=cce_reduce(model,raw,clusters,parameters);
result=test_close(result,'Additive polynomial full order',additive.value,raw(end).value,3e-13,3e-13,...
                  'linear Mobius coefficients recover the exact quadratic set function');
result=test_close(result,'Additive quadratic saturation',additive.by_order(:,:,3),raw(end).value,3e-13,3e-13,...
                  'all connected increments above polynomial degree two cancel');

% Empty-bath order zero has no artificial division by observed signal
empty=clusters([]);
raw=raw(1);
raw.value(:)=0;
observed=cce_reduce(model,raw,empty,parameters);
result=test_true(result,'observed-spin-only exact zero',all(observed.valid,'all')&&all(observed.value==0,'all'),...
                 'a direct baseline requires no artificial multiplicative normalisation');

end


