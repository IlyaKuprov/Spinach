% Tests propagator cache separation by numerical policy. Syntax:
%
%                     result=test_prop_cache()
%
% Outputs:
%
%     result - checks cleanup, storage, backend-key partitioning,
%              repeated hits, and client/worker cache agreement
%
% talos@spindynamics.org

function result=test_prop_cache()

% State the cache-policy target
result=new_test_result('kernel/prop_cache','Propagator cache policy',...
                       'cache hits must match fresh computation under the requested numerical policy.');

% Create a real shared store through Spinach
sys.magnet=1; sys.isotopes={'1H'}; sys.output='hush';
sys.enable={'prop_cache'}; sys.disable={'hygiene'};
sys.parallel={'processes',1}; sys.parprops={};
inter.zeeman.scalar={0}; spin_system=create(sys,inter);
pool=gcp('nocreate'); store=pool.ValueStore;
spin_system.tols.small_matrix=2; spin_system.tols.prop_chop=1e-3;
spin_system.tols.dense_matrix=0.5; dt=0.13;
L=spdiags((1:20)'/20,1,20,20);

% Warm cleanup policies in both orders on finite nonnormal dynamics
for ordering=1:2
    store.remove(store.keys);
    for stage=1:2
        if ordering==stage
            spin_system.sys.disable={'hygiene'};
        else
            spin_system.sys.disable={'hygiene','clean-up'};
        end
        ref_system=spin_system; ref_system.sys.enable={};
        P=propagator(spin_system,L,dt); R=propagator(ref_system,L,dt);
        result=test_close(result,'cleanup policy hit',P,R,0,0,...
                          'changing cleanup must not retrieve the other policy result');
        Q=propagator(spin_system,L,dt);
        result=test_close(result,'repeat policy hit',Q,P,0,0,...
                          'a repeated request must retain its own result');
    end
    result=test_true(result,'separate cleanup records',numel(store.keys)==2,...
                     'both cleanup policies require their own store record');
end

% Exercise storage policy changes in both orders with the same generator
spin_system.sys.disable={'hygiene'}; spin_system.tols.prop_chop=1e-10;
for policy=1:2
    for ordering=1:2
        store.remove(store.keys);
        for stage=1:2
            spin_system.tols.small_matrix=2; spin_system.tols.dense_matrix=0.5;
            if policy==1
                generator=L; spin_system.tols.dense_matrix=double(ordering==stage);
            else
                generator=full(L);
                if ordering==stage, spin_system.tols.small_matrix=20; end
            end
            ref_system=spin_system; ref_system.sys.enable={};
            P=propagator(spin_system,generator,dt); R=propagator(ref_system,generator,dt);
            result=test_close(result,'storage policy values',P,R,0,0,...
                              'storage settings must preserve the requested fresh result');
            result=test_true(result,'storage policy representation',issparse(P)==issparse(R),...
                             'cache hits must honour sparse/full storage policy');
        end
        result=test_true(result,'separate storage records',numel(store.keys)==2,...
                         'each effective storage policy must have a separate record');
    end
end

% Partition backend policy using a case that requires no GPU operations
store.remove(store.keys); spin_system.tols.small_matrix=2;
spin_system.sys.enable={'prop_cache'}; P=propagator(spin_system,L,dt);
spin_system.sys.enable={'prop_cache','gpu'}; Q=propagator(spin_system,L,dt);
result=test_close(result,'small backend-policy parity',Q,P,0,0,...
                  'small unscaled generators use the CPU Taylor path under both flags');
result=test_true(result,'separate backend records',numel(store.keys)==2,...
                 'CPU and GPU policy must not share a cache identity');

% Retrieve the same policy through a pool worker
future=parfeval(pool,@propagator,1,spin_system,L,dt); R=fetchOutputs(future);
result=test_close(result,'worker policy hit',R,Q,0,0,...
                  'workers must retrieve the same policy record as the client');

end


