% Validate deferred powder composition, averaging order, and input rejection.
% Syntax:
%
%                    result=test_clexpn_orientation()
%
% Outputs:
%
%     result - regression result with recorded checks and failures
%
% ilya.kuprov@weizmann.ac.il

function result=test_clexpn_orientation()

% A small anisotropic system needs no global basis for either averaging order
result=new_test_result('kernel/clexpn_orientation','Deferred CCE powder',...
                       'Independent orientation dynamics and explicit averaging order.');
sys.magnet=0.01; sys.isotopes={'E','2H','2H'};
sys.output='hush'; sys.parallel={'local',1};
inter.zeeman.scalar={0,0,0}; inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=diag([1e5 1e5 4e5]);
inter.coupling.matrix{1,3}=diag([2e4 2e4 8e4]);
inter.coupling.matrix{2,2}=diag([-1 -1 2])*2e5;
spin_system=create(sys,inter);
p.observed=1; p.transition=struct('spins',1,'ms',[0.5 -0.5]);
p.preparation='transverse'; p.times=[0 13e-6 41e-6];
p.cce=struct('bath_spins',2:3,'order',1,'approximation','generalized');
p.cce.graph.mode='complete'; p.cce.numerics=struct('atol',1e-6,'rtol',1e-6);
grid=struct('alphas',[0 0],'betas',[0 pi/2],'gammas',[0 0],'weights',[1 3]);
outer=struct('grid',grid); task=context_task(@clexpn,@hahn_echo,p,'labframe');
[after,diagnostics]=orientation_average(spin_system,task,outer);
[via_powder,sph_grid,powder_diagnostics]=powder(spin_system,task,outer,'labframe');
result=test_close(result,'Deferred powder dispatch',via_powder,after,0,0,...
                  'Powder dispatches before basis or full-generator construction.');
result=test_true(result,'Powder quadrature contract',...
                 isequal(sph_grid.weights,[0.25;0.75])&&...
                 isequal(powder_diagnostics.grid,sph_grid),...
                 'The conventional grid output and third diagnostics output agree.');

% Expansion-first averaging equals independently completed orientations
expected=zeros(size(after));
for n=1:2
    local=p; local.orientation=[grid.alphas(n) grid.betas(n) grid.gammas(n)];
    [current,evidence]=clexpn(spin_system,@hahn_echo,local,'labframe');
    result=test_close(result,'Orientation-resolved signal',diagnostics.values{n},current,0,0,...
                      'Each orientation retains the independently calculated complex amplitude.');
    result=test_close(result,'Orientation-resolved error',...
                      diagnostics.orientations{n}.reduction.abs_error_bound,...
                      evidence.reduction.abs_error_bound,0,0,...
                      'Numerical bounds are retained before quadrature reduction.');
    result=test_true(result,'Orientation-resolved validity',...
                     isequal(diagnostics.orientations{n}.reduction.valid,evidence.reduction.valid),...
                     'Validity masks are retained before quadrature reduction.');
    expected=expected+grid.weights(n)/sum(grid.weights)*current;
end
result=test_close(result,'Expansion before powder average',after,expected,0,0,...
                  'Every orientation has its own dynamics and cluster expansion.');

% Cluster parallelism preserves each complete orientation and uses the existing pool
parallel_task=task; parallel_task.parameters.cce.parallel.mode='clusters';
[parallel_answer,parallel_data]=orientation_average(spin_system,parallel_task,outer);
result=test_close(result,'Parallel orientation signal',parallel_answer,after,0,0,...
                  'Serial and cluster-parallel composition have identical gather order.');
result=test_close(result,'Parallel orientation error',parallel_data.reduction.abs_error_bound,...
                  diagnostics.reduction.abs_error_bound,0,0,...
                  'Parallel execution preserves the numerical error ledger.');
result=test_true(result,'Parallel orientation masks',...
                 isequal(parallel_data.reduction.valid,diagnostics.reduction.valid)&&...
                 isequal(parallel_data.values,diagnostics.values),...
                 'Per-orientation amplitudes and final masks are bitwise identical.');

% Calling the composition on a worker must not request another pool
pool=gcp('nocreate');
future=parfeval(pool,@orientation_average,2,spin_system,parallel_task,outer);
[worker_answer,worker_data]=fetchOutputs(future);
result=test_close(result,'Outer-worker orientation signal',worker_answer,after,0,0,...
                  'An outer worker executes nested cluster calculations serially.');
result=test_true(result,'Nested-worker suppression',...
                 strcmp(worker_data.parameters.cce.parallel.mode,'serial')&&...
                 gcp('nocreate').NumWorkers==pool.NumWorkers,...
                 'Only the already existing outer pool owns workers.');

% Rotating an anisotropic Zeeman tensor changes the physical pump selection
selected_system=spin_system; selected_system.inter.zeeman.matrix{2}=diag([500 500 1000]);
pumped=p; pumped.pump=struct('sites',2,'center',1000,'halfwidth',100,...
                            'units','rad/s','fraction',0.25,'allow_empty',true);
pump_task=context_task(@clexpn,@bath_pumped_hahn,pumped,'labframe');
[~,pump_data]=orientation_average(selected_system,pump_task,outer);
selected=cell(1,2);
for n=1:2
    resolved=pump_data.orientations{n}.parameters;
    model=cce_compile(selected_system,resolved,'labframe');
    [local_system,~]=cce_subsystem(model,model.bath_ids,struct(),[]);
    [~,selection]=bath_pumped_hahn(local_system,resolved,zeros(local_system.dimension),0,0);
    selected{n}=selection.selected_sites;
    expected=clexpn(selected_system,@bath_pumped_hahn,resolved,'labframe');
    result=test_close(result,'Orientation-selected pump response',pump_data.values{n},expected,0,0,...
                      'The composed response equals the separately rotated physical experiment.');
end
result=test_true(result,'Orientation-dependent pump selection',...
                 isequal(selected{1},2)&&isempty(selected{2}),...
                 'Selection is recomputed after rotation, not shared across orientations.');

% Internal averaging includes the empty-cluster normalisation before products
p.cce.averaging='internal_then_expansion';
task=context_task(@clexpn,@hahn_echo,p,'labframe');
[before,internal]=orientation_average(spin_system,task,outer);
records=diagnostics.orientations{1}.raw{1}; averaged=zeros(numel(p.times),numel(records));
for n=1:2
    records=diagnostics.orientations{n}.raw{1};
    for k=1:numel(records)
        averaged(:,k)=averaged(:,k)+sph_grid.weights(n)*records(k).value;
    end
end
expected=averaged(:,1).*prod(averaged(:,2:end)./averaged(:,1),2);
result=test_close(result,'Internal averaging before expansion',before,expected,2e-14,2e-14,...
                  'Raw cluster responses and the observed-only reference are averaged first.');
result=test_true(result,'Averaging orders do not commute',...
                 max(abs(before-after))>1e-6&&numel(internal.orientations)==2,...
                 'This anisotropic truncated expansion distinguishes the two estimators.');

% Common sampled spin states are averaged across orientations before their error estimate
p.cce.averaging='after_expansion';
p.cce.bath_state=struct('mode','sampled','samples',4,'seed',19,...
                        'basis','uniform','probabilities','uniform');
task=context_task(@clexpn,@hahn_echo,p,'labframe');
[~,sampled]=orientation_average(spin_system,task,outer);
state_values=zeros(numel(p.times),4);
for n=1:2
    for k=1:4
        state_values(:,k)=state_values(:,k)+sph_grid.weights(n)*...
                         sampled.orientations{n}.reduction.samples{k}.value;
    end
end
result=test_close(result,'Powder sampled-state error',sampled.reduction.sampling_standard_error,...
                  std(state_values,0,2)/2,0,0,...
                  'A common spin draw across orientations is one statistical sample.');
p.cce.averaging='internal_then_expansion';
task=context_task(@clexpn,@hahn_echo,p,'labframe');
[internal_sampled,~]=orientation_average(spin_system,task,outer);
averaged=zeros(numel(p.times),3);
for n=1:2
    for k=1:4
        records=sampled.orientations{n}.raw{k};
        for j=1:3
            averaged(:,j)=averaged(:,j)+sph_grid.weights(n)*records(j).weight*records(j).value;
        end
    end
end
expected=averaged(:,1).*prod(averaged(:,2:end)./averaged(:,1),2);
result=test_close(result,'Internal spin and powder average',internal_sampled,expected,2e-14,2e-14,...
                  'Both quadrature and spin-state weights precede inclusion-exclusion.');

% Integer quadratures are promoted before normalisation, not rounded as integers
integer_grid=grid; integer_grid.weights=uint8([1 3]);
[integer_answer,~]=orientation_average(spin_system,task,struct('grid',integer_grid));
result=test_close(result,'Integer quadrature weights',integer_answer,internal_sampled,0,0,...
                  'Numeric quadrature inputs are normalized in double precision.');

% Named grids use the same kernel/grids lookup and field names as powder
named=task; named.parameters.times=0;
named.parameters.cce=rmfield(named.parameters.cce,'bath_state');
[named_answer,named_diagnostics]=orientation_average(spin_system,named,...
                                  struct('grid','icos_2ang_12pts.mat'));
result=test_close(result,'Named powder quadrature',named_answer,1,2e-14,0,...
                  'A shipped powder grid resolves through the native context path.');
result=test_true(result,'Named powder grid size',numel(named_diagnostics.grid.weights)==12,...
                 'The installed quadrature is not replaced by a hidden grid.');

% A zero-weight orientation does not run a second physical calculation
outer.grid.weights=[1 0]; [single,one]=orientation_average(spin_system,task,outer);
result=test_close(result,'Zero quadrature weight',single,one.orientations{1}.reduction.value,0,0,...
                  'Zero-measure orientations are omitted without changing the response.');

% Single stochastic draws remain unavailable through both averaging orders
request=task.parameters; request.cce.order=1;
request.cce.interlaced=struct('reference_state',1,'max_states',9);
probe_grid=grid; probe_grid.betas=[0 0];
for mode={'sampled','interlaced'}
    request.cce.bath_state=struct('mode',mode{1},'samples',1,'seed',19,...
                                  'basis','uniform','probabilities','uniform');
    for averaging={'after_expansion','internal_then_expansion'}
        if strcmp(mode{1},'interlaced')&&strcmp(averaging{1},'internal_then_expansion'), continue; end
        request.cce.averaging=averaging{1};
        local_task=context_task(@clexpn,@hahn_echo,request,'labframe');
        saved_rng=rng; saved_task=local_task;
        [~,single_draw]=orientation_average(spin_system,local_task,struct('grid',probe_grid));
        result=test_true(result,'Unavailable single-draw uncertainty',...
                         all(isnan(single_draw.reduction.sampling_standard_error),'all'),...
                         'Quadrature nodes do not turn one spin draw into independent samples.');
        result=test_true(result,'Orientation sampling state isolation',...
                         isequal(saved_rng,rng)&&isequal(saved_task,local_task)&&...
                         isequal(single_draw.orientations{1}.ensemble,single_draw.orientations{2}.ensemble),...
                         'The task and global RNG are unchanged; orientations reuse the same spin draws.');
    end
end

% A complete bath still counts independent spin draws rather than quadrature nodes
request.cce.order=2; request.cce.averaging='internal_then_expansion';
for count=[1 4]
    request.cce.bath_state=struct('mode','sampled','samples',count,'seed',19,...
                                  'basis','uniform','probabilities','uniform');
    local_task=context_task(@clexpn,@hahn_echo,request,'labframe');
    [~,complete]=orientation_average(spin_system,local_task,struct('grid',probe_grid));
    if count==1
        result=test_true(result,'Complete-bath single-draw uncertainty',...
                         all(isnan(complete.reduction.sampling_standard_error),'all'),...
                         'Repeating an orientation cannot make an unavailable variance zero.');
    else
        values=zeros(numel(request.times),count);
        for n=1:2
            for k=1:count
                values(:,k)=values(:,k)+complete.grid.weights(n)*...
                            complete.orientations{n}.raw{k}(end).value;
            end
        end
        result=test_close(result,'Internal complete-bath sampling uncertainty',...
                          complete.reduction.sampling_standard_error,std(values,0,2)/sqrt(count),0,0,...
                          'Every spin draw is averaged over orientations before estimating uncertainty.');
    end
    result=test_true(result,'Orientation value diagnostics remain cells',iscell(complete.values),...
                     'Uncertainty work does not replace the orientation-resolved signals.');
end

% Complete interlacing preserves per-order reason shapes during diagnostic merging
request.cce.averaging='after_expansion'; request.cce.bath_state.mode='interlaced';
request.cce.bath_state.samples=1;
local_task=context_task(@clexpn,@hahn_echo,request,'labframe');
[~,complete]=orientation_average(spin_system,local_task,struct('grid',probe_grid));
result=test_true(result,'Complete interlaced orientation diagnostics',...
                 all(complete.reduction.sampling_standard_error==0,'all')&&...
                 isequal(size(complete.reduction.samples{1}.by_order_reason),...
                         size(complete.reduction.samples{1}.by_order_valid)),...
                 'Exact full-bath interlacing has zero sampling error and aligned diagnostic arrays.');

% A shared ledger is rejected before any nested calculation writes a checkpoint
request.cce=rmfield(request.cce,'bath_state'); request.cce.order=1;
checkpoint=[tempname '.mat']; request.cce.checkpoint=struct('file',checkpoint,'resume',false);
for resume=[false true]
    request.cce.checkpoint.resume=resume;
    local_task=context_task(@clexpn,@hahn_echo,request,'labframe'); refused=false;
    try
        orientation_average(spin_system,local_task,struct('grid',grid));
    catch exception
        refused=strcmp(exception.identifier,'orientation_average:checkpoint');
    end
    result=test_true(result,'Multi-orientation checkpoint preflight',refused&&~isfile(checkpoint),...
                     'Both write and resume requests fail before creating the first ledger.');
end
request.cce.checkpoint.resume=false;
local_task=context_task(@clexpn,@hahn_echo,request,'labframe');
one_grid=grid; one_grid.weights=[1 0];
[written,~]=orientation_average(spin_system,local_task,struct('grid',one_grid));
cleanup=onCleanup(@()delete(checkpoint));
local_task.parameters.cce.checkpoint.resume=true;
[resumed,~]=orientation_average(spin_system,local_task,struct('grid',one_grid));
result=test_close(result,'Single-orientation checkpoint resume',resumed,written,0,0,...
                  'Zero-weight nodes do not prevent a single physical ledger from being resumed.');

% The ordinary five-argument pulse-sequence branch retains its grid contract
bas.formalism='zeeman-hilb'; bas.approximation={'none'};
full_system=basis(spin_system,bas);
ordinary_parameters=struct('grid',grid); ordinary_parameters.spins={'E'};
[ordinary,ordinary_grid,ordinary_diagnostics]=powder(full_system,@(~,~,~,~,~)1,...
                                                  ordinary_parameters,'labframe');
result=test_close(result,'Ordinary powder pulse sequence',ordinary,sum(grid.weights),2e-15,0,...
                  'Deferred dispatch does not change the ordinary sequence branch.');
result=test_true(result,'Ordinary powder outputs',isempty(ordinary_diagnostics)&&...
                 isequal(ordinary_grid.weights,[1 3]),...
                 'The existing grid output retains its original shape and values.');

% Every new public field is validated before numerical evaluation
invalid={rmfield(task,'context'),task,task,task,task,task,task,task};
inputs=repmat({struct('grid',grid)},size(invalid));
invalid{2}.context=@powder; invalid{3}.pulse_sequence=1;
invalid{4}.parameters=[]; invalid{5}.assumptions=1;
inputs{6}.grid.weights=[-1 2]; inputs{7}.grid.betas=[0 NaN];
inputs{8}.grid.gammas=0;
for n=1:numel(invalid)
    rejected=false;
    try
        orientation_average(spin_system,invalid{n},inputs{n});
    catch
        rejected=true;
    end
    result=test_true(result,['Invalid deferred input ' num2str(n)],rejected,...
                     'Malformed task fields and quadratures are rejected.');
end
rejected=false;
try
    powder(spin_system,task,struct('grid',grid),'nmr');
catch
    rejected=true;
end
result=test_true(result,'Assumptions agreement',rejected,...
                 'The outer context cannot silently change the nested assumptions.');

end


