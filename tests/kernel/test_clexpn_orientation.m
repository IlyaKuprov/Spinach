% Test: kernel/clexpn_orientation | Native CCE orientation averaging
% Validate native CCE quadrature, averaging order, and input rejection.
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
result=new_test_result('kernel/clexpn_orientation','Native CCE quadrature',...
                       'Independent orientation dynamics and explicit averaging order.');
sys.magnet=0.01;
sys.isotopes={'E','2H','2H'};
sys.output='hush';
sys.parallel={'processes',1};
inter.zeeman.scalar={0,0,0};
inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=diag([1e5 1e5 4e5]);
inter.coupling.matrix{1,3}=diag([2e4 2e4 8e4]);
inter.coupling.matrix{2,2}=diag([-1 -1 2])*2e5;
sys.tols.cce_atol=1e-6;
sys.tols.cce_rtol=1e-6;
if poolsize()>1
    delete(gcp('nocreate'));
end
spin_system=create(sys,inter);
p.observed=1;
p.transition=struct('spins',1,'ms',[0.5 -0.5]);
p.preparation='transverse';
p.times=[0 13e-6 41e-6];
p.cce=struct('bath_spins',2:3,'order',1,'approximation','generalized');
grid=struct('alphas',[0 0],'betas',[0 pi/2],'gammas',[0 0],'weights',[1 3]);
p.grid=grid;
[after,diagnostics]=clexpn(spin_system,@hahn_echo,p,'labframe');
result=test_true(result,'Normalised quadrature',isequal(diagnostics.grid.weights,[0.25;0.75]),...
                 'Relative weights sum to one in the physical average.');

% Pin the independently captured complex amplitudes from the original powder experiment
reference=[1                                    1;...
           0.691207898735802-0.086938954256301i 0.836085214520205-0.006540327054106i;...
           0.830400754296497+0.015773469814684i 0.602738367058514+0.002291040734643i];
result=test_close(result,'Reference orientations',cat(2,diagnostics.values{:}),reference,2e-12,0,...
                  'Both separately rotated echoes reproduce the original complex quadrature components.');
result=test_close(result,'Reference weighted signal',after,reference*[0.25;0.75],2e-12,0,...
                  'The signal has the original amplitude, phase, ordering, and normalisation.');

% Expansion-first averaging equals independently completed orientations
expected=zeros(size(after));
for n=1:2
    local=rmfield(p,'grid');
    local.orientation=[grid.alphas(n) grid.betas(n) grid.gammas(n)];
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

% Unsummed output must retain physical signals, not weighted contributions
resolved=p;
resolved.sum_up=false;
resolved.grid=struct('alphas',[0 0 0],'betas',[pi/2 1 0],'gammas',[0 0 0],'weights',[3 0 1]);
for mode={'after_expansion','internal_then_expansion'}
    resolved.cce.averaging=mode{1};
    [signals,data]=clexpn(spin_system,@hahn_echo,resolved,'labframe');
    result=test_true(result,'Unsummed orientation metadata',isequal(size(signals),[2 1])&&...
                     isequal(signals,data.values)&&isequal(data.grid_indices,[1;3])&&...
                     isequal(data.grid.weights,[0.75;0.25]),...
                     'Only positive nodes remain, in input order, with original indices.');
    result=test_close(result,'Unsummed orientation amplitudes',cat(2,signals{:}),...
                      reference(:,[2 1]),2e-12,0,'Each returned echo has unit quadrature weight.');
end

% A failing orientation invalidates only its affected times and expansion orders
limited=spin_system;
limited.tols.cce_atol=2.55e-8;
limited.tols.cce_rtol=0;
[~,masked]=clexpn(limited,@hahn_echo,p,'labframe');
a=masked.orientations{1}.reduction;
b=masked.orientations{2}.reduction;
result=test_true(result,'Mixed orientation validity',~a.valid(end)&&b.valid(end)&&...
                 isequal(masked.reduction.valid,a.valid&b.valid)&&...
                 isequal(masked.reduction.by_order_valid,a.by_order_valid&b.by_order_valid)&&...
                 isequal(masked.reduction.reason(end),a.by_order_reason(end,:,end)),...
                 'Only the failing node supplies the invalid final point and reason.');
result=test_close(result,'Weighted orientation bounds',diagnostics.reduction.by_order_error,...
                  0.25*diagnostics.orientations{1}.reduction.by_order_error+...
                  0.75*diagnostics.orientations{2}.reduction.by_order_error,0,0,...
                  'Absolute uncertainty is averaged conservatively at each expansion order.');

% Cluster parallelism preserves each complete orientation and uses the existing pool
delete(gcp('nocreate'));
parallel_sys=sys;
parallel_sys.parallel={'processes',2};
parallel_system=create(parallel_sys,inter);
[parallel_answer,parallel_data]=clexpn(parallel_system,@hahn_echo,p,'labframe');
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
future=parfeval(pool,@clexpn,2,parallel_system,@hahn_echo,p,'labframe');
[worker_answer,worker_data]=fetchOutputs(future);
result=test_close(result,'Outer-worker orientation signal',worker_answer,after,0,0,...
                  'An outer worker executes nested cluster calculations serially.');
result=test_true(result,'Nested-worker suppression',...
                 isequal(worker_data.reduction.valid,diagnostics.reduction.valid)&&...
                 gcp('nocreate').NumWorkers==2,...
                 'Only the already existing outer pool owns workers.');

% Rotating an anisotropic Zeeman tensor changes the physical pump selection
selected_system=spin_system;
selected_system.inter.zeeman.matrix{2}=diag([500 500 1000]);
pumped=p;
pumped.pump=struct('sites',2,'center',1000,'halfwidth',100,...
                   'units','rad/s','fraction',0.25,'allow_empty',true);
[~,pump_data]=clexpn(selected_system,@bath_pumped_hahn,pumped,'labframe');
selected=cell(1,2);
for n=1:2
    resolved=pump_data.orientations{n}.parameters;
    model_p=cce_options(selected_system,resolved,@cce_compile,'labframe');
    model=cce_compile(selected_system,model_p,'labframe');
    [local_system,~]=cce_subsystem(model,model_p,model_p.cce.bath_spins,struct(),[]);
    [~,selection]=bath_pumped_hahn(local_system,resolved,zeros(local_system.cce.dimension),0,0);
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
[before,internal]=clexpn(spin_system,@hahn_echo,p,'labframe');
records=diagnostics.orientations{1}.raw{1};
averaged=zeros(numel(p.times),numel(records));
for n=1:2
    records=diagnostics.orientations{n}.raw{1};
    for k=1:numel(records)
        averaged(:,k)=averaged(:,k)+diagnostics.grid.weights(n)*records(k).value;
    end
end
expected=averaged(:,1).*prod(averaged(:,2:end)./averaged(:,1),2);
result=test_close(result,'Internal averaging before expansion',before,expected,2e-14,2e-14,...
                  'Raw cluster responses and the observed-only reference are averaged first.');
reference=[1;0.797052697542094-0.026451925330292i;0.649352634531459-0.002949765953454i];
result=test_close(result,'Reference internal-first signal',before,reference,2e-12,0,...
                  'Averaging the cluster records before expansion reproduces the original complex echo.');
result=test_true(result,'Averaging orders do not commute',...
                 max(abs(before-after))>1e-6&&numel(internal.orientations)==2,...
                 'This anisotropic truncated expansion distinguishes the two estimators.');

% A patch tolerance between the two measured truncation errors selects different families
adaptive=p;
adaptive.cce.expansion='partition';
state=diagnostics.orientations{1}.ensemble.spin_config;
adaptive.cce.partition=struct('mode','adaptive','max_dimension',9,'patch_atol',0.175,...
                              'patch_rtol',0,'max_rounds',2,'max_probes',12,'patches',{{[2 3]}},...
                              'pilot_configs',state,'heldout_configs',state);
adaptive.sum_up=false;
[~,families]=clexpn(spin_system,@hahn_echo,adaptive,'labframe');
result=test_true(result,'Different retained families',...
                 isequal({families.orientations{1}.clusters.sites},{[2 3]})&&...
                 isequal({families.orientations{2}.clusters.sites},{2,3}),...
                 'Only the first crystallite needs a merged pair to meet the patch tolerance.');
adaptive.sum_up=true;
refused=false;
try
    clexpn(spin_system,@hahn_echo,adaptive,'labframe');
catch exception
    refused=strcmp(exception.identifier,'clexpn:clusters');
end
result=test_true(result,'No averaging of incompatible families',refused,...
                 'Raw responses of different physical clusters cannot be averaged together.');

% Common sampled spin states are averaged across orientations before their error estimate
p.cce.averaging='after_expansion';
p.bath_samples=4;
p.bath_seed=19;
[~,sampled]=clexpn(spin_system,@hahn_echo,p,'labframe');
state_values=zeros(numel(p.times),4);
for n=1:2
    for k=1:4
        state_values(:,k)=state_values(:,k)+diagnostics.grid.weights(n)*...
                         sampled.orientations{n}.reduction.samples{k}.value;
    end
end
result=test_close(result,'Powder sampled-state error',sampled.reduction.sampling_standard_error,...
                  std(state_values,0,2)/2,0,0,...
                  'A common spin draw across orientations is one statistical sample.');
p.cce.averaging='internal_then_expansion';
[internal_sampled,~]=clexpn(spin_system,@hahn_echo,p,'labframe');
averaged=zeros(numel(p.times),3);
for n=1:2
    for k=1:4
        records=sampled.orientations{n}.raw{k};
        for j=1:3
            averaged(:,j)=averaged(:,j)+diagnostics.grid.weights(n)*records(j).weight*records(j).value;
        end
    end
end
expected=averaged(:,1).*prod(averaged(:,2:end)./averaged(:,1),2);
result=test_close(result,'Internal spin and powder average',internal_sampled,expected,2e-14,2e-14,...
                  'Both quadrature and spin-state weights precede inclusion-exclusion.');

% Integer quadratures are promoted before normalisation, not rounded as integers
integer_request=p;
integer_request.grid.weights=uint8([1 3]);
[integer_answer,~]=clexpn(spin_system,@hahn_echo,integer_request,'labframe');
result=test_close(result,'Integer quadrature weights',integer_answer,internal_sampled,0,0,...
                  'Numeric quadrature inputs are normalised in double precision.');

% Named grids use the same kernel/grids lookup and field names as powder
named=p;
named.times=0;
named=rmfield(named,{'bath_samples','bath_seed'});
named.grid='icos_2ang_12pts.mat';
[named_answer,named_diagnostics]=clexpn(spin_system,@hahn_echo,named,'labframe');
result=test_close(result,'Named powder quadrature',named_answer,1,2e-14,0,...
                  'A shipped powder grid resolves through the native context path.');
result=test_true(result,'Named powder grid size',numel(named_diagnostics.grid.weights)==12,...
                 'The installed quadrature is not replaced by a hidden grid.');

% A zero-weight orientation does not run a second physical calculation
p.grid.weights=[1 0];
[single,one]=clexpn(spin_system,@hahn_echo,p,'labframe');
result=test_close(result,'Zero quadrature weight',single,one.orientations{1}.reduction.value,0,0,...
                  'Zero-measure orientations are omitted without changing the response.');

% Single stochastic draws remain unavailable through both averaging orders
request=p;
request.cce.order=1;
interlaced=struct('reference_state',1,'max_states',9);
request.bath_samples=1;
probe_grid=grid;
probe_grid.betas=[0 0];
for mode={'sampled','interlaced'}
    if strcmp(mode{1},'interlaced')
        request.cce.interlaced=interlaced;
    end
    for averaging={'after_expansion','internal_then_expansion'}
        if strcmp(mode{1},'interlaced')&&strcmp(averaging{1},'internal_then_expansion')
            continue;
        end
        request.cce.averaging=averaging{1};
        request.grid=probe_grid;
        saved_rng=rng;
        saved_request=request;
        [~,single_draw]=clexpn(spin_system,@hahn_echo,request,'labframe');
        result=test_true(result,'Unavailable single-draw uncertainty',...
                         all(isnan(single_draw.reduction.sampling_standard_error),'all'),...
                         'Quadrature nodes do not turn one spin draw into independent samples.');
        result=test_true(result,'Orientation sampling state isolation',...
                         isequal(saved_rng,rng)&&isequal(saved_request,request)&&...
                         isequal(single_draw.orientations{1}.ensemble,single_draw.orientations{2}.ensemble),...
                         'The request and global RNG are unchanged; orientations reuse the same spin draws.');
    end
end

% A complete bath still counts independent spin draws rather than quadrature nodes
request.cce=rmfield(request.cce,'interlaced');
request.cce.order=2;
request.cce.averaging='internal_then_expansion';
for count=[1 4]
    request.bath_samples=count;
    request.grid=probe_grid;
    [~,complete]=clexpn(spin_system,@hahn_echo,request,'labframe');
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
request.cce.averaging='after_expansion';
request.cce.interlaced=interlaced;
request.bath_samples=1;
[~,complete]=clexpn(spin_system,@hahn_echo,request,'labframe');
result=test_true(result,'Complete interlaced orientation diagnostics',...
                 all(complete.reduction.sampling_standard_error==0,'all')&&...
                 isequal(size(complete.reduction.samples{1}.by_order_reason),...
                         size(complete.reduction.samples{1}.by_order_valid)),...
                 'Exact full-bath interlacing has zero sampling error and aligned diagnostic arrays.');

% A shared ledger is rejected before any nested calculation writes a checkpoint
request=rmfield(request,{'bath_samples','bath_seed'});
request.cce=rmfield(request.cce,'interlaced');
request.cce.order=1;
checkpoint=[tempname '.mat'];
request.cce.checkpoint=struct('file',checkpoint,'resume',false);
for resume=[false true]
    request.cce.checkpoint.resume=resume;
    request.grid=grid;
    refused=false;
    try
        clexpn(spin_system,@hahn_echo,request,'labframe');
    catch exception
        refused=strcmp(exception.identifier,'clexpn:checkpoint');
    end
    result=test_true(result,'Multi-orientation checkpoint preflight',refused&&~isfile(checkpoint),...
                     'Both write and resume requests fail before creating the first ledger.');
end
request.cce.checkpoint.resume=false;
request.grid.weights=[1 0];
[written,~]=clexpn(spin_system,@hahn_echo,request,'labframe');
cleanup=onCleanup(@()delete(checkpoint));
request.cce.checkpoint.resume=true;
[resumed,~]=clexpn(spin_system,@hahn_echo,request,'labframe');
result=test_close(result,'Single-orientation checkpoint resume',resumed,written,0,0,...
                  'Zero-weight nodes do not prevent a single physical ledger from being resumed.');

% The ordinary five-argument pulse-sequence branch retains its grid contract
bas.formalism='zeeman-hilb';
bas.approximation={'none'};
full_system=basis(spin_system,bas);
ordinary_parameters=struct('grid',grid);
ordinary_parameters.spins={'E'};
[ordinary,ordinary_grid]=powder(full_system,@(~,~,~,~,~)1,...
                                ordinary_parameters,'labframe');
result=test_close(result,'Ordinary powder pulse sequence',ordinary,sum(grid.weights),2e-15,0,...
                  'Ordinary powder integration retains the supplied, unnormalised weights.');
result=test_true(result,'Ordinary powder outputs',isequal(ordinary_grid.weights,[1 3]),...
                 'The existing grid output retains its original shape and values.');

% Test malformed quadratures, mutually exclusive orientation inputs, and output selection
invalid=repmat({p},1,8);
invalid{1}.grid.weights=[-1 2];
invalid{2}.grid.betas=[0 NaN];
invalid{3}.grid.gammas=0;
invalid{4}.grid.weights=[0 0];
invalid{5}.grid=rmfield(grid,'alphas');
invalid{6}.orientation=[0 0 0];
invalid{7}.sum_up=[true false];
invalid{8}.grid='nonexistent_cce_test_grid.mat';
for n=1:numel(invalid)
    rejected=false;
    try
        clexpn(spin_system,@hahn_echo,invalid{n},'labframe');
    catch
        rejected=true;
    end
    result=test_true(result,['Invalid quadrature input ' num2str(n)],rejected,...
                     'Malformed or conflicting orientation requests cannot run a calculation.');
end

end


