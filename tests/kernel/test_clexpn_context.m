% Test: kernel/clexpn_context | Native cluster context
% Validate native context dispatch, natural output, and capability rejection.
% Syntax:
%
%                       result=test_clexpn_context()
%
% Outputs:
%
%     result - regression result with recorded checks and failures
%
% ilya.kuprov@weizmann.ac.il

function result=test_clexpn_context()

% An Ising-coupled pair has an exactly refocused Hahn response
result=new_test_result('kernel/clexpn_context','Native cluster context',...
                       'Registered five-argument dispatch and natural numeric output.');
sys.isotopes={'1H','13C'};
sys.magnet=0;
sys.output='hush';
sys.disable={'hygiene'};
inter.zeeman.matrix={zeros(3),zeros(3)};
inter.coupling.matrix=cell(2);
inter.coupling.matrix{1,2}=diag([0 0 7]);
spin_system=create(sys,inter);
p.observed=1;
p.times=[0 0.01 0.03];
p.transition.spins=1;
p.transition.alpha=[1;0];
p.transition.beta=[0;1];
p.cce.numerics=struct('atol',1e-12,'rtol',1e-10,'max_memory_bytes',1e8);
[answer,diagnostics]=clexpn(spin_system,@hahn_echo,p,'labframe');
result=test_close(result,'Natural echo output',answer,ones(3,1),1e-13,1e-13,...
                  'Ideal midpoint refocusing cancels the static conditional field.');
result=test_true(result,'Natural numeric rather than wrapper',isnumeric(answer)&&...
                  isequal(size(answer),[3 1]),'The context returns the sequence array.');

% Observable output may use the detector supplied by the selected transition
observable=p;
observable.output='observable';
implicit=clexpn(spin_system,@free_evolution,observable,'labframe');
observable.coil=[0 1;0 0];
explicit=clexpn(spin_system,@free_evolution,observable,'labframe');
result=test_close(result,'Transition-derived observable detector',implicit,explicit,0,0,...
                  'An omitted coil retains one transition-derived detection column.');

% The physical bath radius reaches graph selection and cluster enumeration
radial=spin_system;
radial.inter.coordinates={[0 0 0],[0 0 2]};
radius=p;
radius.cce.bath_radius=1;
[~,outside]=clexpn(radial,@free_evolution,radius,'labframe');
radius.cce.bath_radius=2;
[~,inside]=clexpn(radial,@free_evolution,radius,'labframe');
result=test_true(result,'Native bath-radius truncation',isempty(outside.graph.bath_ids)&&...
                 isequal(outside.graph.omitted_bath_ids,2)&&isempty(outside.clusters)&&...
                 isequal(inside.graph.bath_ids,2)&&isempty(inside.graph.omitted_bath_ids),...
                 'The explicit Angstrom radius includes its boundary and excludes distant sites.');

% Direct standard-signature evaluation must match the full retained-cluster task
model=cce_compile(spin_system,diagnostics.parameters,'labframe');
[cluster,~]=cce_subsystem(model,2,struct(),[]);
[H,R,K,parameters]=cce_generators(cluster,diagnostics.parameters,[]);
reference=hahn_echo(cluster,parameters,H,R,K);
result=test_close(result,'Standard five-argument dispatch',...
                  diagnostics.raw{1}(end).value,reference,0,0,...
                  'The scheduler invokes the registered experiment without a result adapter.');
contract=cce_capabilities(@hahn_echo);
result=test_true(result,'Explicit registered evaluator',isa(contract.evaluator,'function_handle')&&...
                  isequal(contract.pulse_sequence,@hahn_echo),'Registry owns dispatch.');

% Unsupported handles are rejected before model construction
caught=false;
try
    clexpn(spin_system,@sin,p,'labframe');
catch exception
    caught=strcmp(exception.identifier,'cce_capabilities:unregistered');
end
result=test_true(result,'Unregistered sequence rejected',caught,...
                  'No substring or method-name inference is permitted.');

% Projected and generalised closed-system branches agree for this invariant transition
p.cce.approximation='projected';
projected=clexpn(spin_system,@hahn_echo,p,'labframe');
result=test_close(result,'Projected invariant transition',projected,answer,1e-13,1e-13,...
                  'Explicit branch generators reproduce the physical echo.');

% Exterior pumping follows the selected physical ensemble
p.pump=struct('sites',2,'center',0,'halfwidth',0,...
              'units','Hz','fraction',0.25,'allow_empty',false);
[pumped,pump_diagnostics]=clexpn(spin_system,@bath_pumped_hahn,p,'labframe');
result=test_close(result,'Exterior spin omitted in the empty cluster',pumped,...
                  cos(2*pi*7*p.times(:)/4),1e-12,1e-12,...
                  'The empty cluster omits unconditioned exterior rotations.');
result=test_true(result,'Projected pump diagnostics retained',...
                 isfield(pump_diagnostics.raw{1},'eigen_residual'),...
                 'Branch spectra and arithmetic bounds accompany the signal.');

% Both pulse orderings and coincident pulses have the same Ising modulation
for fraction=[0 0.25 0.5 0.75 1]
    varied=p;
    varied.pump.fraction=fraction;
    observed=clexpn(spin_system,@bath_pumped_hahn,varied,'labframe');
    result=test_close(result,sprintf('Pump fraction %g',fraction),observed,...
                      cos(2*pi*7*p.times(:)*min(fraction,1-fraction)),1e-12,1e-12,...
                      'Endpoint and coincident ideal pulses retain physical ordering.');
end

% A sampled no-field ensemble also allows selected exterior spin rotations
sampled=p;
sampled.cce.bath_state=struct('mode','sampled','basis','uniform',...
                              'probabilities','uniform','samples',4,'seed',19);
[sampled_signal,sampled_diagnostics]=clexpn(spin_system,@bath_pumped_hahn,sampled,'labframe');
assignment=[sampled_diagnostics.ensemble.spin_config.assignment];
expected=mean(exp(-1i*2*pi*7*p.times(:)*(1.5-assignment)/2),2);
result=test_close(result,'Sampled exterior pump ensemble',sampled_signal,expected,1e-12,1e-12,...
                  'Conditioning is absent and completed expansions are averaged.');

p.cce.meanfield='static';
failed=false;
try
    clexpn(spin_system,@bath_pumped_hahn,p,'labframe');
catch exception
    failed=strcmp(exception.identifier,'bath_pumped_hahn:exterior');
end
result=test_true(result,'Static exterior pumping rejected',failed,...
                 'A fixed mean field cannot silently represent a rotated exterior spin.');
p.cce.meanfield='none';
for field={'center','halfwidth','fraction','units','allow_empty','sites','species'}
    invalid=p;
    switch field{1}
        case 'center'
            invalid.pump.center=NaN;
        case 'halfwidth'
            invalid.pump.halfwidth=-1;
        case 'fraction'
            invalid.pump.fraction=2;
        case 'units'
            invalid.pump.units='ppm';
        case 'allow_empty'
            invalid.pump.allow_empty=1;
        case 'sites'
            invalid.pump.sites=1;
        case 'species'
            invalid.pump.species={'nonsense'};
    end
    failed=false;
    try
        clexpn(spin_system,@bath_pumped_hahn,invalid,'labframe');
    catch exception
        failed=startsWith(exception.identifier,'bath_pumped_hahn:')||...
               startsWith(exception.identifier,'cce_options:');
    end
    result=test_true(result,['Pump validates ' field{1}],failed,...
                     'Every physical selector field must be checked before use.');
end

% Interlaced conditioning retains the original configuration rows
interlace_sys.isotopes={'1H','1H','1H','1H'};
interlace_sys.magnet=0;
interlace_sys.disable={'hygiene'};
interlace_sys.parallel={'processes',1};
interlace_sys.output='hush';
interlace_inter.coupling.matrix=cell(4);
interlace_inter.coupling.matrix{1,2}=diag([2 1 5])/(2*pi);
interlace_inter.coupling.matrix{1,3}=diag([0 0 4])/(2*pi);
interlace_inter.coupling.matrix{1,4}=diag([0 0 7])/(2*pi);
interlace_system=create(interlace_sys,interlace_inter);
q.observed=1;
q.cce.bath_spins=2:4;
q.cce.approximation='generalized';
q.normalized=false;
q.times=[0 0.2 0.4];
q.cce.order=2;
q.rho0=[1 1;1 1]/2;
q.transition.alpha=[1;0];
q.transition.beta=[0;1];
q.orientation=[0 0 0];
q.cce.interaction_cutoff=Inf;
q.cce.graph=struct('mode','complete','score','coupling_norm');
q.cce.bath_state=struct('mode','interlaced','samples',1,'seed',31,'basis','uniform',...
                        'probabilities',{{[0.7 0.3],[0.2 0.8],[0.6 0.4]}});
q.cce.meanfield='static';
q.transition.spins=q.observed;
q.cce.interlaced=struct('max_states',16);
q.cce.numerics=struct('atol',1e-9,'rtol',1e-9,'max_memory_bytes',2^26);
[~,interlace_details]=clexpn(interlace_system,@free_evolution,q,'labframe');
q=interlace_details.parameters;
model=cce_compile(interlace_system,q,'labframe');
ensemble=interlace_details.ensemble;
contract=interlace_details.contract;
q.cce.interlaced=struct('reference_state',1,'max_states',16,'condition_sites',{{[2 3]}});
selection=struct('method','gcce','expansion','standard','graph',q.cce.graph,'order',q.cce.order,...
                 'numerics',q.cce.numerics);
family=cce_clusters(model,cce_graph(model,selection),selection);
model.clusters=family;
selected=family(find(arrayfun(@(item)isequal(item.sites,[2 3]),family),1));
calls=0;
conditioning=true;
parent_assignment=[];
schedule=struct('pulse_sequence',@free_evolution,'contract',contract);
raw=cce_interlace(model,selected,ensemble,schedule,q,@counted);
result=test_true(result,'Recursive registered evaluator',calls>4&&conditioning&&...
                 strcmp(raw.kind,'irreducible')&&isfield(raw,'value'),...
                 'Every denominator keeps the parent conditioning and explicit evaluator.');
q.cce.order=3;
q.cce.interlaced.condition_sites={2:4};
calls=0;
selection=struct('method','gcce','expansion','standard','graph',q.cce.graph,'order',q.cce.order,...
                 'numerics',q.cce.numerics);
family=cce_clusters(model,cce_graph(model,selection),selection);
model.clusters=family;
selected=family(end);
raw=cce_interlace(model,selected,ensemble,schedule,q,@counted);
result=test_true(result,'Full-bath registered evaluator',...
                 calls==prod(model.multiplicities(model.bath_ids))&&raw.full_bath,...
                 'Every complete-bath configuration invokes the supplied evaluator exactly once.');

% Capture physical conditioning rather than merely counting numerical outputs
function item=counted(model,cluster,configuration,sequence,parameters)

% Preserve the parent assignment while recording evaluator calls
calls=calls+1;
if isempty(parent_assignment)
    parent_assignment=configuration.assignment;
end
if ~all(ismember(4,q.cce.interlaced.condition_sites{1}))
    conditioning=conditioning&&configuration.assignment(3)==parent_assignment(3);
end
evaluate=contract.evaluator;
item=evaluate(model,cluster,configuration,sequence,parameters);

end

end


