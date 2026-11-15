% Validate native context dispatch, natural output and capability rejection.
% Syntax:
%
%                       result=test_clexpn_context()
%
% ilya.kuprov@weizmann.ac.il

function result=test_clexpn_context()

% An Ising-coupled pair has an exactly refocused Hahn response
result=new_test_result('kernel/clexpn_context','Native cluster context',...
                       'Registered five-argument dispatch and natural numeric output.');
sys.isotopes={'1H','13C'}; sys.magnet=0; sys.output='hush'; sys.disable={'hygiene'};
inter.zeeman.matrix={zeros(3),zeros(3)}; inter.coupling.matrix=cell(2);
inter.coupling.matrix{1,2}=diag([0 0 7]); spin_system=create(sys,inter);
p.observed=1; p.times=[0 0.01 0.03]; p.transition.spins=1;
p.transition.alpha=[1;0]; p.transition.beta=[0;1];
p.cce.numerics=struct('atol',1e-12,'rtol',1e-10,'max_memory_bytes',1e8);
[answer,diagnostics]=clexpn(spin_system,@hahn_echo,p,'labframe');
result=test_close(result,'Natural echo output',answer,ones(3,1),1e-13,1e-13,...
                  'Ideal midpoint refocusing cancels the static conditional field.');
result=test_true(result,'Natural numeric rather than wrapper',isnumeric(answer)&&...
                  isequal(size(answer),[3 1]),'The context returns the sequence array.');

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
p.cce.approximation='projected'; projected=clexpn(spin_system,@hahn_echo,p,'labframe');
result=test_close(result,'Projected invariant transition',projected,answer,1e-13,1e-13,...
                  'Explicit branch generators reproduce the physical echo.');

% Explicit evaluator identity survives conditional subsets and full-bath paths
[legacy_model,q]=fixture(); ensemble=cce_states(legacy_model,q);
q.interlaced=struct('reference_state',1,'max_states',16,'condition_sites',{{[2 3]}});
family=cce_clusters(legacy_model,cce_graph(legacy_model,q),q);
selected=family(find(arrayfun(@(item)isequal(item.sites,[2 3]),family),1));
calls=0; conditioning=true; parent_assignment=[];
schedule=cce_sequence(legacy_model,q);
raw=cce_interlace(legacy_model,selected,ensemble,schedule,q,@counted);
result=test_true(result,'Recursive registered evaluator',calls>4&&conditioning&&...
                 isfield(raw,'factor')&&isfield(raw,'value'),...
                 'Every denominator keeps the parent conditioning and explicit evaluator.');
q.order=3; q.interlaced.condition_sites={2:4}; calls=0;
family=cce_clusters(legacy_model,cce_graph(legacy_model,q),q);
selected=family(end); raw=cce_interlace(legacy_model,selected,ensemble,schedule,q,@counted);
result=test_true(result,'Full-bath registered evaluator',calls>8&&isfield(raw,'factor'),...
                 'The complete-bath branch does not bypass the supplied evaluator.');

% Capture physical conditioning rather than merely counting numerical outputs
function item=counted(model,cluster,configuration,sequence,parameters)
    calls=calls+1;
    if isempty(parent_assignment), parent_assignment=configuration.assignment; end
    if ~all(ismember(4,q.interlaced.condition_sites{1}))
        conditioning=conditioning&&configuration.assignment(3)==parent_assignment(3);
    end
    item=cce_cluster(model,cluster,configuration,sequence,parameters);
end

end

% Reuse the established physical interlacing fixture without changing its parameters
function [model,p]=fixture()
sys.isotopes={'1H','1H','1H','1H'}; sys.magnet=0; sys.disable={'hygiene'}; sys.parallel={'processes',1}; sys.output='hush';
inter.coupling.matrix=cell(4);
inter.coupling.matrix{1,2}=diag([2 1 5])/(2*pi);
inter.coupling.matrix{1,3}=diag([0 0 4])/(2*pi);
inter.coupling.matrix{1,4}=diag([0 0 7])/(2*pi);
spin_system=create(sys,inter); p.observed_spins=1; p.bath_spins=2:4;
p.method='gcce'; p.expansion='standard'; p.quantity='signal'; p.times=[0 0.2 0.4]; p.order=2;
p.rho0=[1 1;1 1]/2; p.transition.alpha=[1;0]; p.transition.beta=[0;1];
p.assumptions='labframe'; p.orientation=[0 0 0]; p.interaction_cutoff=Inf;
p.graph=struct('mode','complete','score','coupling_norm');
p.bath_state=struct('mode','interlaced','samples',1,'seed',31,'basis','uniform',...
                    'probabilities',{{[0.7;0.3],[0.2;0.8],[0.6;0.4]}});
p.meanfield='static'; p.sequence=[];
p.numerics=struct('atol',1e-9,'rtol',1e-9,'max_memory_bytes',2^26);
model=cce_prepare(spin_system,p);
end


