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
sys.tols.cce_atol=1e-12;
sys.tols.cce_rtol=1e-10;
sys.tols.cce_mem=1e8;
spin_system=create(sys,inter);
p.observed=1;
p.times=[0 0.01 0.03];
p.transition.spins=1;
p.transition.alpha=[1;0];
p.transition.beta=[0;1];
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
model_p=cce_options(spin_system,diagnostics.parameters,@cce_compile,'labframe');
model=cce_compile(spin_system,model_p,'labframe');
[cluster,~]=cce_subsystem(model,model_p,2,struct(),[]);
[H,R,K,cluster]=cce_generators(cluster,diagnostics.parameters,[]);
reference=hahn_echo(cluster,diagnostics.parameters,H,R,K);
result=test_close(result,'Standard five-argument dispatch',...
                  diagnostics.raw{1}(end).value,reference,0,0,...
                  'The scheduler invokes the registered experiment without a result adapter.');

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
sampled.bath_samples=4;
sampled.bath_seed=19;
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
interlace_sys.tols.cce_atol=1e-9;
interlace_sys.tols.cce_rtol=1e-9;
interlace_sys.tols.cce_mem=2^26;
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
q.bath_samples=1;
q.bath_seed=31;
q.bath_probs={[0.7 0.3],[0.2 0.8],[0.6 0.4]};
q.cce.meanfield='static';
q.transition.spins=q.observed;
q.cce.interlaced=struct('max_states',16);
[~,interlace_details]=clexpn(interlace_system,@free_evolution,q,'labframe');
q=interlace_details.parameters;
model_p=cce_options(interlace_system,q,@cce_compile,'labframe');
model=cce_compile(interlace_system,model_p,'labframe');
ensemble=interlace_details.ensemble;
q.cce.interlaced=struct('reference_state',1,'max_states',16,'condition_sites',{{[2 3]}});
model_p.cce.order=q.cce.order;
family=cce_clusters(model,cce_graph(model,model_p),model_p);
model_clusters=family;
selected=family(find(arrayfun(@(item)isequal(item.sites,[2 3]),family),1));
raw=cce_interlace(model,selected,model_clusters,ensemble,@free_evolution,q);
result=test_true(result,'Recursive conditioned factor',...
                 strcmp(raw.kind,'irreducible')&&isfield(raw,'value'),...
                 'A conditioned cluster returns its irreducible factor.');
q.cce.order=3;
q.cce.interlaced.condition_sites={2:4};
model_p.cce.order=q.cce.order;
family=cce_clusters(model,cce_graph(model,model_p),model_p);
model_clusters=family;
selected=family(end);
raw=cce_interlace(model,selected,model_clusters,ensemble,@free_evolution,q);
result=test_true(result,'Full-bath conditioned factor',...
                 raw.assignment_count==prod(model.comp.mults(model_p.cce.bath_spins))&&raw.full_bath,...
                 'Every complete-bath configuration is evaluated exactly once.');

end


