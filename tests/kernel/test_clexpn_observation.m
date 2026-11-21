% Test: kernel/clexpn_observation | Native CCE observation
% Tests native preparation, receiver phase, coils, and reference normalisation.
% Syntax:
%
%                   result=test_clexpn_observation()
%
% Outputs:
%
%     result - regression result with recorded checks and failures
%
% ilya.kuprov@weizmann.ac.il

function result=test_clexpn_observation()

% Use complex observed preparation and a non-Hermitian transition detector
result=new_test_result('kernel/clexpn_observation','Native CCE observation',...
                       'Preparation, non-Hermitian detection and actual reference.');
sys.isotopes={'1H','13C'};
sys.magnet=0;
sys.output='hush';
sys.disable={'hygiene'};
spin_system=create(sys,struct());
p.observed=1;
p.times=[0 0.013];
p.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
p.cce.order=1;
p.cce.numerics.max_memory_bytes=2^26;
p.cce.numerics.atol=1e-12;
p.cce.numerics.rtol=1e-10;
[p,~]=cce_options(spin_system,p,@sin,{'labframe'});
model=cce_compile(spin_system,p,{'labframe'});
[cluster,~]=cce_subsystem(model,2,struct(),[]);
p.rho0=[0.5,0.2+0.1i;0.2-0.1i,0.5];
p.coil=[0 1;0 0];
p.receiver_phase=0.37;
p.normalized=true;
[rho,coil,reference,meta]=cce_observation(cluster,p);
expected=exp(-1i*p.receiver_phase)*(0.2+0.1i);
result=test_close(result,'Complex actual reference',reference,expected,2e-16,2e-15,...
                  'Reference detection retains the receiver phase and complex coherence.');
result=test_close(result,'Reference applied once',trace(coil'*rho)/meta.normalization_divisor,1,2e-15,0,...
                  'The state and detector remain unnormalised until experiment reduction.');
result=test_close(result,'Initial signal metadata',meta.initial_signal,reference,0,0,...
                  'Metadata reports the actual reference rather than a fixed one-half.');
result=test_close(result,'Physical density trace',trace(rho),1,2e-15,0,...
                  'The bath product is trace-normalised.');

% Multiple complex detectors preserve time-by-detector orientation
p.normalized=false;
p.coil={p.coil,[0 -1i;1i 0]/2};
[rho,coils,reference,meta]=cce_observation(cluster,p);
[H,~,~,~,~]=cce_generators(cluster,p,[]);
signal=cce_evolution(cluster,H,coils,rho,0,2,'observable');
result=test_close(result,'Multiple detector orientation',signal,repmat(reference,3,1),2e-16,2e-15,...
                  'Observable output is time-by-detector, including time zero.');
result=test_close(result,'Unnormalised divisor',meta.normalization_divisor,[1 1],0,0,...
                  'Unnormalised observations do not acquire an implicit reference division.');

% Descriptor preparation and detection use physical IDs and canonical scaling
spec=rmfield(p,{'rho0','coil','receiver_phase'});
spec.preparation=struct('states','Lx','spins',1);
spec.detection=struct('operators','L+','spins',1);
[rho,coil,reference,~]=cce_observation(cluster,spec);
result=test_close(result,'Descriptor preparation',rho,kron([0 1;1 0]/2,eye(2)/2),0,0,...
                  'Deviation-density preparation retains dimensionless Spinach normalization.');
result=test_close(result,'Non-Hermitian descriptor',coil,kron([0 1;0 0],eye(2)),0,0,...
                  'The raising detector is not silently conjugated or made Hermitian.');
result=test_close(result,'Descriptor reference',reference,0.5,0,0,...
                  'The actual initial signal supplies the normalization reference.');
liouv=cluster;
liouv.formalism='zeeman-liouv';
[vector,detector,~,~]=cce_observation(liouv,spec);
result=test_close(result,'Liouville preparation vectorisation',vector,rho(:),0,0,...
                  'Preparation is vectorised column-major.');
result=test_close(result,'Liouville coil vectorisation',detector,coil(:),0,0,...
                  'The detector uses the same column-major convention.');

% Projected blocks agree with full-space transition detection before scaling
projected=p;
projected.coil=projected.coil{1};
projected.cce.approximation='projected';
[block,detector,reference,meta]=cce_observation(cluster,projected);
result=test_close(result,'Projected actual detected reference',detector'*block,reference,2e-16,2e-15,...
                  'Projected propagation carries the actual prepared coherence amplitude.');
result=test_true(result,'Output orientation declared',strcmp(meta.output_orientation,'time-by-detector'),...
                 'Experiments can rely on an explicit observation orientation.');

% Isotope descriptors exclude matching bath spins from the observed-only view
same_sys=sys;
same_sys.isotopes={'E','E'};
same=create(same_sys,struct());
selected=rmfield(p,{'rho0','coil','receiver_phase'});
selected.observed=2;
selected.transition.spins=2;
selected.cce.bath_spins=1;
selected.preparation=struct('states','Lx','spins','E');
[selected,~]=cce_options(same,selected,@free_evolution,'labframe');
same_model=cce_compile(same,selected,'labframe');
[same_cluster,~]=cce_subsystem(same_model,1,struct(),[]);
[rho,~,~]=cce_observation(same_cluster,selected);
result=test_close(result,'Observed-only isotope preparation',rho,...
                  kron([0 1;1 0]/2,eye(2)/2),0,0,...
                  'A matching bath isotope cannot enter the observed preparation.');
selected.detection=struct('operators','L+','spins','E');
[~,coil,reference]=cce_observation(same_cluster,selected);
result=test_close(result,'Observed-only isotope detection',coil,...
                  kron([0 1;0 0],eye(2)),0,0,...
                  'The observation service resolves the observed physical IDs only.');
result=test_close(result,'Isotope descriptor reference',reference,0.5,0,0,...
                  'Isotope and physical-ID selection have the same reference.');
selected=rmfield(selected,'detection');
selected.normalized=true;
actual=clexpn(same,@free_evolution,selected,'labframe');
result=test_close(result,'Public isotope preparation',actual,ones(numel(selected.times),1),0,0,...
                  'Both the empty and retained-bath tasks accept the observed isotope descriptor.');

% Ambiguous normalisation and unsupported projected observables fail explicitly
bad=p;
bad.rho0=eye(2)/2;
bad.coil=[0 1;0 0];
bad.normalized=true;
result=test_true(result,'Vanishing reference rejected',rejects(@()cce_observation(cluster,bad)),...
                 'A zero initial detected signal cannot be a normalization divisor.');
bad=projected;
bad.coil=eye(2);
result=test_true(result,'Projected population detector rejected',rejects(@()cce_observation(cluster,bad)),...
                 'A projected coherence block cannot produce a population observable.');
bad=spec;
bad.preparation='thermal';
result=test_true(result,'Unknown preparation rejected',rejects(@()cce_observation(cluster,bad)),...
                 'No unimplemented preparation is replaced by a guessed state.');

end

% Catch an expected contract failure without changing the public API
function answer=rejects(callback)

% Record whether the expected contract error occurs
answer=false;
try
    callback();
catch
    answer=true;
end

end


