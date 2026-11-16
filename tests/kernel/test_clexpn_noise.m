% Native noise moments, connected ensembles and Fourier convention tests.
% Syntax: result=test_clexpn_noise()
% ilya.kuprov@weizmann.ac.il

function result=test_clexpn_noise()
result=new_test_result('kernel/clexpn_noise','Native fluctuation correlation',...
                       'Unequal means, additive empty reference and quantum ordering.');
sys.magnet=0; sys.isotopes={'E','1H'}; sys.output='hush'; sys.parallel={'local',2};
inter.zeeman.scalar={0,0}; inter.coupling.matrix=cell(2); s=create(sys,inter);
p.observed=1; p.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]); p.times=[0 0.2 0.4];
p.cce.order=1; p.cce.bath_state=struct('mode','sampled','basis','uniform',...
    'probabilities',{{[0.8 0.2]}},'samples',16,'seed',3);
p.cce.numerics=struct('atol',1e-11,'rtol',1e-11,'max_memory_bytes',2^26);
p.noise=struct('ordering','ordered','connected',true,'manifold','bath_only',...
               'terms',struct('support',2,'operator',diag([0.5 -0.5])));
[a,d]=clexpn(s,@noise_correlation,p,'labframe'); means=zeros(1,16);
for n=1:16, means(n)=d.raw{n}(end).mean_zero; end
reference=0.25-mean(means)^2;
result=test_close(result,'Unequal conditional means',a,reference*ones(3,1),1e-12,1e-12,...
                  'Subtract ensemble first moments, not conditional covariances.');
result=test_true(result,'Nonzero ensemble variance',reference>0,'A mixture of eigenstates fluctuates.');
p.cce.averaging='internal_then_expansion'; b=clexpn(s,@noise_correlation,p,'labframe');
result=test_close(result,'Additive averaging order',a,b,1e-12,1e-12,'Connected ensemble reduction is linear.');
p.cce.bath_state=struct('mode','explicit_product','rho',{{diag([0.8 0.2])}});
p.noise.connected=false; p.noise.terms.operator=[0 1;1 0]/2;
s.inter.zeeman.matrix{2}(3,3)=1.3;
[c,dc]=clexpn(s,@noise_correlation,p,'labframe');

% Descriptor and matrix contracts agree before nonlinear ensemble connection
q=p; q.noise=rmfield(q.noise,'terms');
q.noise.descriptors=struct('operators','Lx','spins',2,'coefficient',1);
[described,dd]=clexpn(s,@noise_correlation,q,'labframe');
result=test_close(result,'Noise descriptor physical identity',described,c,0,0,...
                  'Spin labels reproduce the explicit Hermitian fluctuation.');
result=test_close(result,'Noise descriptor raw moments',dd.raw{1}(end).moment,dc.raw{1}(end).moment,0,0,...
                  'Descriptor restriction preserves the unconnected moments.');

% The optional spectrum is computed after the completed noise ensemble
q.spectrum=struct('convention','exp(+i*omega*t)','sidedness','two-sided',...
    'window','hann','frequency_grid',[-2 0 2],'time_domain','positive',...
    'weights','trapezoid','stationary',true);
[~,transformed]=clexpn(s,@noise_correlation,q,'labframe');
expected=cce_spectrum(q.times,described,struct('spectrum',q.spectrum));
result=test_close(result,'Native noise Fourier composition',transformed.spectrum.value,...
                  expected.value,0,0,'The spectrum uses the completed correlation.');

% A complete deterministic connected correlation retains its direct raw bound
single=q; single.noise.connected=true;
[~,complete]=clexpn(s,@noise_correlation,single,'labframe');
result=test_close(result,'Deterministic full-bath noise estimate',...
    complete.reduction.abs_error_bound,complete.raw{1}(end).abs_error_bound,0,0,...
    'The full-bath reference is not charged a redundant cluster recombination.');

% Every descriptor field is checked, including selections outside the cluster
bad={struct('operators','L+','spins',2,'coefficient',1),...
     struct('operators','Lz','spins',1,'coefficient',1),...
     struct('operators','Lz','spins',[2 2],'coefficient',1),...
     struct('operators','Lz','spins',2,'coefficient',1i),...
     struct('operators','Lz','spins',2,'coefficient',NaN),...
     struct('operators','Lz','spins',[],'coefficient',1)};
for k=1:numel(bad)
    invalid=q; invalid.noise.descriptors=bad{k}; rejected=false;
    try
        clexpn(s,@noise_correlation,invalid,'labframe');
    catch err
        rejected=strcmp(err.identifier,'noise_correlation:descriptors');
    end
    result=test_true(result,['Reject noise descriptor ' num2str(k)],rejected,...
                     'Invalid spin labels, physical IDs and coefficients are rejected.');
end
invalid=q; invalid.noise.descriptors.typo=1; rejected=false;
try
    clexpn(s,@noise_correlation,invalid,'labframe');
catch err
    rejected=strcmp(err.identifier,'cce_options:unknown');
end
result=test_true(result,'Reject unknown descriptor field',rejected,...
                 'Descriptor vocabulary is closed.');

reference=(cos(1.3*p.times(:))+0.6i*sin(1.3*p.times(:)))/4;
result=test_close(result,'Ordered quantum reference',c,reference,1e-12,1e-12,...
                  'An independent single-spin trace fixes the complex ordering sign.');
result=test_close(result,'Unconnected moment retained',c,dc.raw{1}(end).moment,1e-12,1e-12,...
                  'The context requests both experiment outputs.');
p.noise.ordering='symmetrized'; sym=clexpn(s,@noise_correlation,p,'labframe');
result=test_close(result,'Symmetrized real correlation',imag(sym),zeros(size(sym)),1e-12,1e-12,...
                  'Hermitian anticommutator expectation is real.');
result=test_true(result,'Ordering metadata',isfield(dc.raw{1},'mean_time')&&...
                 isfield(dc.raw{1},'mean_zero'),'First moments survive scheduling.');

% Observable cross terms close disconnected sparse-graph supports
sys.isotopes={'E','1H','1H'}; inter.zeeman.scalar={0,0,0}; inter.coupling.matrix=cell(3);
three=create(sys,inter); sparse=p; sparse.cce.order=2;
sparse.cce.graph=struct('mode','coupling','threshold',0);
sparse.cce.bath_state=struct('mode','explicit_product','rho',{{diag([1 0]),diag([1 0])}});
sparse.noise.connected=false; sparse.noise.terms(1)=struct('support',2,'operator',diag([0.5 -0.5]));
sparse.noise.terms(2)=struct('support',3,'operator',diag([0.5 -0.5]));
uncentered=clexpn(three,@noise_correlation,sparse,'labframe');
descriptor=sparse; descriptor.noise=rmfield(descriptor.noise,'terms');
descriptor.noise.descriptors=struct('operators','Lz','spins',[2 3],'coefficient',1);
described=clexpn(three,@noise_correlation,descriptor,'labframe');
result=test_close(result,'Descriptor support graph cross terms',described,uncentered,0,0,...
                  'A spin sum retains cross terms and excludes exterior operators.');

result=test_close(result,'Disconnected noise cross terms',uncentered,ones(size(uncentered)),...
                  1e-12,1e-12,'Sparse graph construction retains the observable support union.');

% Fourier transformation remains a separate explicit positive-sign convention
times=(-8:7)'/16; omega=2*pi*(-8:7)'; frequency=2*pi*3;
options.spectrum=struct('convention','exp(+i*omega*t)','sidedness','two-sided',...
    'window','none','frequency_grid',omega,'time_domain','full','weights','rectangle');
spectrum=cce_spectrum(times,exp(-1i*frequency*times),options);
[~,maximum]=max(abs(spectrum.value));
result=test_close(result,'Fourier positive sign',omega(maximum),frequency,0,0,...
                  'A negative time phase peaks at positive angular frequency.');

% The empty bath contributes its explicitly returned additive reference
records=dc.raw; records{1}(1).value(:)=2; records{1}(2).value(:)=5;
local=dc.parameters; reduced=cce_reduce(records,dc.clusters,local);
result=test_close(result,'Empty reference inclusion-exclusion',reduced.value,5*ones(3,1),...
                  1e-12,1e-12,'The reference is neither omitted nor counted twice.');
end


