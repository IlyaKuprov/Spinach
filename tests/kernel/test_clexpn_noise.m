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


