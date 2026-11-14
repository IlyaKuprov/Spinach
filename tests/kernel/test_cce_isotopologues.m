% Finite and sampled H/D isotopologue probability tests.
%
%                         result=test_cce_isotopologues()
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_isotopologues()

% Enumerate unequal independent isotope probabilities exactly
result=new_test_result('kernel/cce_isotopologues','CCE isotopologues',...
                       'H/D probabilities, actual Hilbert dimensions and local sampling.');
s=struct('mode','explicit','isotopes',{{'1H','1H'}},'magnetic',0,...
         'coordinates_A',[0 0 0;1 0 0],'species',{{'1H','2H'}},...
         'occupancy',[0.8 0.2;0.3 0.7],...
         'isotopologues',struct('mode','enumerate','max_configurations',4));
[~,~,bath]=cce_bath(s); configs=bath.configurations;
result=test_close(result,'Unequal independent weights',[configs.weight],[0.24 0.56 0.06 0.14],1e-15,1e-15,'Unequal independent weights');
result=test_true(result,'H/D physical multiplicities',isequal(vertcat(configs.multiplicities),[2 2;2 3;3 2;3 3]),'H/D physical multiplicities');

% Compare weighted complex observables with independently enumerated species values
values=arrayfun(@(c)prod(c.multiplicities)+1i*sum(c.multiplicities),configs);
exact=sum([configs.weight].*values);
independent=0.24*(4+4i)+0.56*(6+5i)+0.06*(6+5i)+0.14*(9+6i);
result=test_close(result,'Complex configuration averaging occurs after each result',exact,independent,1e-14,1e-14,'Complex configuration averaging occurs after each result');

% Compare each species-dependent Hamiltonian against an analytic FID product
sys.isotopes={'E','1H','1H'}; sys.magnet=0; sys.output='hush';
sys.disable={'hygiene'}; sys.parallel={'local',1};
inter.zeeman.matrix=repmat({zeros(3)},1,3); ss=create(sys,inter);
p.central_spins=1; p.bath_spins=[2 3]; p.method='cce'; p.expansion='standard';
p.quantity='coherence'; p.order=2; p.central.alpha=[1;0]; p.central.beta=[0;1];
p.assumptions='labframe'; p.orientation=[0 0 0]; p.interaction_cutoff=Inf;
p.graph=struct('mode','complete','score','coupling_norm'); p.meanfield='none';
p.times=[0 0.001 0.003]; p.sequence=[]; p.bath_state.mode='mixed_product';
p.numerics=struct('atol',1e-10,'rtol',1e-10,'max_memory_bytes',2^24,'reduction','direct');
weighted=zeros(3,1); analytic=zeros(3,1);
for k=1:numel(configs)
    ss.comp.isotopes=[{'E'} configs(k).isotopes]; ss.comp.mults=[2 configs(k).multiplicities];
    reference=ones(3,1); ss.inter.coupling.matrix=cell(3);
    for j=1:2
        coupling=2*pi*(100+50*j)/configs(k).multiplicities(j);
        ss.inter.coupling.matrix{1,j+1}=diag([0 0 coupling]);
        multiplicity=configs(k).multiplicities(j); projections=(multiplicity-1)/2:-1:-(multiplicity-1)/2;
        reference=reference.*mean(exp(-1i*coupling*p.times(:)*projections),2);
    end
    actual=cce(ss,p); weighted=weighted+configs(k).weight*actual.value;
    analytic=analytic+configs(k).weight*reference;
end
result=test_close(result,'Weighted species-specific exact FIDs',weighted,analytic,1e-12,1e-12,...
                  'Rebuild physical spin dimensions and average independently propagated configurations.');

% Preserve declared joint probabilities rather than factorizing them
joint=rmfield(s,{'species','occupancy','isotopologues'});
joint.configurations=struct('isotopes',{{'1H','1H'},{'2H','2H'}},'weight',{0.25,0.75});
[~,~,bath]=cce_bath(joint);
result=test_close(result,'Correlated joint isotope weights',[bath.configurations.weight],[0.25 0.75],0,0,'Correlated joint isotope weights');

% Test a fixed sampling confidence bound and global-RNG independence
s.isotopologues=struct('mode','sample','count',6000); s.seed=431;
state=rng; [~,~,sample]=cce_bath(s);
sampled=arrayfun(@(c)prod(c.multiplicities),sample.configurations);
reference=real(exact); variance=sum([configs.weight].*(real(values)-reference).^2);
error=abs(mean(sampled)-reference); bound=5*sqrt(variance/numel(sampled));
result=test_true(result,'Sampling agrees within a declared five-sigma bound',error<bound&&isequal(state,rng),'Sampling agrees within a declared five-sigma bound');
fprintf('CCE isotopologues: mean error %.6g; five-sigma bound %.6g; 6000 draws.\n',error,bound);

end
