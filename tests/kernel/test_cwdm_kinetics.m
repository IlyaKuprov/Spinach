% Direct-sum mass-action and spin-transport tests T8-T11. Syntax:
%
%                    result=test_cwdm_kinetics()
%
% Outputs:
%
%    result - analytic mass-action, exchange, closure, and convergence checks
%
% T8 uses arbitrary spin orders, zero concentrations, and spin-free pools.
% T9 compares unit-coordinate exchange with expm of the chemical rate matrix.
% T10 compares the production RKMK4 route with ode45 at 1e-12/1e-14 tolerances.
% T11 resolves the cross-reactant product order with an unweighted coil.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_kinetics()

% Construct A+B to C with a spin-free product sink
fprintf('TESTING: CWDM mass action and reaction maps T8-T11\n');
result=new_test_result('kernel/cwdm_kinetics','CWDM chemistry T8-T11',...
                      'Direct-sum reaction generators against analytic references.');
sys.magnet=1; sys.isotopes={'1H','1H','1H','1H'};
sys.output='hush'; sys.disable={'hygiene'}; sys.parallel={'local',1};
inter.chem.parts={1,2,3:4,[]}; inter.chem.concs=[0.7 0.3 0 0];
reaction=struct('reactants',[1 2],'products',3,'matching',[1 3;2 4],'rate',50);
inter.chem.reactions={reaction};
bas.formalism='sphten-liouv'; bas.approximation={'none','none','none','none'};
s=basis(create(sys,inter),bas); K=kinetics(s,'report');
eta=unit_state(s)+0.1*state(s,'Lz',1)+0.2*state(s,'Lz',2);
product=coil_state(s,'Lz',3,'exact')+2*coil_state(s,'Lz',4,'exact');
reference=-50*0.3*(0.7*coil_state(s,'E',1,'exact')+0.1*state(s,'Lz',1))...
          -50*0.7*(0.3*coil_state(s,'E',2,'exact')+0.2*state(s,'Lz',2))...
          +10.5*(coil_state(s,'E',3,'exact')+0.1*product);
result=test_close(result,'T8 additive full source',K(0,eta)*eta,reference,1e-12,0,...
                  'loss follows the reactant spin state; unit arrival is counted once');
result=test_close(result,'T8 units',chem_concs(s,K(0,eta)*eta),[-10.5 -10.5 10.5 0],1e-12,0,...
                  'mass action gives k*cA*cB for each event');

% Sample states independently of their spin orders and initial populations
rng(8); balance=[1;1;2;2];
for n=1:8
    eta=randn(s.bas.offsets(end),1);
    concs=rand(1,4); if n==1, concs(1)=0; end
    eta(s.bas.offsets(1:end-1)+1)=concs;
    deriv=chem_concs(s,K(0,eta)*eta);
    expected=50*concs(1)*concs(2)*[-1 -1 1 0];
    result=test_close(result,['T8 random ' num2str(n)],deriv,expected,0,1e-12,...
                      'spin orders do not enter nonselective mass action');
    result=test_close(result,['T8 balance ' num2str(n)],deriv*balance,0,1e-12,0,...
                      'the product contains two reactant atom equivalents');
end

% Report mixed row and column memberships without changing the generator
mixed=s; mixed.chem.parts={1,2,[3;4],[]};
mixed.chem.reactions={struct('reactants',[3 1],'products',[],...
                            'matching',zeros(0,2),'rate',2,'closure','additive')};
ordinary=kinetics(mixed); mixed.sys.output=1;
[text,reported]=evalc('kinetics(mixed,''report'');');
result=test_true(result,'mixed membership reporting',contains(text,'traced spins [1 3 4]'),...
                 'reporting accepts row and column spin memberships in the same reaction');
result=test_close(result,'report generator invariance',reported(0,eta),ordinary(0,eta),0,0,...
                  'formatting reaction membership does not change the assembled generator');

% Check a reverse first-order reaction and a tracked spin-free sink
reverse=struct('reactants',3,'products',[1 2],'matching',[3 1;4 2],'rate',2);
sink=struct('reactants',3,'products',4,'matching',zeros(0,2),'rate',3);
inter.chem.reactions={reaction,reverse,sink};
s=basis(create(sys,inter),bas); K=kinetics(s);
eta=unit_state(s); eta(s.bas.offsets(3)+1)=0.4;
result=test_close(result,'T8 reverse and sink',chem_concs(s,K(0,eta)*eta),...
                  [-9.7 -9.7 8.5 1.2],1e-12,0,'reverse and sink stoichiometries add');

% Read and evolve different concentrations in two independent voxels
stack=[eta;2*eta]; concs=chem_concs(s,stack);
result=test_close(result,'voxel unit extraction',concs,[0.7 0.3 0.4 0;1.4 0.6 0.8 0],0,0,...
                  'space-times-spin layout carries an independent concentration row per voxel');
deriv=K(0,stack)*stack;
result=test_close(result,'voxel reaction independence',deriv,...
                  [K(0,eta)*eta;K(0,2*eta)*(2*eta)],1e-12,0,'no inter-voxel reaction entries');

% Repeated spin-free reactants have polynomial rates at zero population
inter.chem.reactions={struct('reactants',[4 4],'products',4,'matching',zeros(0,2),'rate',3)};
s=basis(create(sys,inter),bas); K=kinetics(s);
for conc=[0 0.4]
    eta=unit_state(s); eta(end)=conc;
    result=test_close(result,'repeated spin-free loss',chem_concs(s,K(0,eta)*eta),...
                      [0 0 0 -3*conc^2],1e-12,0,'two drains and one fill give net loss k*c^2');
end

% First-order exchange transports every coordinate including the unit
forward=struct('reactants',1,'products',2,'matching',[1 2],'rate',2);
reverse=struct('reactants',2,'products',1,'matching',[2 1],'rate',3);
inter.chem.reactions={forward,reverse};
s=basis(create(sys,inter),bas); K=kinetics(s); eta=unit_state(s);
expected=expm([-2 3;2 -3]*0.73)*[0.7;0.3];
result=test_true(result,'T9 sparse first order',issparse(K),'constant first order returns a sparse matrix');
result=test_close(result,'T9 exchange expm',chem_concs(s,expm(full(K)*0.73)*eta),...
                  [expected' 0 0],1e-11,0,'unit coordinates follow the two-state chemical rate matrix');
result=test_close(result,'T9 polarisation map',K*state(s,'Lz',1),...
                  -2*state(s,'Lz',1)+1.4*coil_state(s,'Lz',2,'exact'),1e-12,0,...
                  'matched polarisation moves with its source concentration');

% Time-dependent rates are evaluated at the requested stage time
inter.chem.reactions={forward}; inter.chem.reactions{1}.rate=@(t)2+t;
s=basis(create(sys,inter),bas); K=kinetics(s); eta=unit_state(s);
result=test_close(result,'time rate',chem_concs(s,K(0.5,eta)*eta),[-1.75 1.75 0 0],1e-12,0,...
                  'the first-order rate at t=0.5 is 2.5');

% Resolve a time-only callback once for a multi-voxel generator evaluation
rate_calls=0; timed=s; timed.chem.reactions{1}.rate=@counted_rate;
timed_gen=kinetics(timed); stack=[eta;2*eta;3*eta];
deriv=timed_gen(0.5,stack)*stack;
result=test_true(result,'one callback per stage',rate_calls==1,...
                 'a time-only rate is evaluated once regardless of voxel count');
reference=[K(0.5,eta)*eta;K(0.5,2*eta)*(2*eta);K(0.5,3*eta)*(3*eta)];
result=test_close(result,'shared voxel rate',deriv,reference,0,0,...
                  'all voxels use the same resolved rate at the same stage time');

% Share each time schedule across unequal voxel populations
profile clear; profile on;
spatial=K(0.5,[eta;2*eta;3*eta]);
profile off; timing=profile('info');
callback=contains({timing.FunctionTable.FunctionName},'@(t)2+t');
calls=sum([timing.FunctionTable(callback).NumCalls]);
result=test_true(result,'one schedule call',calls==1,...
                 'the shared stage time requires one callback evaluation, not one per voxel');
result=test_close(result,'time rate across voxels',...
                  chem_concs(s,spatial*[eta;2*eta;3*eta]),...
                  [1;2;3]*[-1.75 1.75 0 0],1e-12,0,...
                  'the same schedule rate acts on each local concentration');

% Static zero-rate higher-order records remain usable in linear contexts
inactive=s; inactive.chem.reactions={reaction};
inactive.chem.reactions{1}.rate=0;
inactive.chem.reactions{1}.closure='additive';
zero_gen=kinetics(inactive);
result=test_true(result,'zero higher-order matrix',issparse(zero_gen)&&nnz(zero_gen)==0,...
                 'numeric zero rates give the exact static zero generator');
parameters.spins={'1H'}; parameters.offset=0;
parameters.sweep=100; parameters.npoints=4; parameters.decouple={};
parameters.rho0=state(inactive,'L+','1H');
parameters.coil=coil_state(inactive,'L+','1H','exact');
zero_fid=liquid(inactive,@acquire,parameters,'nmr');
inactive.chem.reactions={};
empty_fid=liquid(inactive,@acquire,parameters,'nmr');
result=test_close(result,'zero higher-order acquisition',zero_fid,empty_fid,1e-12,0,...
                  'a standard linear acquisition agrees with absent reactions');
inactive.chem.reactions={reaction}; inactive.chem.reactions{1}.rate=@(t)0*t;
inactive.chem.reactions{1}.closure='additive';
result=test_true(result,'zero callback stays dynamic',isa(kinetics(inactive),'function_handle'),...
                 'a callback cannot be classified from a single sampled rate');

% Product closure adds cross-reactant order without altering concentrations
inter.chem.reactions={reaction}; inter.chem.reactions{1}.closure='product';
s=basis(create(sys,inter),bas); K=kinetics(s);
eta=unit_state(s)+0.1*state(s,'Lz',1)+0.2*state(s,'Lz',2);
coil=coil_state(s,{'Lz','Lz'},{3,4},'exact');
expected=50*0.7*0.3*0.1*0.2;
result=test_close(result,'T11 product order',(coil'*(K(0,eta)*eta))/(coil'*coil),...
                  expected,1e-10,0,'the matched order receives k*cA*cB*pA*pB');
s.chem.reactions{1}.closure='additive'; K=kinetics(s);
result=test_close(result,'T11 additive order',coil'*(K(0,eta)*eta),0,1e-10,0,...
                  'additive closure omits cross-reactant polarisation products');

% Partial tracing destroys unmatched source orders and creates only identity
traced=struct('reactants',3,'products',1,'matching',[3 1],'rate',2);
partial=s; partial.chem.reactions={traced}; partial.chem.reactions{1}.closure='additive';
trace_gen=kinetics(partial); source=coil_state(s,{'Lz','Lz'},{3,4},'exact');
result=test_close(result,'unmatched source trace',trace_gen*source,-2*source,1e-12,0,...
                  'an order involving a traced spin has no product arrival');

% Report missing source orders when a reactant basis is truncated
restricted=bas; restricted.projections={0,[],[],[]};
small=basis(create(sys,inter),restricted); small.sys.output=1;
[text,maps]=evalc('react_gen(small,small.chem.reactions{1});');
result=test_true(result,'truncated source report',contains(text,'8 product rows without a source')&&...
                 size(maps{1},1)==8,'truncated source descriptors are counted instead of invented');

% Reject ambiguous repeated-spin matching and invalid callback rates
repeated=s.chem.reactions{1}; repeated.reactants=[1 1];
repeated.matching=[1 3]; rejected=false;
try
    react_gen(s,repeated);
catch err
    rejected=strcmp(err.identifier,'Spinach:react_gen:repeatedMatching');
end
result=test_true(result,'repeated matching guard',rejected,'global spin labels cannot identify molecular occurrences');
bad=s; bad.chem.reactions{1}.rate=@(t)-1; bad_gen=kinetics(bad); rejected=false;
try
    bad_gen(0,eta);
catch err
    rejected=strcmp(err.identifier,'Spinach:kinetics:rateValue');
end
result=test_true(result,'time rate guard',rejected,'a callback cannot return a negative rate');
rejected=false;
try
    chem_concs(s,eta(1:end-1));
catch err
    rejected=strcmp(err.identifier,'Spinach:chem_concs:state');
end
result=test_true(result,'voxel length guard',rejected,'each voxel must contain a complete spin block');

% Reject ambiguous matched product copies but retain spin-free stoichiometry
repeated=s.chem.reactions{1}; repeated.products=[3 3]; rejected=false;
try
    react_gen(s,repeated);
catch err
    rejected=strcmp(err.identifier,'Spinach:react_gen:repeatedProductMatching');
end
result=test_true(result,'repeated product matching guard',rejected,...
                 'global destination labels cannot distinguish molecular product occurrences');
repeated.products=[4 4]; repeated.matching=zeros(0,2);
free=s; free.chem.reactions={repeated}; free_gen=kinetics(free); free_eta=unit_state(free);
result=test_close(result,'repeated spin-free products',chem_concs(free,free_gen(0,free_eta)*free_eta),...
                  [-10.5 -10.5 0 21],1e-12,0,'two unlabelled product occurrences carry twice the event population');

% Fourth-order convergence through the shipped state-dependent stepper
eta=unit_state(s); rhs=@(t,y)K(t,y)*y;
[~,trajectory]=ode45(rhs,[0 0.02],full(eta),odeset('RelTol',1e-12,'AbsTol',1e-14));
reference=trajectory(end,:)'; errors=zeros(1,3);
for n=1:3
    dt=0.002/2^(n-1); current=eta;
    for k=1:round(0.02/dt)
        current=step(s,{@(t,y)1i*K(t,y),(k-1)*dt,'RKMK4'},current,dt);
    end
    errors(n)=norm(current-reference);
end
fprintf('CWDM_T10_ERRORS %.12g %.12g %.12g RATIOS %.8g %.8g\n',...
        errors,errors(1:2)./errors(2:3));
result=test_true(result,'T10 fourth order',all(errors(1:2)./errors(2:3)>15)&&...
                 all(errors(1:2)./errors(2:3)<17),'halving the step decreases error by approximately sixteen');
result=test_true(result,'T10 absolute error',errors(3)<1e-10,'the finest-step error is below 1e-10');

% Count production rate evaluations without changing their physical value
function rate=counted_rate(t)
    rate_calls=rate_calls+1; rate=2+t;
end

end


