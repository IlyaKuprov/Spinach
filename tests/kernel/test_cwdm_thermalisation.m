% Tests block-local IME recovery and the shared-identity artefact (T5-T7).
% Syntax:
%
%                    result=test_cwdm_thermalisation()
%
% Outputs:
%
%     result - trace, thermal recovery, and unit-column checks
%
% Two proton pairs recover for 20 s at 14.1 T and 298 K. T5 and T6 use
% an absolute infinity-norm tolerance of 1e-10. T7 explicitly constructs
% the old shared-identity IME representation from the local blocks: the
% Hamiltonian and nonunit relaxation columns intertwine exactly, whereas
% the common source drives both substances at unit concentration. The
% measured stock-kernel comparison is a separate external validation.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_thermalisation()

% Announce the thermalisation invariants
fprintf('TESTING: CWDM thermalisation (T5-T7)\n');
result=new_test_result('kernel/cwdm_thermalisation','CWDM thermalisation',...
                      'IME preserves local populations and has only local unit-column sources.');

% Build the prototype without chemical exchange to isolate the IME source
sys.magnet=14.1; sys.isotopes={'1H','1H','1H','1H'};
inter.chem.parts={1:2,3:4,[]}; inter.chem.concs=[0.7 0.3 0];
inter.zeeman.scalar={1,2,3,4}; inter.coupling.scalar=cell(4);
inter.coupling.scalar{1,2}=10; inter.coupling.scalar{3,4}=8;
inter.relaxation={'t1_t2'}; inter.r1_rates={2,2,2,2}; inter.r2_rates={5,5,5,5};
inter.equilibrium='IME'; inter.temperature=298; inter.rlx_keep='secular';
bas.formalism='sphten-liouv'; bas.approximation={'none','none','none'};
s=assume(test_spin_system(sys,inter,bas),'nmr');
H=hamiltonian(s); R=relaxation(s); target=equilibrium(s);
units=s.bas.offsets(1:end-1)+1; nonunit=setdiff(1:s.bas.offsets(end),units);

% Check every sampled trace and the recovered thermal state after 20 seconds
rho=full(unit_state(s)); unit_error=0;
for t=[0 0.1 1 5 20]
    evolved=expm(full((-1i*H+R)*t))*rho;
    unit_error=max(unit_error,norm(evolved(units)-s.chem.concs(:),inf));
end
thermal_error=norm(evolved-target,inf);
result=test_true(result,'T5 conserved populations',unit_error<1e-10,...
                 'Hamiltonian and IME evolution preserve every concentration');
result=test_true(result,'T6 recovered equilibrium',thermal_error<1e-10,...
                 'twenty seconds recovers the weighted thermal target');
fprintf('CWDM_T5_T6 unit_error=%.16g thermal_error=%.16g\n',unit_error,thermal_error);

% Check concentration independence, empty blocks, and the public target contract
empty=s; empty.chem.concs=[0 0 1];
result=test_true(result,'zero-population generator',isequal(relaxation(empty),R),...
                 'IME contains unit-concentration sources even for absent substances');
result=test_true(result,'spin-free recovery',...
                 isequal(step(empty,1i*R,full(unit_state(empty)),20),full(unit_state(empty))),...
                 'a populated spin-free block cannot acquire spin order');
zero=s; zero.rlx.equilibrium='zero'; R0=relaxation(zero);
one=s; one.chem.concs(:)=1;
explicit=thermalize(s,R0,[],[],equilibrium(one),'IME');
explicit=complex(clean_up(s,explicit,s.tols.rlx_zero));
result=test_close(result,'explicit IME targets',explicit,R,0,0,...
                  'the public thermaliser agrees after the same relaxation cleanup');
rejected=false;
try
    thermalize(s,R0,[],[],target,'IME');
catch err
    rejected=strcmp(err.identifier,'Spinach:thermalize:targetConcentration');
end
result=test_true(result,'weighted target rejected',rejected,...
                 'a preweighted target must not silently acquire a second concentration');
result=test_true(result,'unit columns only',nnz(R(:,nonunit)-R0(:,nonunit))==0,...
                 'IME modifies no nonunit relaxation column');

% Collapse all local identities to one common coordinate without averaging
embedding=sparse(numel(nonunit)+1,s.bas.offsets(end));
embedding(1,units)=1;
embedding(2:end,nonunit)=speye(numel(nonunit));
merged_h=embedding*H*embedding'; merged_r=embedding*R0*embedding';
merged_r(:,1)=embedding*sum(R(:,units),2);
residual=embedding*R-merged_r*embedding;
h_error=norm(embedding*H-merged_h*embedding,inf);
diss_error=norm(residual(:,nonunit),inf);
source_error=norm(residual(:,units),inf);
result=test_true(result,'T7 Hamiltonian intertwining',h_error==0,...
                 'collapsing the identity does not change coherent dynamics');
result=test_true(result,'T7 dissipative intertwining',diss_error==0,...
                 'the discrepancy has no nonunit columns');
result=test_true(result,'T7 source discrepancy',source_error>0,...
                 'one common identity cannot carry independent thermal source populations');

% Measure the longitudinal recovery of the common-source representation
merged_final=expm(full(20*merged_r))*(embedding*rho);
coil=coil_state(s,'Lz','all','exact'); merged_coil=embedding*coil;
direct_signal=real(coil'*evolved); merged_signal=real(merged_coil'*merged_final);
ratio=merged_signal/direct_signal;
result=test_true(result,'T7 factor two',abs(ratio-2)<1e-5,...
                 'these near-identical proton pairs show the prototype factor-two recovery');
fprintf('CWDM_T7 H=%.16g nonunit=%.16g units=%.16g direct=%.16g merged=%.16g ratio=%.16g\n',...
        h_error,diss_error,source_error,direct_signal,merged_signal,ratio);

% Enforce the same target trace contract in Zeeman-Liouville space
sys=struct('magnet',14.1,'isotopes',{{'1H'}}); inter=struct();
inter.chem.parts={1}; inter.chem.concs=0.3; inter.temperature=298;
inter.relaxation={'damp'}; inter.damp_rate=2;
inter.equilibrium='zero'; inter.rlx_keep='labframe';
bas=struct('formalism','zeeman-liouv','approximation',{{'none'}});
s=assume(test_spin_system(sys,inter,bas),'nmr');
R0=relaxation(s); target=equilibrium(s); one=s; one.chem.concs=1;
R=thermalize(s,R0,[],[],equilibrium(one),'IME');
result=test_close(result,'Zeeman IME stationary target',R*target,zeros(4,1),...
                  1e-14,0,'unit-trace targets thermalise concentration-weighted states');
for concentration=[0 0.3 2]
    one.chem.concs=concentration; rejected=false;
    try
        thermalize(s,R0,[],[],equilibrium(one),'IME');
    catch err
        rejected=strcmp(err.identifier,'Spinach:thermalize:targetConcentration');
    end
    result=test_true(result,'Zeeman weighted target rejected',rejected,...
                     'zero and nonunit target traces are rejected without normalisation');
end
fprintf('CWDM_ZEEMAN_TARGET stationary=%.16g rejected=%d\n',norm(R*target),rejected);

end


