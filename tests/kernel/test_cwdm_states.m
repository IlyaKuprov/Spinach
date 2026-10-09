% Tests concentration-weighted states and unweighted detection vectors (T6).
% Syntax:
%
%                       result=test_cwdm_states()
%
% Outputs:
%
%     result - exact state-weighting and thermal-equilibrium checks
%
% Independent single-substance constructions supply the reference shapes.
% Exact and cheap states, identities, level projectors, zero populations,
% and spin-free substances are covered. Thermal vectors use a 1e-12 bound.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_states()

% Announce the concentration contract
fprintf('TESTING: CWDM weighted states (T6)\n');
result=new_test_result('kernel/cwdm_states','CWDM weighted states',...
                      'States carry concentrations; coils carry operator shapes.');

% Build two proton pairs and a spin-free substance
sys.magnet=14.1; sys.isotopes={'1H','1H','1H','1H'};
inter.chem.parts={1:2,3:4,[]}; inter.chem.concs=[0.7 0.3 0];
inter.zeeman.scalar={1,2,3,4}; inter.coupling.scalar=cell(4);
inter.coupling.scalar{1,2}=10; inter.coupling.scalar{3,4}=8;
inter.temperature=298;
bas.formalism='sphten-liouv'; bas.approximation={'none','none','none'};
s=test_spin_system(sys,inter,bas); units=s.bas.offsets(1:end-1)+1;

% Reject malformed caller descriptions at the public state boundary
for args={{s,'Lz',[1 1],'exact'},{s,{'Lz','Lx'},{1},'exact'},...
          {s,'Lz',1,'unknown'},{s,[0 0],[],'exact'}}
    rejected=false;
    try
        state(args{1}{:});
    catch err
        fprintf('WRAPPER_REJECTION %s %s\n',err.stack(1).name,err.message);
        rejected=strcmp(err.stack(1).name,'grumble')&&...
                 strcmp(err.stack(1).file,which('state'));
    end
    result=test_true(result,'wrapper-local argument rejection',rejected,...
                     'invalid descriptions are rejected by the public wrapper grumbler');
end

% Compare operator shapes with independent unit-concentration substances
reference=cell(3,1);
for n=1:2
    spins=inter.chem.parts{n}; rows=(s.bas.offsets(n)+1):s.bas.offsets(n+1);
    local_sys=sys; local_sys.isotopes=sys.isotopes(spins);
    local_inter=inter; local_inter.chem.parts={1:2}; local_inter.chem.concs=1;
    local_inter.zeeman.scalar=inter.zeeman.scalar(spins);
    local_inter.coupling.scalar=inter.coupling.scalar(spins,spins);
    local_bas=bas; local_bas.approximation={'none'};
    local=test_spin_system(local_sys,local_inter,local_bas);
    reference{n}=s.chem.concs(n)*equilibrium(local);
    for method={'exact','cheap'}
        for label={'Lz','L+','E','ZL1'}
            expected=zeros(s.bas.offsets(end),1);
            expected(rows)=coil_state(local,label{1},1,method{1});
            actual=coil_state(s,label{1},spins(1),method{1});
            name=sprintf('%s %s block %d',method{1},label{1},n);
            result=test_true(result,['coil ' name],isequal(actual,expected),...
                             'an independent molecule supplies the same unweighted shape');
            result=test_true(result,['state ' name],...
                             isequal(state(s,label{1},spins(1),method{1}),s.chem.concs(n)*actual),...
                             'concentration multiplies the operator shape exactly once');
        end
    end
end
reference{3}=s.chem.concs(3);
result=test_close(result,'thermal direct sum',equilibrium(s),vertcat(reference{:}),...
                  1e-12,0,'equilibrium equals independently weighted thermal blocks');
result=test_true(result,'unit coordinates',...
                 isequal(full(unit_state(s)),full(sparse(units,1,s.chem.concs,s.bas.offsets(end),1))),...
                 'the identity contains exactly one concentration per substance');

% Check an all-spin sum is weighted separately in each hosting block
coil=coil_state(s,'Lz','all','exact'); expected=full(coil);
for n=1:s.bas.nsubst
    rows=(s.bas.offsets(n)+1):s.bas.offsets(n+1);
    expected(rows)=s.chem.concs(n)*expected(rows);
end
result=test_true(result,'all-spin weighting',isequal(state(s,'Lz','all'),expected),...
                 'an isotope sum does not share one concentration across blocks');

% Retain finite geometric coils when both spin-bearing populations vanish
empty=s; empty.chem.concs=[0 0 1];
result=test_true(result,'zero-population state',nnz(state(empty,'Lz','all'))==0,...
                 'absent substances carry no spin order');
result=test_true(result,'zero-population coil',isequal(coil_state(empty,'Lz','all','exact'),coil),...
                 'detection is independent of concentration');
result=test_true(result,'spin-free equilibrium',isequal(equilibrium(empty),unit_state(empty)),...
                 'the populated spin-free substance carries only its unit coordinate');
fprintf('CWDM_T6 exact_weighting=1 zero_population=1 spin_free=1 thermal_error=%.16g\n',...
        norm(equilibrium(s)-vertcat(reference{:}),inf));

% Keep storage-only wavefunction probabilities independent of concentration
sys=struct('magnet',0,'isotopes',{{'1H'}}); inter=struct();
inter.chem.parts={1};
bas=struct('formalism','zeeman-wavef','approximation',{{'none'}});
for concentration=[0 0.3 2]
    inter.chem.concs=concentration;
    s=test_spin_system(sys,inter,bas); psi=state(s,0.5);
    result=test_true(result,'unweighted wavefunction',isequal(psi,[1;0]),...
                     'storage-only kets retain unit probability at every concentration');
end
fprintf('CWDM_WAVEFUNCTION norm_squared=%.16g\n',norm(psi)^2);

% Require all four arguments on the new unweighted primitive
rejected=false;
try
    coil_state(s,0.5,[]);
catch err
    rejected=strcmp(err.identifier,'MATLAB:minrhs');
end
result=test_true(result,'fixed coil signature',rejected,...
                 'coil_state has no implicit method or spin-list defaults');
result=test_true(result,'explicit wavefunction coil',...
                 isequal(coil_state(s,0.5,[],'exact'),[1;0]),...
                 'the explicit four-argument wavefunction API remains supported');

end


