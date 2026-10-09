% Tests substance-local relaxation against independently built spin systems.
% Syntax:
%
%                       result=test_cwdm_relaxation()
%
% Outputs:
%
%     result - direct-sum Hamiltonian and relaxation checks
%
% Two different heteronuclear pairs carry CSA and dipolar interactions,
% distinct correlation times, and distinct phenomenological rates. Full
% matrices, including every unit row and column, are compared at 1e-12
% relative tolerance. Both Redfield integration paths are exercised.
% Nottingham requires a two-electron manifold in every substance; a
% nucleus-only partner is an unsupported specification, not a zero block.
% There is no relaxation theory named weak in the supported input API.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_relaxation()

% Announce the direct-sum reference
fprintf('TESTING: CWDM relaxation direct sums\n');
result=new_test_result('kernel/cwdm_relaxation','CWDM relaxation direct sums',...
                      'Each relaxation block equals an independent substance.');

% Specify two physically distinct pairs without cross-substance interactions
sys.magnet=1; sys.isotopes={'1H','13C','19F','15N'};
inter.chem.parts={1:2,3:4}; inter.chem.concs=[0.7 0.3];
inter.zeeman.matrix={diag([1 2 6]),diag([10 20 30]),...
                     diag([-4 3 8]),diag([15 25 50])};
inter.coupling.matrix=cell(4);
inter.coupling.matrix{1,2}=diag([-100 -100 200]);
inter.coupling.matrix{3,4}=diag([-200 -200 400]);
inter.rlx_keep='secular'; inter.equilibrium='zero';
bas.formalism='sphten-liouv'; bas.approximation={'none','none'};
theories={{'redfield'},{'redfield'},{'lindblad'},{'SRFK'},{'t1_t2','redfield'}};
policies={'kite','secular','secular','secular','secular'};

% Build every reference from its own physical input rather than slicing R
for n=1:numel(theories)
    current=inter; current.relaxation=theories{n};
    current.rlx_keep=policies{n};
    if ismember('redfield',theories{n})
        current.tau_c={1e-9,3e-9};
    end
    if ismember('t1_t2',theories{n})
        current.r1_rates={1,2,3,4}; current.r2_rates={3,4,5,6};
    end
    if ismember('lindblad',theories{n})
        current.lind_r1_rates=[1 2 3 4]; current.lind_r2_rates=[3 4 5 6];
    end
    if ismember('SRFK',theories{n})
        current.srfk_tau_c={[1 1e-6]}; current.srfk_mdepth=cell(4);
        current.srfk_mdepth{1,2}=100; current.srfk_mdepth{3,4}=200;
    end
    s=assume(test_spin_system(sys,current,bas),'nmr');
    R=cell(1,2); H=cell(1,2);
    for k=1:2
        spins=inter.chem.parts{k}; local_sys=sys;
        local_sys.isotopes=sys.isotopes(spins);
        local_inter=current; local_inter.chem.parts={1:2};
        local_inter.chem.concs=1;
        local_inter.zeeman.matrix=current.zeeman.matrix(spins);
        local_inter.coupling.matrix=current.coupling.matrix(spins,spins);
        if isfield(current,'tau_c'), local_inter.tau_c=current.tau_c(k); end
        if isfield(current,'r1_rates')
            local_inter.r1_rates=current.r1_rates(spins);
            local_inter.r2_rates=current.r2_rates(spins);
        end
        if isfield(current,'lind_r1_rates')
            local_inter.lind_r1_rates=current.lind_r1_rates(spins);
            local_inter.lind_r2_rates=current.lind_r2_rates(spins);
        end
        if isfield(current,'srfk_mdepth')
            local_inter.srfk_mdepth=current.srfk_mdepth(spins,spins);
        end
        local_bas=bas; local_bas.approximation={'none'};
        local=assume(test_spin_system(local_sys,local_inter,local_bas),'nmr');
        R{k}=relaxation(local); H{k}=hamiltonian(local);
    end
    reference=blkdiag(R{:}); actual=relaxation(s);
    label=[strjoin(theories{n},'+') ' ' policies{n}];
    rel_error=norm(actual-reference,'fro')/norm(reference,'fro');
    result=test_true(result,label,isfinite(rel_error)&&rel_error<=1e-12,...
                     'the complete matrix equals the independent direct sum');
    fprintf('CWDM_RELAXATION %s relative_error=%.16g\n',label,rel_error);
    units=s.bas.offsets(1:end-1)+1;
    result=test_true(result,[label ' units'],...
                     nnz(actual(units,:))==0&&nnz(actual(:,units))==0,...
                     'every unthermalised unit row and column is exactly zero');
    if n==1
        result=test_close(result,'Hamiltonian direct sum',hamiltonian(s),...
                          blkdiag(H{:}),0,1e-12,...
                          'the same two substances have independent Hamiltonians');
    end
    if ismember('redfield',theories{n})
        s.sys.disable=unique([s.sys.disable {'asyredf'}]);
        actual=relaxation(s);
        rel_error=norm(actual-reference,'fro')/norm(reference,'fro');
        result=test_true(result,[label ' serial'],...
                         isfinite(rel_error)&&rel_error<=1e-12,...
                         'serial integration preserves the same independent blocks');
        fprintf('CWDM_RELAXATION %s serial_relative_error=%.16g\n',label,rel_error);
    end
end

% Record the absence of a weak relaxation model in the input contract
bad=inter; bad.relaxation={'weak'}; rejected=false;
try
    test_spin_system(sys,bad,bas);
catch err
    rejected=contains(err.message,'unrecognised relaxation theory');
end
result=test_true(result,'weak unavailable',rejected,...
                 'weak is not an implemented relaxation theory');

% Exercise the Nottingham two-electron manifold with a separate nucleus
sys.isotopes={'E','E','1H'}; nott.chem.parts={1:2,3};
nott.chem.concs=[1 1]; nott.relaxation={'nottingham'};
nott.rlx_keep='secular'; nott.equilibrium='zero';
nott.nott_r1e=100; nott.nott_r2e=1000;
nott.nott_r1n=1; nott.nott_r2n=2;
s=assume(test_spin_system(sys,nott,bas),'nmr'); rejected=false;
try
    relaxation(s);
catch err
    rejected=strcmp(err.identifier,'Spinach:relaxation:nottinghamSubstance');
end
result=test_true(result,'Nottingham nucleus-only substance',rejected,...
                 'every Nottingham substance must contain its own two-electron manifold');
fprintf('CWDM_NOTTINGHAM nucleus_only_named_rejection=%d\n',rejected);

% Reject electrons split across substances before building their products
nott.chem.parts={1,[2 3]};
s=assume(test_spin_system(sys,nott,bas),'nmr'); rejected=false;
try
    relaxation(s);
catch err
    rejected=strcmp(err.identifier,'Spinach:relaxation:nottinghamSubstance');
end
result=test_true(result,'Nottingham split electrons',rejected,...
                 'the two-electron manifold cannot span substances');

% Retain the ordinary single-substance Nottingham theory
nott.chem.parts={1:3}; nott.chem.concs=1; bas.approximation={'none'};
s=assume(test_spin_system(sys,nott,bas),'nmr'); R=relaxation(s);
result=test_true(result,'Nottingham single substance',...
                 norm(R,'fro')>0&&nnz(R(1,:))==0&&nnz(R(:,1))==0,...
                 'a supported two-electron substance retains nonzero trace-preserving relaxation');
fprintf('CWDM_NOTTINGHAM single_substance_norm=%.16g\n',norm(R,'fro'));

end


