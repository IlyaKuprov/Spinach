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
% Nottingham requires a single two-electron substance; segmented input
% is unsupported even when each substance has its own electron pair.
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
                 'Nottingham does not support a separate nucleus-only substance');
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

% Reject compiled two-pair Nottingham input regardless of global spin ordering
local_sys=struct('magnet',1,'isotopes',{{'E','E','E','E'}});
local_bas=bas; local_bas.approximation={'none','none'};
partitions={{1:2,3:4},{[1 3],[2 4]}};
for n=1:numel(partitions)
    local_inter=struct(); local_inter.chem.parts=partitions{n};
    local_inter.chem.concs=[1 1];
    multi=assume(test_spin_system(local_sys,local_inter,local_bas),'nmr');
    multi.rlx=s.rlx; rejected=false;
    try
        relaxation(multi);
    catch err
        rejected=strcmp(err.identifier,'Spinach:relaxation:nottinghamSubstance');
    end
    result=test_true(result,sprintf('Nottingham two pairs %d',n),rejected,...
                     'the relaxation consumer rejects segmented Nottingham descriptors');
    fprintf('CWDM_NOTTINGHAM two_pairs_order=%d named_rejection=%d\n',n,rejected);
end

% Compare two thermalised steady states against independently built substances
sys=struct('magnet',1,'isotopes',{{'1H','13C'}});
inter=struct(); inter.chem.parts={1,2}; inter.chem.concs=[0.7 0.3];
inter.relaxation={'t1_t2'}; inter.r1_rates={1,2}; inter.r2_rates={3,4};
inter.equilibrium='IME'; inter.temperature=298; inter.rlx_keep='secular';
bas=struct('formalism','sphten-liouv','approximation',{{'none','none'}});
s=assume(test_spin_system(sys,inter,bas),'nmr');
P=expm(full(relaxation(s))); units=s.bas.offsets(1:end-1)+1;
for method={'newton','squaring'}
    reference=cell(2,1);
    for n=1:2
        local_sys=sys; local_sys.isotopes=sys.isotopes(n);
        local_inter=inter; local_inter.chem.parts={1}; local_inter.chem.concs=s.chem.concs(n);
        local_inter.r1_rates=inter.r1_rates(n); local_inter.r2_rates=inter.r2_rates(n);
        local_bas=bas; local_bas.approximation={'none'};
        local=assume(test_spin_system(local_sys,local_inter,local_bas),'nmr');
        reference{n}=steady(local,expm(full(relaxation(local))),[],method{1});
    end
    actual=steady(s,P,[],method{1}); reference=vertcat(reference{:});
    result=test_close(result,['steady direct sum ' method{1}],actual,reference,1e-12,1e-12,...
                      'both substance blocks equal their independently concentration-weighted steady states');
    result=test_true(result,['steady units ' method{1}],all(actual(units)==s.chem.concs(:)),...
                     'each substance carries its specified concentration');
    guess=full(unit_state(s)); guess(3)=0.1; guess(7)=-0.2;
    result=test_close(result,['steady initial guess ' method{1}],...
                      steady(s,P,guess,method{1}),reference,1e-12,1e-12,...
                      'a non-equilibrium initial guess converges with both unit constraints');
    fprintf('CWDM_STEADY %s error=%.16g units=%s\n',...
            method{1},norm(actual-reference),mat2str(actual(units)'));
end

% Reject either unthermalised block even when its partner is thermalised
for n=1:2
    block=(s.bas.offsets(n)+1):s.bas.offsets(n+1);
    bad=P; bad(block,block)=eye(numel(block));
    for method={'newton','squaring'}
        rejected=false;
        try
            steady(s,bad,[],method{1});
        catch err
            rejected=strcmp(err.identifier,'Spinach:steady:unthermalisedSubstance')&&...
                     contains(err.message,sprintf('substance %d',n));
        end
        result=test_true(result,sprintf('steady mixed %d %s',n,method{1}),rejected,...
                         'every substance must have its own thermalisation source');
        fprintf('CWDM_STEADY_MIXED block=%d method=%s named_rejection=%d\n',...
                n,method{1},rejected);
    end
end

% Reject segmented NGCE before its single-unit projection can couple substances
H0=1e-3*operator(s,'Lz','all'); H1=repmat({sparse(size(H0,1),size(H0,2))},2001,1);
for rate=[0 2]
    rejected=false;
    try
        ngce(s,H0,H1,1,10,rate);
    catch err
        rejected=strcmp(err.identifier,'Spinach:ngce:segmentedSubstances');
    end
    result=test_true(result,sprintf('NGCE segmented rate %g',rate),rejected,...
                     'NGCE rejects multiple substances with and without regularisation');
    fprintf('CWDM_NGCE segmented_reg=%g named_rejection=%d\n',rate,rejected);
end

% Retain regularised NGCE for a supported single-substance zero trajectory
H0=1e-3*operator(local,'Lz','all');
H1=repmat({sparse(size(H0,1),size(H0,2))},2001,1);
[R,dR]=ngce(local,H0,H1,1,10,2);
reference=-2*speye(size(H0)); reference(1,1)=0;
result=test_close(result,'NGCE single regularisation',R,reference,0,1e-14,...
                  'regularisation damps active states but leaves the single unit direction undamped');
result=test_close(result,'NGCE single zero uncertainty',dR,sparse(size(H0,1),size(H0,2)),0,0,...
                  'a zero stochastic trajectory has zero uncertainty');
fprintf('CWDM_NGCE single_regularisation_error=%.16g uncertainty=%.16g\n',...
        norm(R-reference,'fro'),norm(dR,'fro'));

% Refuse loss of trace conservation or normalisation in a later substance
bad=P; bad(units(2),units(2)+1)=0.1; rejected=false;
try
    steady(s,bad,[],'newton');
catch err
    rejected=contains(err.message,'conserve every substance unit coordinate');
end
result=test_true(result,'steady later trace row',rejected,...
                 'every conserved unit row is validated, not only the first');
guess=full(unit_state(s)); guess(units(2))=0; rejected=false;
try
    steady(s,P,guess,'newton');
catch err
    rejected=contains(err.message,'every substance unit coordinate of rho');
end
result=test_true(result,'steady later normalisation',rejected,...
                 'an initial guess must carry the specified substance concentrations');

% Drive each pumped state from only its own substance population
R=sparse(s.bas.offsets(end),s.bas.offsets(end));
rho=coil_state(s,'Lz',2,'exact'); pumped=magpump(s,R,rho,2);
reference=R; reference(:,units(2))=2*rho;
result=test_close(result,'pump second substance',pumped,reference,0,0,...
                  'a spin selected in the second substance is sourced by its own unit column');
source=full(unit_state(s)); source(units)=[0.2;0.8];
result=test_close(result,'pump unequal populations',pumped*source,1.6*rho,0,1e-14,...
                  'the second target is driven by the second population, not the first');
rho=coil_state(s,'Lz',1,'exact')+coil_state(s,'Lz',2,'exact'); pumped=magpump(s,R,rho,2);
reference=2*(0.2*coil_state(s,'Lz',1,'exact')+0.8*coil_state(s,'Lz',2,'exact'));
result=test_close(result,'pump both substances',pumped*source,reference,0,1e-14,...
                  'a state spanning several substances is sourced independently in each block');
rejected=false; rho(units(2))=1;
try
    magpump(s,R,rho,2);
catch err
    rejected=strcmp(err.message,'unit state cannot be pumped.');
end
result=test_true(result,'pump later identity rejection',rejected,...
                 'an identity component in any substance is rejected');
fprintf('CWDM_MAGPUMP population_error=%.16g later_identity_rejected=%d\n',...
        norm(pumped*source-reference),rejected);

% Reject the scalar solid-effect unit projection for independent substances
sys=struct('magnet',1,'isotopes',{{'E','1H'}});
inter=struct(); inter.chem.parts={1,2}; inter.chem.concs=[.7 .3];
inter.relaxation={'t1_t2'}; inter.r1_rates={10,1}; inter.r2_rates={20,2};
inter.rlx_keep='secular'; inter.equilibrium='zero'; inter.temperature=298;
bas=struct('formalism','sphten-liouv','approximation',{{'none','none'}});
s=test_spin_system(sys,inter,bas);
parameters=struct('mw_pwr',0,'theory','exact','nuclear_frq',1,...
                  'calc_type','steady_state'); rejected=false;
try
    solid_effect(s,parameters);
catch err
    rejected=strcmp(err.identifier,'Spinach:solid_effect:segmentedSubstances');
end
result=test_true(result,'segmented solid effect',rejected,...
                 'steady-state solid effect rejects multiple trace null directions');
fprintf('CWDM_SOLID_EFFECT named_rejection=%d\n',rejected);

% Retain the complete supported single-substance steady-state experiment
inter.chem.parts={1:2}; inter.chem.concs=1; bas.approximation={'none'};
s=test_spin_system(sys,inter,bas); actual=solid_effect(s,parameters);
[I,Q]=hamiltonian(assume(s,'labframe'),'left'); rho=equilibrium(s,I,Q,[0 0 0]);
coils=[state(s,'Lz',1) state(s,'Lz',2)];
result=test_close(result,'single solid effect zero drive',actual,coils'*rho,...
                  1e-12,1e-12,'without microwave drive the steady state is thermal');

% Reject converted segmented Zeeman steady solves before default initialisation
sys=struct('magnet',1,'isotopes',{{'1H','1H'}});
inter=struct(); inter.chem.parts={1,2}; inter.chem.concs=[1 1];
bas=struct('formalism','zeeman-hilb','approximation',{{'none','none'}});
s=test_spin_system(sys,inter,bas);
[s,~,~]=sim2liouv(s,struct(),sparse(4,4),[],[]);
for method={'newton','squaring'}
    for guess={[],ones(8,1)}
        rejected=false;
        try
            steady(s,speye(8),guess{1},method{1});
        catch err
            rejected=strcmp(err.identifier,'Spinach:steady:segmentedZeeman');
        end
        result=test_true(result,['Zeeman steady ' method{1} ' ' int2str(numel(guess{1}))],...
                         rejected,'segmented Zeeman initialisation raises the named boundary');
        fprintf('CWDM_ZEEMAN_STEADY %s guess=%d named_rejection=%d\n',...
                method{1},numel(guess{1}),rejected);
    end
end

% DNP scans reject multiple independent trace null directions before solving
sys=struct('magnet',1,'isotopes',{{'E','1H'}});
inter=struct(); inter.chem.parts={1,2}; inter.chem.concs=[.7 .3];
inter.relaxation={'t1_t2'}; inter.r1_rates={10,1}; inter.r2_rates={20,2};
inter.rlx_keep='secular'; inter.equilibrium='zero';
bas=struct('formalism','sphten-liouv','approximation',{{'none','none'}});
s=assume(test_spin_system(sys,inter,bas),'esr');
H=hamiltonian(s); R=relaxation(s); K=sparse(size(R,1),size(R,2));
parameters=struct('mw_pwr',0,'mw_frq',0,'g_ref',s.tols.freeg,...
                  'rho0',unit_state(s),'coil',state(s,'Lz',2),...
                  'mw_oper',operator(s,'Lx',1),'ez_oper',operator(s,'Lz',1),...
                  'method','lvn-backs','fields',0);
for scan={@dnp_freq_scan,@dnp_field_scan}
    if isequal(scan{1},@dnp_field_scan), parameters.method='backslash'; end
    rejected=false;
    try
        scan{1}(s,parameters,H,R,K);
    catch err
        rejected=strcmp(err.identifier,['Spinach:' func2str(scan{1}) ':segmentedSubstances']);
    end
    result=test_true(result,['segmented ' func2str(scan{1})],rejected,...
                     'the single-trace scan algorithm rejects segmented input');
    fprintf('CWDM_DNP_SCAN %s named_rejection=%d\n',func2str(scan{1}),rejected);
end

% A supported single-substance field scan retains its zero-drive equilibrium
inter.chem.parts={1:2}; inter.chem.concs=1; bas.approximation={'none'};
s=assume(test_spin_system(sys,inter,bas),'esr');
H=hamiltonian(s); R=relaxation(s); K=sparse(size(R,1),size(R,2));
parameters.rho0=unit_state(s)+.1*state(s,'Lz',2);
parameters.coil=state(s,'Lz',2); parameters.mw_oper=operator(s,'Lx',1);
parameters.ez_oper=operator(s,'Lz',1); parameters.fields=[-.001 0 .001];
actual=dnp_field_scan(s,parameters,H,R,K);
expected=repmat(parameters.coil'*parameters.rho0,3,1);
result=test_close(result,'single field scan zero drive',actual,expected,...
                  1e-12,1e-12,'longitudinal equilibrium survives a zero-drive field scan');

% A spin-free pool has no active rows and must not be mistaken for an
% unthermalised substance: its steady state is its unit coordinate
sys=struct('magnet',1,'isotopes',{{'1H'}});
inter=struct(); inter.chem.parts={1,[]}; inter.chem.concs=[0.7 0.3];
inter.relaxation={'t1_t2'}; inter.r1_rates={1}; inter.r2_rates={3};
inter.equilibrium='IME'; inter.temperature=298; inter.rlx_keep='secular';
bas=struct('formalism','sphten-liouv','approximation',{{'none','none'}});
pool=assume(test_spin_system(sys,inter,bas),'nmr');
P=expm(full(relaxation(pool))); units=pool.bas.offsets(1:end-1)+1;
for method={'newton','squaring'}
    actual=steady(pool,P,[],method{1});
    result=test_true(result,['steady spin-free pool ' method{1}],...
                     all(actual(units)==pool.chem.concs(:))&&(numel(actual)==pool.bas.offsets(end)),...
                     'a spin-free pool keeps its unit coordinate without a thermalisation source');
    fprintf('CWDM_STEADY_POOL method=%s units=%s\n',method{1},mat2str(actual(units)'));
end

% IME requires independent relaxation blocks even when cross terms preserve trace
sys=struct('magnet',1,'isotopes',{{'1H','1H'}});
inter=struct(); inter.chem.parts={1,2}; inter.chem.concs=[1 1];
bas=struct('formalism','sphten-liouv','approximation',{{'none','none'}});
s=test_spin_system(sys,inter,bas); unit=unit_state(s);
rho_eq=unit+.1*state(s,'Lz',1)+.3*state(s,'Lz',2);
R=-speye(s.bas.offsets(end)); units=s.bas.offsets(1:end-1)+1;
R(units,units)=0; cross=R; cross(3,7)=.2; cross(7,3)=.2;
result=test_true(result,'cross relaxation zero unit action',norm(cross*unit)==0,...
                 'the fixture passes the existing already-thermalised check');
rejected=false;
try
    thermalize(s,cross,[],[],rho_eq,'IME');
catch err
    rejected=strcmp(err.identifier,'Spinach:thermalize:crossSubstanceRelaxation');
end
result=test_true(result,'cross relaxation IME rejection',rejected,...
                 'a trace-preserving cross-substance relaxation block is unsupported');
fprintf('CWDM_IME_CROSS named_rejection=%d unit_action=%.16g\n',rejected,norm(cross*unit));

% Supported IME direct sums annihilate the requested equilibrium at round-off
actual=thermalize(s,R,[],[],rho_eq,'IME');
result=test_close(result,'IME equilibrium stationary',actual*rho_eq,zeros(size(rho_eq)),...
                  10*eps*norm(actual,'fro')*norm(rho_eq),0,...
                  'independent IME blocks annihilate the unweighted target state');
fprintf('CWDM_IME_STATIONARY residual=%.16g\n',norm(actual*rho_eq));

end


