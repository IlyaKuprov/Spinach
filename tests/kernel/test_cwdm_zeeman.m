% Tests Zeeman direct-sum operators and thermal states (T16 and T17).
% Syntax:
%
%                        result=test_cwdm_zeeman()
%
% Outputs:
%
%    result - independent representation and concentration checks
%
% Zeeman T1/T2 construction remains unsupported: the bare spherical-tensor
% relaxation matrix is explicitly transformed before production IME is
% applied. The measured formalism probe motivates absolute bounds 1e-14
% for equilibrium, 1e-10 for spectral-norm intertwining errors, and 1e-12
% for the concentration-generator identity. Hilbert matrix exchange and
% independent-reactant tensor products use algebraic references.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_zeeman()

% Announce independent representation checks
result=new_test_result('kernel/cwdm_zeeman','CWDM Zeeman T16-T17',...
                       'Local Zeeman blocks preserve thermal states, spin dynamics, and traces.');

% Build the two-species proton-pair formalism-probe fixture
sys.magnet=14.1; sys.isotopes={'1H','1H','1H','1H'};
inter.temperature=298; inter.zeeman.scalar={1,2,3,4};
inter.coupling.scalar=cell(4); inter.coupling.scalar{1,2}=10;
inter.coupling.scalar{3,4}=8; inter.chem.parts={1:2,3:4};
inter.chem.concs=[0.4 0.6]; inter.relaxation={'t1_t2'};
inter.chem.reactions={struct('reactants',1,'products',2,'matching',[1 3;2 4],'rate',20),...
                      struct('reactants',2,'products',1,'matching',[3 1;4 2],'rate',10)};
inter.r1_rates={1,1.5,1,1.5}; inter.r2_rates={3,4,3,4};
inter.rlx_keep='labframe'; inter.equilibrium='IME';
bas.formalism='sphten-liouv'; bas.approximation={'none','none'};
s=test_spin_system(sys,inter,bas); s=assume(s,'nmr');
bas.formalism='zeeman-liouv'; z=test_spin_system(sys,inter,bas); z=assume(z,'nmr');

% Convert spherical-tensor coefficients to trace-one physical coordinates
P=full(sphten2zeeman(s))/4;
rho_s=equilibrium(s); rho_z=equilibrium(z);
H=hamiltonian(s); hz=hamiltonian(z);
R=relaxation(s); bare=s; bare.rlx.equilibrium='zero';
bare=relaxation(bare); unit_z=z; unit_z.chem.concs(:)=1;
rz=thermalize(z,P*bare/P,[],[],equilibrium(unit_z),'IME');
eq_error=norm(rho_z-P*rho_s); h_error=norm(hz*P-P*H);
r_error=norm(rz*P-P*R);
result=test_close(result,'T16 equilibrium',rho_z,P*rho_s,1e-14,0,...
                  'the physical converter includes the local inverse Hilbert dimension');
result=test_close(result,'T16 Hamiltonian',h_error,0,1e-10,0,...
                  'independently assembled Hamiltonians intertwine');
result=test_close(result,'T16 thermalised relaxation',r_error,0,1e-10,0,...
                  'production Zeeman IME matches transformed spherical-tensor IME');

% Test the probe exchange on every component including the trace
exchange=[-20 10;20 -10]; K=kinetics(z); ks=kinetics(s);
result=test_close(result,'T16 production exchange map',K,kron(exchange,speye(16)),1e-12,0,...
                  'reaction records transport every local Zeeman density-matrix component');
tau=kron(speye(2),reshape(speye(4),1,16));
generator=-1i*hz+rz+K; generator_error=norm(tau*generator-exchange*tau);
result=test_close(result,'T16 concentration generator',tau*generator,exchange*tau,1e-12,0,...
                  'spin dynamics and IME preserve traces while exchange transports them');
result=test_close(result,'T16 weighted traces',tau*rho_z,inter.chem.concs',1e-14,0,...
                  'thermal block traces are concentrations');
fprintf('CWDM_T16 eq=%.16g H=%.16g R=%.16g generator=%.16g\n',...
        eq_error,h_error,r_error,generator_error);

% Compare pulse excitation and dual detection without empirical rescaling
pulse_s=operator(s,'Ly','1H'); pulse_z=operator(z,'Ly','1H');
rho_s=expm(full(-1i*pulse_s*pi/2))*rho_s;
rho_z=expm(full(-1i*pulse_z*pi/2))*rho_z;
coil_s=coil_state(s,'L+','1H','exact'); coil_z=coil_state(z,'L+','1H','exact');
prop_s=expm(full(-1i*H+R+ks)*0.0002); prop_z=expm(full(generator)*0.0002);
fid_s=zeros(256,1); fid_z=zeros(256,1);
for n=1:256
    fid_s(n)=coil_s'*rho_s; fid_z(n)=coil_z'*rho_z;
    rho_s=prop_s*rho_s; rho_z=prop_z*rho_z;
end
result=test_close(result,'T16 FID',fid_z,fid_s,1e-12,0,...
                  'unweighted coils detect the same concentration-weighted signal');
fprintf('CWDM_T16_FID error=%.16g peak=%.16g\n',max(abs(fid_z-fid_s)),max(abs(fid_z)));

% Compare Hilbert thermal blocks and block-aware conversion
inter=rmfield(inter,{'r1_rates','r2_rates'});
inter.relaxation={}; bas.formalism='zeeman-hilb';
h=test_spin_system(sys,inter,bas); rho=equilibrium(h);
blocks={rho(1:4,1:4),rho(5:8,5:8)};
result=test_close(result,'T17 Hilbert thermal conversion',hilb2liouv(blocks,'statevec'),...
                  equilibrium(z),1e-14,0,'Hilbert and Zeeman-Liouville equilibria agree block-wise');
result=test_close(result,'T17 Hilbert traces',cellfun(@trace,blocks),inter.chem.concs,1e-14,0,...
                  'each Hilbert block has its own partition function and concentration');
H=hamiltonian(assume(h,'nmr'));
reference=blkdiag(kron(speye(4),H(1:4,1:4))-kron(H(1:4,1:4).',speye(4)),...
                  kron(speye(4),H(5:8,5:8))-kron(H(5:8,5:8).',speye(4)));
result=test_close(result,'block commutator conversion',...
                  hilb2liouv({H(1:4,1:4),H(5:8,5:8)},'comm'),reference,0,0,...
                  'cell conversion is the direct sum of local commutator Kronecker maps');

% Exercise strongly polarised exchange and the independent-reactant product
rho_a=diag([0.7 0.3])*0.4; rho_b=diag([0.2 0.8])*0.6;
dot_a=-20*rho_a+10*rho_b; dot_b=20*rho_a-10*rho_b;
trace_error=abs(trace(dot_a+dot_b));
rate=7; dot_a=-rate*trace(rho_b)*rho_a; dot_b=-rate*trace(rho_a)*rho_b;
U=sparse([1 2 3 4],[1 3 2 4],1,4,4); dot_c=rate*U*kron(rho_a,rho_b)*U';
balance=[trace(dot_a)+trace(dot_c);trace(dot_b)+trace(dot_c)];
result=test_close(result,'T17 exchange trace',trace_error,0,1e-14,0,...
                  'first-order matrix exchange preserves total trace');
result=test_close(result,'T17 product atom balance',balance,zeros(2,1),1e-14,0,...
                  'permuted tensor-product arrival conserves each atom equivalent');
result=test_close(result,'T17 extent',trace(dot_c),1.68,1e-14,0,...
                  'strong polarisation does not change the mass-action extent');
fprintf('CWDM_T17 trace=%.16g atom_balance=%.16g extent=%.16g\n',...
        trace_error,norm(balance),trace(dot_c));

% Include unequal blocks, a spin-free state, and a zero concentration
sys.isotopes={'1H','1H','1H'}; inter=struct('temperature',298);
inter.zeeman.scalar={1,2,3}; inter.chem.parts={1,2:3,[]}; inter.chem.concs=[0.4 0 0.2];
bas.approximation={'none','none','none'};
for formalism={'zeeman-liouv','zeeman-hilb'}
    bas.formalism=formalism{1}; z=test_spin_system(sys,inter,bas);
    rho=equilibrium(z); coil=coil_state(z,'Lz',2,'exact');

    % Check product placement, identity support, and both sparse formats
    product=kron(diag([0.5 -0.5]),[0 0.5;0.5 0]);
    for kind={'left','right','comm','acomm'}
        expected={sparse(2,2),product,sparse(1,1)};
        if strcmp(formalism{1},'zeeman-liouv')
            expected=hilb2liouv(expected,kind{1});
        else
            expected=blkdiag(expected{:});
        end
        actual=operator(z,{'Lz','Lx'},{2,3},kind{1});
        result=test_close(result,['local product ' formalism{1} ' ' kind{1}],actual,expected,1e-14,0,...
                          'a product operator acts only on its hosting molecule');
        xyz=operator(z,{'Lz','Lx'},{2,3},kind{1},'xyz');
        result=test_close(result,['triplets ' formalism{1} ' ' kind{1}],...
                          sparse(xyz(:,1),xyz(:,2),xyz(:,3),size(actual,1),size(actual,2)),...
                          actual,0,0,'CSC and XYZ outputs use identical global offsets');
    end
    rejected=false;
    try
        operator(z,{'Lz','Lx'},{1,2});
    catch err
        rejected=strcmp(err.identifier,'Spinach:which_subst:crossSubstance');
    end
    result=test_true(result,['cross-substance product ' formalism{1}],rejected,...
                     'different molecules have no tensor-product spin operator');
    identity=operator(z,{'E'},{2},'left');
    idx=(z.bas.offsets(2)+1):z.bas.offsets(3);
    result=test_true(result,['identity support ' formalism{1}],...
                     nnz(identity)==numel(idx)&&isequal(identity(idx,idx),speye(numel(idx))),...
                     'an explicit identity stays within the selected substance');
    result=test_close(result,['zero-population state ' formalism{1}],state(z,'Lz',2),...
                      0*coil,0,0,'weighted state vanishes but its detection shape remains');
    result=test_true(result,['zero-population coil ' formalism{1}],nnz(coil)>0,...
                     'concentrations do not enter unweighted detection');
    units=unit_state(z);
    result=test_true(result,['compiled shape ' formalism{1}],size(units,1)==z.bas.offsets(end),...
                     'geometric units occupy the direct sum rather than the global tensor product');
    for n=1:3
        idx=(z.bas.offsets(n)+1):z.bas.offsets(n+1);
        dim=prod(z.comp.mults(z.chem.parts{n}));
        if strcmp(formalism{1},'zeeman-liouv')
            population=reshape(speye(dim),1,[])*rho(idx);
        else
            population=trace(rho(idx,idx));
        end
        result=test_close(result,[formalism{1} ' population ' num2str(n)],population,...
                          inter.chem.concs(n),1e-14,0,'spin-free and empty-population blocks retain their traces');
    end
end

% Keep molar equilibrium magnetisation within its single-substance contract
sys.magnet=1; sys.isotopes={'E','E','E','E'}; inter=struct();
inter.temperature=298; inter.zeeman.scalar={2,2,2,2};
inter.chem.parts={1:2,3:4}; inter.chem.concs=[0.4 0.6];
bas.formalism='zeeman-hilb'; bas.approximation={'none','none'};
z=test_spin_system(sys,inter,bas); parameters.grid='single_crystal';
rejected=false;
try
    eqmag(z,parameters);
catch err
    rejected=strcmp(err.identifier,'Spinach:eqmag:multipleSubstances')&&...
             strcmp(err.message,'eqmag requires a single substance; mixture molar normalisation is not supported.');
end
result=test_true(result,'eqmag mixture boundary',rejected,...
                 'two four-state molecules require an explicit mixture molar-normalisation contract');
inter.chem.parts={1:4}; inter.chem.concs=1; bas.approximation={'none'};
z=test_spin_system(sys,inter,bas); rho=equilibrium(z);
expected=zeros(1,3); labels={'Lx','Ly','Lz'};
for n=1:3
    expected(n)=-2*real(trace(rho*operator(z,labels{n},'E')))/real(trace(rho));
end
result=test_close(result,'eqmag singleton',eqmag(z,parameters),expected,1e-12,0,...
                  'the isotropic single-substance magnetisation retains the explicit thermal trace');
fprintf('CWDM_ZEEMAN_FAILURES %d\n',numel(result.failures));
end


