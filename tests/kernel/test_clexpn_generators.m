% Test: kernel/clexpn_generators | Native CCE generators
% Tests native generators, branch propagation, dissipation, and finite pulses.
% Syntax:
%
%                    result=test_clexpn_generators()
%
% Outputs:
%
%     result - regression result with recorded checks and failures
%
% ilya.kuprov@weizmann.ac.il

function result=test_clexpn_generators()

% Keep a deterministic retained bath with nonconsecutive physical-spin order
result=new_test_result('kernel/clexpn_generators','Native CCE generators',...
                       'Generator signs, projected branches and propagation.');
sys.isotopes={'1H','13C','14N'};
sys.magnet=0;
sys.output='hush';
sys.disable={'hygiene'};
inter.zeeman.matrix={diag([0 0 11]),diag([0 0 5]),diag([0 0 3])};
inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=diag([0 0 7]);
inter.coupling.matrix{1,3}=diag([0 0 2]);
inter.coupling.matrix{2,3}=diag([1 2 4]);
spin_system=create(sys,inter);
p.observed=1;
p.times=[0 0.013];
p.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
p.cce.order=2;
p.cce.numerics.max_memory_bytes=2^28;
p.cce.numerics.atol=1e-12;
p.cce.numerics.rtol=1e-10;
[p,~]=cce_options(spin_system,p,@sin,{'labframe'});
model=cce_compile(spin_system,p,{'labframe'});
[cluster,~]=cce_subsystem(model,[3 2],struct(),[]);
[H,R,K,aux,cache]=cce_generators(cluster,p,[]);
result=test_close(result,'Hilbert relaxation zero',norm(R)+norm(K),0,0,0,...
                  'A coherent Hilbert generator has no dissipative contribution.');
[H2,~,~,~,cache]=cce_generators(cluster,aux,cache);
result=test_close(result,'Exact generator cache',H2,H,0,0,...
                  'Identical retained inputs must return an identical generator.');
result=test_true(result,'Generator cache hit',cache.hits==1,...
                 'The generator cache records exact native-input reuse.');
stock=p;
stock.cce.backend='stock';
[Hs,~,~,stock_aux,~]=cce_generators(cluster,stock,[]);
result=test_true(result,'Native stock backend selected',strcmp(stock_aux.cce_aux.backend,'stock'),...
                 'The native backend preference reaches generator construction.');
result=test_close(result,'Direct and stock Hamiltonians',H,Hs,2e-10,2e-13,...
                  'Both construction paths implement the same retained Hamiltonian.');
result=test_true(result,'Retained metadata order',isequal(aux.cce_aux.physical_ids,[1 3 2]),...
                 'Observed-first physical IDs are preserved.');

% A changed representation is a cache miss, not a reused Hilbert matrix
changed=cluster;
changed.formalism='zeeman-liouv';
[changed_H,~,~,~,changed_cache]=cce_generators(changed,p,cache);
result=test_true(result,'Representation cache invalidation',...
                 changed_cache.misses==2&&isequal(size(changed_H),[cluster.dimension^2 cluster.dimension^2]),...
                 'An exact cache key includes the retained representation.');

% Compare conventional evolution against independent dense matrix exponentials
[rho,coil,~,~]=cce_observation(cluster,p);
time=0.013;
U=expm(-1i*H*time);
reference=U*rho*U';
actual=cce_evolution(cluster,H,coil,rho,time,1,'final');
result=test_close(result,'Hilbert dense exponential',actual,reference,2e-13,2e-13,...
                  'The retained Hilbert service preserves density-matrix propagation.');
trajectory=cce_evolution(cluster,H,coil,rho,time,2,'trajectory');
signal=cce_evolution(cluster,H,coil,rho,time,2,'observable');
result=test_close(result,'Initial trajectory state',trajectory(:,:,1),rho,0,0,...
                  'Trajectory output includes the unmodified initial state.');
result=test_close(result,'Non-Hermitian detection',signal(2),trace(coil'*reference),2e-13,2e-13,...
                  'Detection is conjugate-linear in the coil.');
control=2*pi*17*cce_operator(cluster,'Lx',1);
actual=cce_evolution(cluster,H+control,[],rho,time,1,'final');
U=expm(-1i*(H+control)*time);
result=test_close(result,'Finite coherent pulse',actual,U*rho*U',2e-12,2e-13,...
                  'Finite pulses add their angular-frequency control Hamiltonian.');

% Compare column-major Liouville propagation and stock evolution explicitly
liouv=cluster;
liouv.formalism='zeeman-liouv';
[L,R,K,~,~]=cce_generators(liouv,p,[]);
actual=cce_evolution(liouv,L+1i*R+1i*K,[],rho(:),time,1,'final');
result=test_close(result,'Hilbert and Liouville identity',actual,reference(:),2e-12,2e-13,...
                  'Liouville vectorisation is column-major.');
small=create(struct('isotopes',{{'1H'}},'magnet',0,'output','hush','disable',{{'hygiene'}}),struct());
bas.formalism='zeeman-liouv';
bas.approximation={'none'};
small=basis(small,bas);
small=assume(small,'labframe');
isolated=p;
isolated.cce.bath_spins=[];
isolated.cce.order=0;
one=cce_compile(small,isolated,{'labframe'});
[one,~]=cce_subsystem(one,[],struct(),[]);
[ho,ro,ko,~,~]=cce_generators(one,p,[]);
initial=[0.7;0.1i;-0.1i;0.3];
detector=[0;0;1;0];
ours=cce_evolution(one,ho+1i*ro+1i*ko,detector,initial,time,2,'observable');
stock_signal=evolution(small,ho+1i*ro+1i*ko,detector,initial,time,2,'observable');
result=test_close(result,'Stock Liouville observable semantics',ours,stock_signal,1e-13,1e-13,...
                  'Supported observable trajectories agree with stock evolution.');

% Conditional alpha/beta generators preserve complex branch phase and swapping
projected=p;
projected.cce.approximation='projected';
[branches,~,~,~,~]=cce_generators(cluster,projected,[]);
[block,detector,~,~]=cce_observation(cluster,projected);
actual=cce_evolution(cluster,branches,detector,block,time,1,'final');
db=cluster.bath_dimension;
expected=expm(-1i*branches.left*time)*reshape(block,db,db)*expm(1i*branches.right*time);
result=test_close(result,'Projected left and right branches',actual,expected(:),2e-12,2e-13,...
                  'The left and right conditional generators retain the detected phase.');
swapped=branches;
swapped.left=branches.right;
swapped.right=branches.left;
actual=cce_evolution(cluster,swapped,detector,actual,time,1,'final');
expected=expm(-1i*branches.right*time)*expected*expm(1i*branches.left*time);
result=test_close(result,'Projected pulse branch swap',actual,expected(:),2e-12,2e-13,...
                  'A refocusing pulse swaps conditional branches without conjugating the state.');

% The raising operator pumps the lower magnetic sublevel with rate gamma
rate=3.7;
dissipative=isolated;
dissipative.cce.bath_spins=2;
dissipative.cce.order=1;
system=create(struct('isotopes',{{'1H','1H'}},'magnet',0,'output','hush','disable',{{'hygiene'}}),struct());
system.bas.formalism='zeeman-liouv';
bath_model=cce_compile(system,dissipative,{'labframe'});
[bath_cluster,~]=cce_subsystem(bath_model,2,struct(),[]);
dissipative.bath_jumps=struct('support',2,'operator','L+','rate',rate);
[hd,rd,kd,~,~]=cce_generators(bath_cluster,dissipative,[]);
initial=kron([1 0;0 0],[0 0;0 1]);
initial=initial(:);
final=cce_evolution(bath_cluster,hd+1i*rd+1i*kd,[],initial,time,1,'final');
expected=kron([1 0;0 0],diag([1-exp(-rate*time),exp(-rate*time)]));
result=test_close(result,'Raising-operator pumping rate',final,expected(:),2e-14,2e-13,...
                  'L=H+iR gives decay, not exponential population growth.');
final=reshape(final,4,4);
result=test_close(result,'Dissipative trace',trace(final),1,2e-14,0,...
                  'A complete GKSL generator preserves trace.');
result=test_close(result,'Dissipative Hermiticity',final,final',2e-14,0,...
                  'A complete GKSL generator preserves Hermiticity.');

% Compare projected dissipation against a full-space coherence block
projected.bath_jumps=struct('support',2,'operator','L+','rate',rate);
[branches,~,~,~,~]=cce_generators(cluster,projected,[]);
actual=cce_evolution(cluster,branches,detector,block,time,1,'final');
fullp=projected;
fullp.cce.approximation='generalized';
[hf,rf,kf,~,~]=cce_generators(liouv,fullp,[]);
fullstate=kron([0 1;0 0],reshape(block,db,db));
fullstate=cce_evolution(liouv,hf+1i*rf+1i*kf,[],fullstate(:),time,1,'final');
fullstate=reshape(fullstate,2*db,2*db);
expected=fullstate(1:db,db+1:end);
result=test_close(result,'Projected and full dissipative block',actual,expected(:),2e-12,2e-13,...
                  'Projected jumps preserve the complete left/right GKSL contraction.');

% The projected swap matrix matches independently reversed GKSL branches
swap=branches.swap_matrix;
reversed=branches;
reversed.left=branches.right;
reversed.right=branches.left;
reversed.dissipator=swap*conj(branches.dissipator)*swap;
reversed_state=cce_evolution(cluster,reversed,detector,block,time,1,'final');
expected=reshape(conj(cce_evolution(cluster,branches,detector,...
                    reshape(reshape(block,db,db)',db^2,1),time,1,'final')),db,db).';
result=test_close(result,'Projected dissipator branch reversal',reversed_state,expected(:),2e-12,2e-13,...
                  'A coherence-adjoint transformation reverses both conditional branches.');

% Ordered multi-spin jump embedding agrees with a physical-ID operator product
ordered=p;
local=pauli(2);
nitrogen=pauli(3);
ordered.bath_jumps=struct('support',[2 3],'operator',kron(local.p,nitrogen.z),'rate',rate);
[~,multi,~,~,~]=cce_generators(liouv,ordered,[]);
op=cce_operator(cluster,{'L+','Lz'},{2,3});
d=cluster.dimension;
expected=rate*(kron(conj(op),op)-0.5*kron(eye(d),op'*op)-...
               0.5*kron((op'*op).',eye(d)));
result=test_close(result,'Ordered multi-spin jump support',multi,expected,2e-13,2e-13,...
                  'Jump support order is distinct from retained tensor order.');

% Static virtual-state corrections retain the established coherent formula
second_sys.isotopes={'1H','1H','1H'};
second_sys.magnet=0;
second_sys.output='hush';
second_sys.disable={'hygiene'};
second_system=create(second_sys,struct());
second_system.inter.zeeman.matrix{1}=1000*eye(3);
second_system=cce_test_write(second_system,[1 2],diag([100 0 0]));
second_system=cce_test_write(second_system,[1 3],diag([80 0 0]));
second=p;
second.cce.approximation='projected';
second.cce.second_order=struct('enabled',true,'gap_tolerance',1,'overlap_tolerance',1e-10);
second_model=cce_compile(second_system,second,{'labframe'});
[second_cluster,~]=cce_subsystem(second_model,[2 3],struct(),[]);
[branches,~,~,~,~]=cce_generators(second_cluster,second,[]);
x=[0 1;1 0]/2;
field=100*kron(x,eye(2))+80*kron(eye(2),x);
result=test_close(result,'Static alpha virtual-state correction',branches.left,...
                  500*eye(4)+field*field/4000,2e-12,2e-13,...
                  'The susceptibility keeps the observed matrix elements and signed gaps.');
result=test_close(result,'Static beta virtual-state correction',branches.right,...
                  -500*eye(4)-field*field/4000,2e-12,2e-13,...
                  'The opposite observed manifold has the opposite energy denominator.');

% Unsupported representations and projection-changing jumps must be rejected
bad=projected;
bad.bath_jumps.support=1;
result=test_true(result,'Observed pumping rejected in projected mode',...
                 rejects(@()cce_generators(cluster,bad,[])),...
                 'Projection-changing jumps cannot silently leave the retained manifolds.');
result=test_true(result,'Hilbert dissipation rejected',...
                 rejects(@()cce_generators(cluster,dissipative,[])),...
                 'Conventional dissipation requires a Liouville representation.');
bad=projected;
bad.bath_jumps.rate=-1;
result=test_true(result,'Negative jump rate rejected',...
                 rejects(@()cce_generators(cluster,bad,[])),...
                 'A negative rate cannot enter the GKSL dissipator.');
result=test_true(result,'Unsupported output rejected',...
                 rejects(@()cce_evolution(cluster,H,[],rho,time,1,'refocus')),...
                 'Undeclared evolution output semantics must fail explicitly.');
result=test_true(result,'Untagged branch struct rejected',...
                 rejects(@()cce_evolution(cluster,struct('left',H),[],rho,time,1,'final')),...
                 'An arbitrary structure is not a projected generator.');

end

% Catch expected contract errors without accepting a silent fallback
function answer=rejects(callback)

% Record whether the expected contract error occurs
answer=false;
try
    callback();
catch
    answer=true;
end

end


