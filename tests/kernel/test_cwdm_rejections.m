% Tests rejection of cross-substance physical specifications. Syntax:
%
%                       result=test_cwdm_rejections()
%
% Outputs:
%
%     result - T3 checks of coupling, product, and symmetry rejections
%
% Direct-sum substances cannot carry pair interactions or product operators
% between blocks. Symmetry inputs use local columns and reject an index
% outside their substance. Error identifiers and messages are asserted.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_rejections()

% Announce the specification boundary
fprintf('TESTING: CWDM cross-substance rejections (T3)\n');
result=new_test_result('kernel/cwdm_rejections','CWDM rejections',...
                      'Cross-substance physics is an error, never a zero.');

% Construct independent two-spin substances
sys.magnet=1; sys.isotopes={'1H','1H','1H','1H'};
inter.chem.parts={1:2,3:4}; inter.chem.concs=[1 0];
bas.formalism='sphten-liouv'; bas.approximation={'none','none'};
s=test_spin_system(sys,inter,bas);

% Every spin must belong to one of the compiled substances
bad=inter; bad.chem.parts={1:2,3}; rejected=false;
try
    test_spin_system(sys,bad,bas);
catch err
    rejected=strcmp(err.identifier,'Spinach:basis:incompletePartition');
end
result=test_true(result,'incomplete partition',rejected,...
                 'an unassigned spin is rejected before direct-sum compilation');

% Column-vector parts are valid complete partitions too
columns=s; columns.chem.parts={[1;2],[3;4]};
by_columns=basis(columns,bas);
result=test_true(result,'column partition',isequal(by_columns.bas.basis,s.bas.basis),...
                 'partition coverage does not depend on spin-index vector orientation');

% Reject a scalar coupling between different substances
bad=inter; bad.coupling.scalar=cell(4); bad.coupling.scalar{2,3}=10;
rejected=false;
try
    test_spin_system(sys,bad,bas);
catch err
    rejected=strcmp(err.identifier,'Spinach:create:crossSubstanceCoupling')&&...
             contains(err.message,'couplings detected between spins in different chemical species');
end
result=test_true(result,'cross coupling',rejected,...
                 'create raises the named cross-substance coupling error');

% Reject products even when one explicitly requested factor is identity
for labels={{'Lz','Lz'},{'E','Lz'}}
    for constructor={@operator,@state}
        rejected=false;
        try
            constructor{1}(s,labels{1},{1,3});
        catch err
            rejected=strcmp(err.identifier,'Spinach:which_subst:crossSubstance')&&...
                     contains(err.message,'spin list crosses chemical boundaries');
        end
        result=test_true(result,['product ' func2str(constructor{1}) ' ' labels{1}{1}],...
                         rejected,'cross-substance product specifications raise a named error');
    end
end

% Reject a symmetry label that would address another substance globally
bad=bas; bad.sym_group={{'S2'},{}}; bad.sym_spins={{[1 3]},{}};
rejected=false;
try
    basis(s,bad);
catch err
    rejected=contains(err.message,'incorrect spin labels in bas.sym_spins');
end
result=test_true(result,'cross symmetry',rejected,...
                 'symmetry indices must be local to the declared substance');

% Segmented Zeeman filters reject unsupported constructions
sys.isotopes={'1H','1H'}; inter.chem.parts={1,2};
for formalism={'zeeman-hilb','zeeman-liouv'}
    bas.formalism=formalism{1}; s=test_spin_system(sys,inter,bas);
    rho=ones(s.bas.offsets(end),1);
    if strcmp(formalism{1},'zeeman-hilb'), rho=diag(rho); end
    for selector={'correlation','decouple','homospoil'}
        rejected=false;
        try
            switch selector{1}
                case 'correlation'
                    correlation(s,rho,0,1);
                case 'decouple'
                    [~,rho]=decouple(s,[],rho,1);
                case 'homospoil'
                    homospoil(s,rho,'destroy');
            end
        catch err
            rejected=strcmp(err.identifier,['Spinach:' selector{1} ':segmentedZeeman']);
        end
        result=test_true(result,[selector{1} ' ' formalism{1}],rejected,...
                         'unsupported direct-sum Zeeman filtering raises the named error');
    end
    bad=bas; bad.sym_group={{'S2'},{}}; bad.sym_spins={{1},{}};
    rejected=false;
    try
        basis(s,bad);
    catch err
        rejected=strcmp(err.identifier,'Spinach:basis:segmentedZeeman');
    end
    result=test_true(result,['Zeeman symmetry rejection ' formalism{1}],rejected,...
                     'only genuinely segmented Zeeman symmetry is rejected');
end

% Segmented wavefunctions cannot represent the compiled direct sum
sys.isotopes={'1H','1H','1H'};
inter=struct(); inter.chem.parts={1,2:3}; inter.chem.concs=[1 0];
bas=struct('formalism','zeeman-wavef','approximation',{{'none','none'}});
s=test_spin_system(sys,inter,bas); rejected=false;
try
    state(s,[0.5 0.5 0.5]);
catch err
    rejected=strcmp(err.identifier,'Spinach:state:segmentedZeeman');
end
result=test_true(result,'segmented wavefunction',rejected,...
                 'unsupported segmented wavefunctions raise a named error');
inter=struct(); bas.approximation={'none'};
s=test_spin_system(sys,inter,bas);
result=test_close(result,'single-substance wavefunction',...
                  state(s,[0.5 0.5 0.5]),sparse(1,1,1,8,1),0,0,...
                  'the existing tensor-product wavefunction is retained for one substance');
sys.isotopes={'1H','1H'};

% Single-substance Zeeman units and Boltzmann states retain stock normalisation
inter=struct('temperature',298); bas.approximation={'none'};
for formalism={'zeeman-hilb','zeeman-liouv'}
    bas.formalism=formalism{1}; s=test_spin_system(sys,inter,bas);
    beta=s.tols.hbar/(s.tols.kbol*inter.temperature);
    H=diag([-3 -1 1 3])/beta;
    expected=diag(exp(-beta*diag(H))); expected=expected/trace(expected);
    unit=speye(4);
    if strcmp(formalism{1},'zeeman-liouv')
        H=kron(speye(4),H); expected=expected(:);
        unit=unit(:)/norm(unit(:));
    end
    result=test_close(result,['single-substance unit ' formalism{1}],...
                      unit_state(s),unit,0,0,'single-substance unit normalisation is unchanged');
    result=test_close(result,['single-substance equilibrium ' formalism{1}],...
                      equilibrium(s,H),expected,1e-10,1e-10,...
                      'the state equals the trace-normalised Boltzmann exponential');
end

% Coherent states cannot use a tensor product in a segmented Zeeman space
sys.isotopes={'1H','C3'}; inter=struct();
inter.chem.parts={1,2}; inter.chem.concs=[1 1]; bas.approximation={'none','none'};
for formalism={'zeeman-hilb','zeeman-liouv'}
    bas.formalism=formalism{1}; s=test_spin_system(sys,inter,bas); rejected=false;
    try
        coherent(s,2,.5);
    catch err
        rejected=strcmp(err.identifier,'Spinach:coherent:segmentedZeeman');
    end
    result=test_true(result,['segmented coherent ' formalism{1}],rejected,...
                     'the coherent tensor product is rejected before construction');
    fprintf('CWDM_COHERENT %s named_rejection=%d\n',formalism{1},rejected);
    local_bas=bas; local_bas.approximation={'none'};
    local=test_spin_system(sys,struct(),local_bas);
    amps=[1 .5 .5^2/sqrt(2)]; amps=amps/norm(amps);
    expected=kron(eye(2),amps'*amps);
    if strcmp(formalism{1},'zeeman-liouv'), expected=expected(:); end
    result=test_close(result,['single coherent ' formalism{1}],coherent(local,2,.5),...
                      expected,1e-14,1e-14,'the supported truncated coherent product is unchanged');
end

% Caller-supplied generators must not transfer between compiled substances
sys=struct('magnet',1,'isotopes',{{'1H','1H'}}); inter=struct();
inter.chem.parts={1,2}; inter.chem.concs=[1 1];
bas=struct('formalism','sphten-liouv','approximation',{{'none','none'}});
s=test_spin_system(sys,inter,bas); rho=state(s,'Lz',1)+state(s,'Lz',2);
for n=1:2
    cross=sparse(3,7,1,8,8); if n==2, cross=cross'; end
    for caller={'reduce','evolution'}
        rejected=false;
        try
            if strcmp(caller{1},'reduce')
                reduce(s,cross,rho);
            else
                evolution(s,cross,[],rho,0.1,1,'final');
            end
        catch err
            rejected=strcmp(err.identifier,'Spinach:reduce:crossSubstanceGenerator');
        end
        result=test_true(result,['cross generator ' caller{1} ' ' int2str(n)],rejected,...
                         'cross-substance input is rejected before projector construction');
        fprintf('CWDM_REDUCE_CROSS caller=%s direction=%d named_rejection=%d\n',caller{1},n,rejected);
    end
end

% Independent blocks retain full-generator dynamics without numerical rounding
s.sys.disable=[s.sys.disable {'clean-up'}];
L=operator(s,'Lx',1)+2*operator(s,'Lx',2);
actual=evolution(s,L,[],rho,0.1,1,'final'); expected=expm(-0.1i*full(L))*rho;
result=test_close(result,'independent generator evolution',actual,expected,1e-12,1e-12,...
                  'per-substance reduction preserves block-diagonal generator action');
fprintf('CWDM_REDUCE_BLOCKS state_error=%.16g\n',norm(actual-expected));

% Validate the Hamiltonian action independently in both equilibrium blocks
sys=struct('magnet',1,'isotopes',{{'1H','1H'}}); inter=struct();
inter.chem.parts={1,2}; inter.chem.concs=[1 1]; inter.temperature=298;
bas.formalism='sphten-liouv'; bas.approximation={'none','none'};
s=test_spin_system(sys,inter,bas);
left=operator(s,'Lz',1,'left')+operator(s,'Lz',2,'left');
for n=1:2
    idx=(s.bas.offsets(n)+1):s.bas.offsets(n+1);
    mixed=left; comm=operator(s,'Lz',n); mixed(idx,idx)=comm(idx,idx);
    rejected=false;
    try
        equilibrium(s,mixed);
    catch err
        rejected=strcmp(err.identifier,'Spinach:equilibrium:notLeftProduct')&&...
                 contains(err.message,['substance ' int2str(n)]);
    end
    result=test_true(result,['mixed equilibrium block ' int2str(n)],rejected,...
                     'a valid left product cannot hide a commutator in another substance');
    fprintf('CWDM_EQUILIBRIUM_MIXED block=%d named_rejection=%d\n',n,rejected);
end

% Proper left products retain their independently normalised Boltzmann states
rho=equilibrium(s,left); unit=unit_state(s);
beta=s.tols.hbar/(s.tols.kbol*s.rlx.temperature);
for n=1:2
    idx=(s.bas.offsets(n)+1):s.bas.offsets(n+1);
    expected=expm(-beta*full(left(idx,idx)))*unit(idx);
    expected=expected/dot(unit(idx),expected);
    result=test_close(result,['valid equilibrium block ' int2str(n)],rho(idx),expected,...
                      1e-14,1e-14,'each valid block retains its Boltzmann state');
end

% Exactly zero spinful Hamiltonian blocks retain their unit state
for n=1:2
    idx=(s.bas.offsets(n)+1):s.bas.offsets(n+1);
    zero=left; zero(idx,idx)=0; actual=equilibrium(s,zero);
    expected=rho; expected(idx)=unit(idx);
    result=test_close(result,['zero equilibrium block ' int2str(n)],actual,expected,...
                      1e-14,1e-14,'a zero block is a valid left product, not a nonzero commutator');
    fprintf('CWDM_EQUILIBRIUM_ZERO block=%d unit_error=%.16g\n',n,norm(actual(idx)-unit(idx)));
end

% Reject cross-substance terms in either direction, including anisotropic input
for n=1:2
    cross=sparse(3,5,1,s.bas.offsets(end),s.bas.offsets(end));
    if n==2, cross=cross'; end
    for anisotropic=[false true]
        rejected=false;
        try
            if anisotropic
                Q={cell(3)}; Q{1}{2,2}=cross;
                equilibrium(s,left,Q,[0 0 0]);
            else
                equilibrium(s,left+cross);
            end
        catch err
            rejected=strcmp(err.identifier,'Spinach:equilibrium:crossSubstanceHamiltonian');
        end
        result=test_true(result,['cross equilibrium ' int2str(n) ' ' int2str(anisotropic)],...
                         rejected,'the assembled Hamiltonian must not transfer between substances');
        fprintf('CWDM_EQUILIBRIUM_CROSS direction=%d anisotropic=%d named_rejection=%d\n',...
                n,anisotropic,rejected);
    end
end

end
