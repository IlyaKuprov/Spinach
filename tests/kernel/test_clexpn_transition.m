% Native observed-transition selection, units, ordering, and degeneracy tests.
%
%                         result=test_clexpn_transition()
%
% Outputs:
%
%     result - regression result with recorded checks and failures
%
% ilya.kuprov@weizmann.ac.il

function result=test_clexpn_transition()

% State the descriptive-transition contract independently of the old driver
result=new_test_result('kernel/clexpn_transition','Observed transition selection',...
                       'Magnetic quantum numbers, energy levels, explicit states and gaps.');
sys=struct('isotopes',{{'14N','1H'}},'magnet',1,'output','hush',...
           'disable',{{'hygiene'}},'parallel',{{'local',1}});
system=create(sys,struct());
p=struct('observed',1,'transition',struct('spins',1,'ms',[0 1]),...
         'orientation',[0 0 0]);
[transition,meta]=cce_transition(system,p,'labframe');
result=test_true(result,'Spin-one descending-m selection',...
                 isequal(transition.alpha,[0;1;0])&&isequal(transition.beta,[1;0;0]),...
                 'The m=0,+1 selection reproduces the frozen one-nucleus Hahn manifolds.');
result=test_close(result,'Native angular-frequency gap',meta.gap,...
                    -spin('14N')*sys.magnet,1e-8,1e-14,'Analytic native angular-frequency gap.');

% Preserve caller order for unlike observed spins
p.observed=[2 1];
p.transition=struct('spins',[2 1],'ms',[0.5 -0.5;0 1]);
[transition,meta]=cce_transition(system,p,'labframe');
result=test_true(result,'Ordered multi-spin magnetic states',...
                 isequal(transition.alpha,kron([1;0],[0;1;0]))&&...
                 isequal(transition.beta,kron([0;1],[1;0;0]))&&...
                 isequal(meta.multiplicities,[2 3]),...
                 'Physical spin IDs determine tensor-product ordering; sorting is not imposed.');
result=test_close(result,'Ordered two-spin energy gap',meta.gap,...
                    (spin('1H')-spin('14N'))*sys.magnet,1e-7,1e-14,'Analytic native angular-frequency gap.');

% Energy levels select ascending energies with deterministic eigenvector phase
p.observed=2;
p.transition=struct('spins',2,'levels',[1 2]);
[a,meta]=cce_transition(system,p,'labframe');
[b,~]=cce_transition(system,p,'labframe');
result=test_true(result,'Energy order and repeatable phase',...
                 isequal(a,b)&&isequal(a.alpha,[1;0])&&isequal(a.beta,[0;1]),...
                 'A proton in a positive field has an unambiguous ordered transition.');
result=test_close(result,'Eigenstate angular gap',meta.gap,spin('1H')*sys.magnet,...
                    1e-7,1e-14,'Analytic native angular-frequency gap.');

% Reject unresolved zero-field energy indices instead of choosing an axis
sys.magnet=0;
zero=create(sys,struct());
rejected=false;
try
    cce_transition(zero,p,'labframe');
catch exception
    rejected=strcmp(exception.identifier,'cce_transition:degeneracy');
end
result=test_true(result,'Unresolved degenerate states refused',rejected,...
                 'Energy indices do not define a transition within a degenerate manifold.');
p.transition=struct('spins',2,'alpha',[1;1]/sqrt(2),'beta',[1;-1]/sqrt(2));
[a,meta]=cce_transition(zero,p,'labframe');
result=test_true(result,'Explicit states resolve the physical choice',...
                 isequal(a.alpha,p.transition.alpha)&&isequal(a.beta,p.transition.beta)&&meta.gap==0,...
                 'Explicit orthonormal states are a physical selector, not a legacy translator.');

% Validate selector ownership, allowed quantum numbers, shapes, and populations
invalid=cell(1,8);
invalid{1}=p;
invalid{1}.transition.extra=true;
invalid{2}=p;
invalid{2}.transition.spins=1;
invalid{3}=p;
invalid{3}.transition.alpha=[2;0];
invalid{4}=p;
invalid{4}.transition.ms=[0.5 -0.5];
invalid{5}=p;
invalid{5}.transition=struct('spins',2,'ms',[1 0]);
invalid{6}=p;
invalid{6}.transition=struct('spins',2,'levels',[1 1]);
invalid{7}=p;
invalid{7}.transition=struct('spins',2,'ms',[0.5; -0.5]);
invalid{8}=p;
invalid{8}.orientation=[0 0 NaN];
for k=1:numel(invalid)
    rejected=false;
    try
        cce_transition(system,invalid{k},'labframe');
    catch exception
        rejected=startsWith(exception.identifier,'cce_transition:');
    end
    result=test_true(result,sprintf('Transition rejection %d',k),rejected,...
                     'Invalid or ambiguous state choices are never normalized or guessed.');
end

% Check actual orientation with an explicitly prescribed total Zeeman tensor
spec=struct('mode','explicit','isotopes',{{'13C'}},'coordinates_A',[0 0 0],...
            'tensors_Hz',diag([100 200 300]));
[native,inter]=spin_bath(struct('magnet',1),struct(),spec);
native.output='hush';
native.disable={'hygiene'};
native.parallel={'local',1};
anisotropic=create(native,inter);
selection=struct('observed',1,'transition',struct('spins',1,'levels',[1 2]),...
                 'orientation',[0 pi/4 0]);
[rotated,meta]=cce_transition(anisotropic,selection,'labframe');
result=test_close(result,'Frame-dependent Zeeman gap',meta.gap,2*pi*sqrt(100^2+200^2),...
                  1e-6,1e-12,'The rotated tensor row gives the effective-field magnitude.');
result=test_true(result,'Rotated eigenstates remain orthonormal',...
                 norm([rotated.alpha rotated.beta]'*[rotated.alpha rotated.beta]-eye(2),'fro')<1e-14&&...
                 all(meta.residuals<1e-6),'The reported energy states use the requested frame.');

end


