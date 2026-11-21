% Tests input merging against independently specified systems. Syntax:
%
%                         result=test_merge_inp()
%
% Outputs:
%
%     result - input-merging regression test results
%
% Chemical substances are independent direct-sum blocks. Their input
% concentrations and correlation times retain substance order regardless
% of vector orientation. Bosonic mode arrays retain particle order, with
% no interactions between different input parts. The handwritten inputs
% below therefore define the same Hamiltonian and relaxation operators.
% Operator comparisons allow only accumulated floating-point roundoff.
%
% ilya.kuprov@weizmann.ac.il

function result=test_merge_inp()

% Announce the merge contract
fprintf('TESTING: input merging\n');
result=new_test_result('kernel/merge_inp','Input merging',...
                      'Merged inputs equal independently assembled systems.');

% Set common physical and execution settings
sys.magnet=1; sys.isotopes={'1H'};
sys.output='hush'; sys.disable={'hygiene'}; sys.parallel={'local',1};
a.chem.parts={1}; a.chem.concs=2;
bas.formalism='sphten-liouv';

% Check bosonic modes, default concentrations, and column-shaped inputs
for n=1:7
    switch n
        case 1
            % D1: a mode subsystem and a spin-only subsystem
            left=sys; left.magnet=0; left.isotopes={'1H','C3'};
            right=sys; right.magnet=0;
            inter.modes.frqs={[],1000}; inter.modes.carriers={[],1e9};
            inter.modes.anharms={[],5}; inter.modes.lifetimes={[],0.5};
            inter.modes.t2_times={[],0.05};
            inter.modes.longitudinal=cell(2);
            inter.modes.longitudinal{1,2}=3;
            other.temperature=0; inter.temperature=0; expected.temperature=0;
            inter.chem.parts={1:2}; other.chem.parts={1};
            expected_sys=left; expected_sys.isotopes={'1H','C3','1H'};
            expected.modes.frqs={[],1000,[]};
            expected.modes.carriers={[],1e9,[]};
            expected.modes.anharms={[],5,[]};
            expected.modes.lifetimes={[],0.5,[]};
            expected.modes.t2_times={[],0.05,[]};
            expected.modes.longitudinal=cell(3);
            expected.modes.longitudinal{1,2}=3;
            expected.chem.parts={1:2,3}; expected.chem.concs=[1 1];
            bas.formalism='zeeman-liouv'; bas.approximation={'none','none'};
        case 2
            % D2: create() defaults each single-substance concentration to one
            left=sys; right=sys;
            inter.chem.parts={1}; other.chem.parts={1};
            expected_sys=sys; expected_sys.isotopes={'1H','1H'};
            expected.chem.parts={1,2}; expected.chem.concs=[1 1];
            bas.formalism='sphten-liouv'; bas.approximation={'none','none'};
        case {3,4}
            % D3/D4: two substances followed by one substance
            left=sys; left.isotopes={'1H','13C'}; right=sys;
            inter.chem.parts={1;2}; inter.chem.concs=[2;3]; other=a;
            expected_sys=sys; expected_sys.isotopes={'1H','13C','1H'};
            expected.chem.parts={1,2,3}; expected.chem.concs=[2 3 2];
            bas.approximation={'none','none','none'};
            if n==4
                inter.tau_c={1e-9;2e-9}; other.tau_c={3e-9};
                inter.relaxation={'redfield'}; other.relaxation={'redfield'};
                inter.rlx_keep='diagonal'; other.rlx_keep='diagonal';
                inter.equilibrium='zero'; other.equilibrium='zero';
                expected.tau_c={1e-9,2e-9,3e-9};
                expected.relaxation={'redfield'}; expected.equilibrium='zero';
                expected.rlx_keep='diagonal';
            end
        case 5
            % A single mode input must also be an identity merge
            left=sys; left.magnet=0; left.isotopes={'1H','C3'};
            inter.modes.frqs={[],1000}; inter.modes.anharms={[],5};
            expected_sys=left; expected=inter;
            bas.formalism='zeeman-liouv'; bas.approximation={'none'};
        case 6
            % Independent mode networks with spin-indexed tensor derivatives
            left=sys; left.magnet=0; left.isotopes={'1H','C3','C3'};
            right=sys; right.magnet=0; right.isotopes={'13C','C3'};
            inter.modes.frqs={[],1000,1200};
            inter.modes.carriers={[],1e9,2e9}; inter.modes.anharms={[],5,6};
            inter.modes.lifetimes={[],0.5,[]};
            inter.modes.qfactors={[],[],1e12}; inter.modes.t2_times={[],0.05,[]};
            other.modes.frqs={[],1400}; other.modes.anharms={[],7};
            other.modes.linewidths={[],1};
            expected_sys=left; expected_sys.isotopes={'1H','C3','C3','13C','C3'};
            expected.modes.frqs={[],1000,1200,[],1400};
            expected.modes.carriers={[],1e9,2e9,[],[]};
            expected.modes.anharms={[],5,6,[],7};
            expected.modes.lifetimes={[],0.5,[],[],[]};
            expected.modes.qfactors={[],[],1e12,[],[]};
            expected.modes.t2_times={[],0.05,[],[],[]};
            expected.modes.linewidths={[],[],[],[],1};
            inter.modes.exchange=cell(3); inter.modes.exchange{2,3}=8;
            expected.modes.exchange=cell(5); expected.modes.exchange{2,3}=8;
            inter.modes.kerr=cell(3); inter.modes.kerr{2,3}=9;
            expected.modes.kerr=cell(5); expected.modes.kerr{2,3}=9;
            inter.modes.dispersive=cell(3); inter.modes.dispersive{1,2}=10;
            expected.modes.dispersive=cell(5); expected.modes.dispersive{1,2}=10;
            inter.modes.coupling_mod=cell(3);
            inter.modes.coupling_mod{2,2}={cell(3)};
            inter.modes.coupling_mod{2,2}{1}{1,1}=diag([1 2 3]);
            other.modes.coupling_mod=cell(2);
            other.modes.coupling_mod{2,2}={cell(2)};
            other.modes.coupling_mod{2,2}{1}{1,1}=diag([4 5 6]);
            expected.modes.coupling_mod=cell(5);
            expected.modes.coupling_mod{2,2}={cell(5)};
            expected.modes.coupling_mod{2,2}{1}{1,1}=diag([1 2 3]);
            expected.modes.coupling_mod{5,5}={cell(5)};
            expected.modes.coupling_mod{5,5}{1}{4,4}=diag([4 5 6]);
            inter.modes.zeeman_mod=cell(3);
            inter.modes.zeeman_mod{2,2}={{[1 2 3],[],[]}};
            other.modes.zeeman_mod=cell(2);
            other.modes.zeeman_mod{2,2}={{[4 5 6],[]}};
            expected.modes.zeeman_mod=cell(5);
            expected.modes.zeeman_mod{2,2}={{[1 2 3],[],[],[],[]}};
            expected.modes.zeeman_mod{5,5}={{[],[],[],[4 5 6],[]}};
            other.temperature=0; inter.temperature=0; expected.temperature=0;
            inter.chem.parts={1:3}; other.chem.parts={1:2};
            expected.chem.parts={1:3,4:5}; expected.chem.concs=[1 1];
            bas.formalism='zeeman-liouv'; bas.approximation={'none','none'};
        case 7
            % Explicit concentrations imply the local default partitions
            left=sys; right=sys;
            inter.chem.concs=2; other.chem.concs=3;
            expected_sys=sys; expected_sys.isotopes={'1H','1H'};
            expected.chem.parts={1,2}; expected.chem.concs=[2 3];
            bas.formalism='sphten-liouv'; bas.approximation={'none','none'};
    end

    % Validate each part independently before merging
    create(left,inter);
    if n==5
        [merged_sys,merged_inter]=merge_inp({left},{inter});
    else
        create(right,other);
        [merged_sys,merged_inter]=merge_inp({left,right},{inter,other});
    end
    observed=assume(basis(create(merged_sys,merged_inter),bas),'labframe');
    reference=assume(basis(create(expected_sys,expected),bas),'labframe');
    label=['case ' int2str(n)];

    % Compare the physical records and independent operator constructions
    result=test_true(result,[label ' particles'],...
                     isequal(observed.comp.isotopes,reference.comp.isotopes),...
                     'Particle order is the concatenation of the input parts.');
    result=test_true(result,[label ' chemistry'],...
                     isequaln(observed.chem,reference.chem),...
                     'Substance indices and absolute concentrations match the handwritten input.');
    result=test_true(result,[label ' interactions'],...
                     isequaln(observed.inter,reference.inter),...
                     'Interaction tensors, mode parameters, and index offsets match.');
    result=test_true(result,[label ' relaxation records'],...
                     isequaln(observed.rlx,reference.rlx),...
                     'Per-substance relaxation parameters retain their physical meaning.');
    tolerance=100*eps()*sum(reference.bas.nstates);
    result=test_close(result,[label ' Hamiltonian'],...
                      hamiltonian(observed),hamiltonian(reference),...
                      tolerance,tolerance,'Both inputs describe the same independent substances.');
    result=test_close(result,[label ' relaxation'],...
                      relaxation(observed),relaxation(reference),...
                      tolerance,tolerance,'Both inputs describe the same relaxation processes.');
    clear inter other expected
end

% A single global magnetic field cannot represent conflicting input fields
right=sys; right.magnet=2; rejected=false;
try
    merge_inp({sys,right},{a,a});
catch err
    rejected=contains(err.message,'magnet')&&contains(err.message,'same');
end
result=test_true(result,'conflicting magnetic fields',rejected,...
                 'An impossible global merge must identify the conflicting field.');

end


