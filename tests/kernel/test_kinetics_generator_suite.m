% Tests kinetics and flow generator helpers. Syntax:
%
%                    result=test_kinetics_generator_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks linear equilibrium, single-substance reaction and flux
% generators, matched multi-substance chemical transport, and a minimal
% hydrodynamic diffusion generator against conservation invariants.
%
% ilya.kuprov@weizmann.ac.il

function result=test_kinetics_generator_suite()

% Announce the test target
fprintf('TESTING: Kinetics and flow generator functions\n');

% State the kinetics target of the test
result=new_test_result('kernel/kinetics_generator_suite',...
                       'Kinetics and flow generator functions',...
                       'kinetic generators must conserve matter and equilibrate closed systems correctly.');

% A two-state reversible Markov generator has a closed equilibrium ratio
K=[-2 1;2 -1];
c0=[3;0];
ceq=equilibrate(K,c0);
result=test_close(result,'equilibrate two-state detailed balance',ceq,[1;2],1e-14,1e-14,...
                  'at equilibrium k_21*c_1=k_12*c_2 while total concentration is conserved');
result=test_close(result,'equilibrate zero shortcut',equilibrate(K,[0;0]),[0;0],0,0,...
                  'zero initial concentration remains zero');

% Build a single-substance flux system with a local descriptor
sys.magnet=0;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0 0};
inter.chem.parts={[1 2]}; inter.chem.concs=1;
inter.chem.reactions={struct('reactants',1,'products',1,...
    'matching',[1 2;2 1],'rate',1)};
bas.formalism='sphten-liouv'; bas.approximation={'none'};
spin_system=test_spin_system(sys,inter,bas);

% A nonempty local basis with no reaction has no reactant generators
reaction.reactants=[];
reaction.products=[];
reaction.matching=zeros(0,2);
G=react_gen(spin_system,reaction);
result=test_true(result,'react_gen single-substance empty reaction',...
                 iscell(G)&&isempty(G),...
                 'the cell-valued local descriptor is accepted for an empty reaction');

% Intramolecular flux conserves spin order column by column
Kspin=kinetics(spin_system);
result=test_close(result,'kinetics flux column sums',sum(full(Kspin),1),zeros(1,size(Kspin,2)),1e-14,1e-14,...
                  'single-substance intramolecular flux has zero column sums');

% Check matched multi-substance reaction maps and generators
inter.chem.parts={1,2}; inter.chem.concs=[1 1];
inter.chem.reactions={struct('reactants',1,'products',2,...
    'matching',[1 2],'rate',1),...
    struct('reactants',2,'products',1,...
    'matching',[2 1],'rate',1)};
bas.approximation={'none','none'};
spin_system=test_spin_system(sys,inter,bas);
G=react_gen(spin_system,spin_system.chem.reactions{1});
result=test_true(result,'react_gen matched product rows',...
                 isequal(G{1},[(5:8)' (1:4)']),...
                 'all four product descriptors map to the corresponding source descriptors');
Kspin=kinetics(spin_system);
result=test_close(result,'kinetics matched exchange',Kspin,...
                  [-eye(4) eye(4);eye(4) -eye(4)],1e-14,1e-14,...
                  'two matched first-order records give the independent block exchange generator');

% A minimal two-cell diffusion mesh must produce a conservative symmetric generator
mesh.vor.ncells=2;
mesh.vor.weights=[1;1];
mesh.vor.vertices=[0 0; 0 1];
mesh.vor.cells={ [1 2], [1 2] };
mesh.idx.active=[1;2];
mesh.idx.triangles=[1 2 3];
mesh.x=[0;1;0]; mesh.y=[0;0;1];
mesh.u=[0;0;0]; mesh.v=[0;0;0];
flow_system.sys.output='hush';
flow_system.mesh=mesh;
F=flow_gen(flow_system,struct('diff',0.5));
result=test_close(result,'flow_gen minimal diffusion',F,[-0.5 0.5;0.5 -0.5],1e-14,1e-14,...
                  'two identical cells sharing a unit boundary with D=0.5 give a symmetric conservative diffusion generator');
result=test_close(result,'flow_gen column sums',sum(full(F),1),[0 0],1e-14,1e-14,...
                  'hydrodynamic flow/diffusion generator conserves total population');

end
