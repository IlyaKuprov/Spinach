% Tests kinetics and flow generator helpers. Syntax:
%
%                    result=test_kinetics_generator_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks linear equilibrium, chemical reaction generators, full
% chemical kinetics generators, and a minimal hydrodynamic diffusion
% generator against conservation and detailed-balance invariants.
%
% ilya.kuprov@weizmann.ac.il

function result=test_kinetics_generator_suite()

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

% Build a two-site exchange Spinach system used by react_gen and kinetics
sys.magnet=0;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0 0};
inter.chem.parts={1,2};
inter.chem.rates=[-1 1;1 -1];
inter.chem.concs=[1 1];
bas.formalism='sphten-liouv'; bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% react_gen must build a conservative drain/fill mapping for a specified reaction
reaction.reactants=1;
reaction.products=2;
reaction.matching=[1 2];
G=react_gen(spin_system,reaction);
result=test_close(result,'react_gen one reaction count',numel(G),1,0,0,...
                  'one reactant channel produces one reaction-generator matrix');
result=test_close(result,'react_gen conservation',sum(full(G{1}),1),zeros(1,size(G{1},2)),1e-14,1e-14,...
                  'reaction drain and fill terms conserve total population column by column');

% Full kinetics generator for symmetric exchange must conserve population
Kspin=kinetics(spin_system);
result=test_close(result,'kinetics column sums',sum(full(Kspin),1),zeros(1,size(Kspin,2)),1e-14,1e-14,...
                  'closed chemical kinetics generator has zero column sums');

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
result=test_close(result,'flow_gen minimal diffusion',F,[-1 1;1 -1],1e-14,1e-14,...
                  'two identical cells sharing a unit boundary with D=0.5 give a symmetric conservative diffusion generator');
result=test_close(result,'flow_gen column sums',sum(full(F),1),[0 0],1e-14,1e-14,...
                  'hydrodynamic flow/diffusion generator conserves total population');

end
