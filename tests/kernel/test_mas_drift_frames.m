% Test: kernel/mas_drift_frames | MAS drift frames
% Tests numerical frame ordering and phase indexing of MAS drift stacks.
% Syntax: result=test_mas_drift_frames()
%
% Talos <talos@spindynamics.org>

function result=test_mas_drift_frames()

% Build a coupled heteronuclear system with nonzero shifts
sys.isotopes={'15N','13C'};
sys.magnet=9.4;
inter.zeeman.scalar={12,35};
inter.coupling.matrix={[],[100 0 0;0 100 0;0 0 -200];[],[]};
bas.formalism='zeeman-hilb';
bas.approximation={'none'};
spin_system=assume(test_spin_system(sys,inter,bas),'labframe');
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.n_ticks=4;
parameters.n_phases=2;
parameters.n_slices=5;
parameters.grid=struct('alphas',[0 0.2],'betas',[0.1 0.8],...
                       'gammas',[0.2 0.3],'weights',[1 3]/4);
parameters.rframes={{'15N',2},{'13C',1}};
[drifts,weights]=mas_drifts(spin_system,parameters);
result=new_test_result('kernel/mas_drift_frames','MAS drift frames','Checks ordered numerical frames.');
result=test_close(result,'case weights',weights,[1;1;3;3]/8,0,0,'Weights follow orientation, then phase.');

% Independently sample orientation two, phase two, at its first midpoint
[H,Q]=hamiltonian(spin_system);
D=euler2dcm(0,acos(parameters.axis(3)),0)*...
  euler2dcm(0,0,5*pi/4)*euler2dcm(0.2,0.8,0.3);
[alpha,beta,gamma]=dcm2euler(D);
H=H+orientation(Q,[alpha beta gamma]);
Cn=carrier(spin_system,'15N');
Cc=carrier(spin_system,'13C');
H=rotframe(spin_system,Cn,(H+H')/2-Cc,'15N',2);
H=rotframe(spin_system,Cc,(H+H')/2+Cc,'13C',1);
result=test_close(result,'ordered frames',drifts{4}{1},full(H),1e-7,1e-9,...
                  'Sample the nonzero phase at the first tick midpoint.');
result=test_close(result,'cycle wrap',drifts{4}{5},drifts{4}{1},0,0,'A full revolution repeats the tick.');

% Compare legacy and explicit single-frame calls on a uniform grid
parameters.grid.weights=[0.5 0.5];
parameters=rmfield(parameters,'rframes');
parameters.spins={'15N'};
sys.isotopes={'15N'};
single=assume(test_spin_system(sys,struct(),bas),'labframe');
legacy=mas_drifts(single,parameters);
parameters.rframes={{'15N',2}};
explicit=mas_drifts(single,parameters);
result=test_close(result,'legacy frame',legacy{1}{2},explicit{1}{2},0,0,'Legacy order is two.');

end


