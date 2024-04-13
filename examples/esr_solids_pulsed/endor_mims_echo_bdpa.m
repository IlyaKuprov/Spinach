% Stimulated echo stage of the Mims ENDOR pulse sequence on
% BDPA. The nuclear pulse is not applied, this is echo dia-
% gnostics stage. The echo gets sharper when g-tensor aniso-
% tropy is increased.
%
% Run time: seconds.
%
% i.kuprov@soton.ac.uk
% annemarie.kehl@mpinat.mpg.de

function endor_mims_echo_bdpa()

% Isotopes                          
sys.isotopes={'E','1H','1H'};
                          
% Magnet field
sys.magnet=3.35;

% Interactions
inter.zeeman.matrix=cell(1,3);
inter.zeeman.matrix{1}=[2.00263 0.00000 0.00000;
                        0.00000 2.00260 0.00000;
                        0.00000 0.00000 2.00257];
inter.coupling.matrix=cell(3,3);
inter.coupling.matrix{1,2}=[7.70 0.00 0.00
                            0.00 5.30 0.00
                            0.00 0.00 2.00]*1e6;
inter.coupling.matrix{1,3}=[1.00 0.00 0.00
                            0.00 1.00 0.00
                            0.00 0.00 1.26]*1e6;
                       
% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'E'};
parameters.electrons=1;
parameters.tau=200e-9;
parameters.n_dur=50e-6;
parameters.grid='rep_2ang_400pts_sph';
parameters.nsteps=200;

% Simulation
answer=powder(spin_system,@endor_mims_echo,parameters,'esr');

% Plotting
time_axis=linspace(0,2*parameters.tau,parameters.nsteps+1);
figure(); plot(1e9*time_axis,real(answer)); xlim tight;
kgrid; kxlabel('time, ns'); kylabel('intensity, a.u.');

end

