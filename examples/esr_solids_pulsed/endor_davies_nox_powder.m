% Davies ENDOR simulation for a nitroxide radical. Soft pulses
% are simulated using Fokker-Planck formalism. This is a pain-
% fully slow brute-force time-domain simulation with explicit
% soft pulses and full account of the effect of the orientati-
% on selection using a large spherical averaging grid.
%
% Calculation time: hours.
%
% ilya.kuprov@weizmann.ac.il

function endor_davies_nox_powder()

% Isotopes                          
sys.isotopes={'E','14N'};
                          
% Magnet field
sys.magnet=3.5;

% Interactions
inter.zeeman.matrix=cell(1,2);
inter.zeeman.matrix{1}=[2.01045  0.00000  0.00000
                        0.00000  2.00641  0.00000
                        0.00000  0.00000  2.00211];
inter.coupling.matrix=cell(2,2);
inter.coupling.matrix{1,2}=[1.2356  0.0000  0.6322
                            0.0000  1.1266  0.0000
                            0.6322  0.0000  8.2230]*1e7;
                       
% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory
inter.relaxation={'t1_t2'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.r1_rates={20e6 0.5e6};
inter.r2_rates={20e6 0.5e6};

% Disable trajectory-level SSR algorithms
sys.disable={'trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'E','14N'};
parameters.rho0=state(spin_system,'Lz','E');
parameters.coil=state(spin_system,'L+','E');
parameters.grid='rep_2ang_12800pts_sph';
parameters.offset=[-2e8 0];
parameters.method='expm';
parameters.verbose=0;

% Electron pulse parameters
parameters.e_rnk=2;
parameters.e_phi=pi/2;
parameters.e_frq=-300e6;
parameters.e_dur=10e-9;
parameters.e_pwr=2*pi*16.5e7;

% Nucleus pulse parameters
parameters.n_rnk=3;
parameters.n_phi=pi/2;
parameters.n_frq=linspace(-80e6,80e6,200);
parameters.n_dur=1e-7;
parameters.n_pwr=pi*1e7;

% Simulation
answer=powder(spin_system,@endor_davies,parameters,'esr');

% Plotting
figure(); plot(parameters.n_frq/1e6,real(answer));
kgrid; axis tight; kylabel('(RF on)/(RF off)'); 
kxlabel('Nuclear frequency, MHz');

end

