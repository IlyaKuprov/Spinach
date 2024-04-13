% Mims ENDOR pulse sequence on BDPA with ideal electron pulses,
% reproducing Figure 10 from 
%
%            https://doi.org/10.1007/s00723-020-01269-z
%
% Run time: hours.
%
% i.kuprov@soton.ac.uk
% annemarie.kehl@mpinat.mpg.de

function endor_mims_bdpa()

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

% Relaxation theory
inter.relaxation={'t1_t2'};
inter.r1_rates={1e3 1e4 1e4};
inter.r2_rates={1e1 1e2 1e2};
inter.equilibrium='zero';
inter.rlx_keep='diagonal';
                       
% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'E'};
parameters.electrons=1;
parameters.nuclei=[2 3];
parameters.tau=200e-9;
parameters.n_dur=50e-6; % duration of nuclear pi pulse
parameters.rf_b1_field=-pi/(parameters.n_dur*spin('1H'));
parameters.n_frq=linspace(138e6,148e6,100);
parameters.n_rnk=2;
parameters.grid='rep_2ang_400pts_sph';

% Simulation
answer=powder(spin_system,@endor_mims_ideal,parameters,'esr');

% Plotting
figure(); plot(1e-6*parameters.n_frq,abs(answer)); 
kgrid; kylabel('abs. intensity, a.u.'); xlim tight;
kxlabel('Lab frame nuclear frequency, MHz');

end

