% NOESY/EXSY experiment on phenylaziridine, including scalar relaxation
% of the second kind induced by the 14N nucleus, in a situation where
% the chemical exchange is intermediate, lines are broadened and scalar
% relaxation of the first kind must be accounted for. Set to reproduce
% Figures 1a and 4a from 
%
%               http://dx.doi.org/10.1002/ange.201410271
%
% All parameters, except for the isotropic chemical shifts, exchange ra-
% tes and correlation times come from a DFT calculation.
%
% Calculation time: minutes on a Tesla V100 GPU, much longer on CPU
%
% ilya.kuprov@weizmann.ac.il
% barbara.odell@chem.ox.ac.uk
% tim.claridge@chem.ox.ac.uk

function aziridine_exsy_2() 

% Magnet induction 
sys.magnet=11.75;

% Isotopes 
sys.isotopes={'1H','1H','1H','14N','1H','1H','1H','1H','1H','1H',...
              '1H','1H','1H','14N','1H','1H','1H','1H','1H','1H'};  

% Coordinates (Angstrom)  
inter.coordinates={[-3.41285300    0.73156000    1.07036600];
                   [-2.18218300   -0.53340100    1.62534100];
                   [-1.83251000    1.55844200   -0.65127900];
                   [-2.51263500   -0.53063900   -0.50722600];
                   [-3.28847400   -0.19858900   -1.08239500];
                   [ 0.45545168    2.23095146   -0.76407436];
                   [-0.46596257   -1.76861142    0.49251865];
                   [ 2.88380973    1.74814805   -0.56768603];
                   [ 1.97436454   -2.26233655    0.71081156];
                   [ 3.65174100   -0.50071826    0.17866774];
                   [-3.41285300    0.73156000    1.07036600];
                   [-2.18218300   -0.53340100    1.62534100];
                   [-1.83251000    1.55844200   -0.65127900];
                   [-2.51263500   -0.53063900   -0.50722600];
                   [-3.28847400   -0.19858900   -1.08239500];
                   [ 0.45545168    2.23095146   -0.76407436];
                   [-0.46596257   -1.76861142    0.49251865];
                   [ 2.88380973    1.74814805   -0.56768603];
                   [ 1.97436454   -2.26233655    0.71081156];
                   [ 3.65174100   -0.50071826    0.17866774]};
  
% 14N quadrupolar coupling
inter.coupling.matrix=cell(20);  
inter.coupling.matrix{4,4}=1e6*[-2.108e+000  1.299e+000 -1.434e+000
                                 1.299e+000  1.326e+000  2.141e+000
                                -1.434e+000  2.141e+000  7.820e-001];
inter.coupling.matrix{14,14}=1e6*[-2.108e+000  1.299e+000 -1.434e+000
                                   1.299e+000  1.326e+000  2.141e+000
                                  -1.434e+000  2.141e+000  7.820e-001]; 

% Scalar couplings
inter.coupling.scalar=cell(20);  
inter.coupling.scalar{2,1}= -1.122e+001; 
inter.coupling.scalar{3,1}=  5.824e+000; 
inter.coupling.scalar{3,2}=  1.867e+000; 
inter.coupling.scalar{4,1}=  4.404e+000; 
inter.coupling.scalar{4,2}= -5.709e-001; 
inter.coupling.scalar{4,3}=  5.049e+000; 
inter.coupling.scalar{5,1}=  7.729e+000; 
inter.coupling.scalar{5,2}=  1.308e+001; 
inter.coupling.scalar{5,3}=  8.069e+000; 
inter.coupling.scalar{6,3}= -1.481e+000; 
inter.coupling.scalar{7,3}= -1.548e+000; 
inter.coupling.scalar{10,3}=-1.753e+000; 
inter.coupling.scalar{5,4}=  4.456e+001; 
inter.coupling.scalar{8,6}=  1.020e+001; 
inter.coupling.scalar{8,7}=  1.228e+000; 
inter.coupling.scalar{9,6}=  1.249e+000; 
inter.coupling.scalar{9,7}=  1.011e+001; 
inter.coupling.scalar{10,8}= 9.900e+000; 
inter.coupling.scalar{10,9}= 9.964e+000;
inter.coupling.scalar{12,11}= -1.122e+001; 
inter.coupling.scalar{13,11}=  5.824e+000; 
inter.coupling.scalar{13,12}=  1.867e+000; 
inter.coupling.scalar{14,11}=  4.404e+000; 
inter.coupling.scalar{14,12}= -5.709e-001; 
inter.coupling.scalar{14,13}=  5.049e+000; 
inter.coupling.scalar{15,11}=  7.729e+000; 
inter.coupling.scalar{15,12}=  1.308e+001; 
inter.coupling.scalar{15,13}=  8.069e+000; 
inter.coupling.scalar{16,13}= -1.481e+000; 
inter.coupling.scalar{17,13}= -1.548e+000; 
inter.coupling.scalar{20,13}= -1.753e+000; 
inter.coupling.scalar{15,14}=  4.456e+001; 
inter.coupling.scalar{18,16}=  1.020e+001; 
inter.coupling.scalar{18,17}=  1.228e+000; 
inter.coupling.scalar{19,16}=  1.249e+000; 
inter.coupling.scalar{19,17}=  1.011e+001; 
inter.coupling.scalar{20,18}=  9.900e+000; 
inter.coupling.scalar{20,19}=  9.964e+000; 
 
% Chemical shift tensors
inter.zeeman.matrix=cell(1,20);  
inter.zeeman.matrix{1}=  -[3.486e+001  4.615e+000 4.067e+000;  4.615e+000 2.935e+001  1.088e+000; 4.067e+000  1.088e+000 2.604e+001]; 
inter.zeeman.matrix{2}=  -[2.949e+001 -1.201e+000 2.228e+000; -1.201e+000 3.119e+001 -4.875e+000; 2.228e+000 -4.875e+000 2.956e+001];
inter.zeeman.matrix{3}=  -[2.793e+001 -1.604e+000 5.846e+000;  0.124e+000 2.394e+001 -0.814e+000; 2.463e+000 -1.065e+000 3.525e+001];
inter.zeeman.matrix{4}=  -[2.703e+002  2.576e+001 2.657e+001;  2.576e+001 1.598e+002 -9.804e+000; 2.657e+001 -9.804e+000 2.627e+002]; 
inter.zeeman.matrix{5}=  -[3.127e+001  3.101e+000 4.521e+000;  3.101e+000 2.988e+001  5.836e+000; 4.521e+000  5.836e+000 3.643e+001]; 
inter.zeeman.matrix{6}=  -[2.763e+001  3.900e+000 2.213e-001;  3.900e+000 2.405e+001  4.991e-001; 2.213e-001  4.991e-001 2.065e+001]; 
inter.zeeman.matrix{7}=  -[2.772e+001 -3.706e+000 2.839e-001; -3.706e+000 2.450e+001 -5.738e-001; 2.839e-001 -5.738e-001 2.076e+001]; 
inter.zeeman.matrix{8}=  -[2.703e+001 -1.538e+000 9.105e-001; -1.538e+000 2.444e+001 -1.979e-001; 9.105e-001 -1.979e-001 2.145e+001]; 
inter.zeeman.matrix{9}=  -[2.719e+001  1.483e+000 8.995e-001;  1.483e+000 2.445e+001  1.705e-001; 8.995e-001  1.705e-001 2.145e+001]; 
inter.zeeman.matrix{10}= -[2.391e+001  9.055e-002 4.221e-001;  9.055e-002 2.764e+001 -5.050e-003; 4.221e-001 -5.050e-003 2.159e+001];
inter.zeeman.matrix{11}= -[3.486e+001  4.615e+000 4.067e+000;  4.615e+000 2.935e+001  1.088e+000; 4.067e+000  1.088e+000 2.604e+001]; 
inter.zeeman.matrix{12}= -[2.949e+001 -1.201e+000 2.228e+000; -1.201e+000 3.119e+001 -4.875e+000; 2.228e+000 -4.875e+000 2.956e+001];
inter.zeeman.matrix{13}= -[2.793e+001 -1.604e+000 5.846e+000;  0.124e+000 2.394e+001 -0.814e+000; 2.463e+000 -1.065e+000 3.525e+001];
inter.zeeman.matrix{14}= -[2.703e+002  2.576e+001 2.657e+001;  2.576e+001 1.598e+002 -9.804e+000; 2.657e+001 -9.804e+000 2.627e+002]; 
inter.zeeman.matrix{15}= -[3.127e+001  3.101e+000 4.521e+000;  3.101e+000 2.988e+001  5.836e+000; 4.521e+000  5.836e+000 3.643e+001]; 
inter.zeeman.matrix{16}= -[2.763e+001  3.900e+000 2.213e-001;  3.900e+000 2.405e+001  4.991e-001; 2.213e-001  4.991e-001 2.065e+001]; 
inter.zeeman.matrix{17}= -[2.772e+001 -3.706e+000 2.839e-001; -3.706e+000 2.450e+001 -5.738e-001; 2.839e-001 -5.738e-001 2.076e+001]; 
inter.zeeman.matrix{18}= -[2.703e+001 -1.538e+000 9.105e-001; -1.538e+000 2.444e+001 -1.979e-001; 9.105e-001 -1.979e-001 2.145e+001]; 
inter.zeeman.matrix{19}= -[2.719e+001  1.483e+000 8.995e-001;  1.483e+000 2.445e+001  1.705e-001; 8.995e-001  1.705e-001 2.145e+001]; 
inter.zeeman.matrix{20}= -[2.391e+001  9.055e-002 4.221e-001;  9.055e-002 2.764e+001 -5.050e-003; 4.221e-001 -5.050e-003 2.159e+001]; 

% Experimental chemical shifts
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,1:20,[2.237 1.576 2.970 0.0 1.134 7.250 7.250 7.250 7.320 7.320...
                                                        2.060 2.111 2.843 0.0 0.536 7.250 7.250 7.250 7.320 7.320]);
% Chemical kinetics
kplus=1.2e3; kminus=1.2e3;
inter.chem.parts={[1  2  3  4  5  6  7  8  9  10]...
                  [11 12 13 14 15 16 17 18 19 20]};
inter.chem.rates=[-kplus  kminus
                   kplus -kminus];
inter.chem.concs=[kminus kplus];

% Relaxation theory 
inter.relaxation={'redfield','SRFK','SRSK'};
inter.equilibrium='zero';
inter.srsk_sources=[4 14];
inter.rlx_keep='kite';
inter.tau_c={50e-12 50e-12};
inter.srfk_tau_c={[1.0 1/kplus]};
inter.srfk_mdepth=cell(20);
inter.srfk_mdepth{5,1}=15;
inter.srfk_mdepth{5,2}=15;
inter.srfk_mdepth{5,3}=15;
inter.srfk_mdepth{2,3}=3;
inter.srfk_mdepth{1,3}=2.5;
inter.srfk_mdepth{15,11}=15;
inter.srfk_mdepth{15,12}=15;
inter.srfk_mdepth{15,13}=15;
inter.srfk_mdepth{12,13}=3;
inter.srfk_mdepth{11,13}=2.5;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.connectivity='scalar_couplings';
bas.level=4;
bas.space_level=3;

% Disable Krylov algorithm
sys.disable={'krylov'};
sys.enable={'greedy'}; % 'gpu'

% Proximity cut-off
sys.tols.prox_cutoff=10.0;

% Spinach housekeeping 
spin_system=create(sys,inter); 
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.offset=2000;
parameters.sweep=[4000 4000];
parameters.npoints=[128 128];
parameters.zerofill=[512 512];
parameters.spins={'1H'};
parameters.axis_units='ppm';
parameters.tmix=0.800;

% Concentration-aware initial state
parameters.rho0=state(spin_system,'Lz','1H','chem');

% Simulation
fid=liquid(spin_system,@noesy,parameters,'nmr');

% Apodisation
fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}});
fid.sin=apodisation(spin_system,fid.sin,{{'sqcos'},{'sqcos'}});

% F2 Fourier transform
f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1));
f1_sin=real(fftshift(fft(fid.sin,parameters.zerofill(2),1),1));

% States signal
f1_states=f1_cos-1i*f1_sin;

% F1 Fourier transform
spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2);

% Plotting
kfigure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,-real(spectrum),parameters,...
        20,[0.0005 0.01 0.025 0.25],2,256,6,'both');

end

