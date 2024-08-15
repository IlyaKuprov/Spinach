% 2D EXSY of transmembrane exchange of 2,2,3,3-tetrafluoroglucose. See
% the fitting example set for the script that yielded the parameters
% used below.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk
% dmitry.shishmarev@sydney.edu.au
% philip.kuchel@sydney.edu.au

function glucose_exsy_a()

% Magnet field
sys.magnet=9.4;

% Isotopes
sys.isotopes={'19F','19F','19F','19F',...
              '19F','19F','19F','19F',...
              '19F','19F','19F','19F',...
              '19F','19F','19F','19F'};

% Chemical shifts
inter.zeeman.scalar={-120.5380 -133.9429 -129.3169 -129.5320 ...    % Alpha, inside
                     -120.8361 -134.2763 -129.5988 -129.7528 ...    % Alpha, outside
                     -136.4639 -139.5340 -131.8536 -132.1196 ...    % Beta,  inside
                     -136.7790 -139.6982 -132.2901 -132.1922};      % Beta,  outside

% J-couplings
inter.coupling.scalar=cell(16,16);

inter.coupling.scalar{1,2}=  271.2924;
inter.coupling.scalar{3,4}=  271.2924;
inter.coupling.scalar{1,3}=    0.5401;
inter.coupling.scalar{1,4}=  -25.9884;
inter.coupling.scalar{2,3}=    9.9625;
inter.coupling.scalar{2,4}=  -40.7675;       % Alpha, inside

inter.coupling.scalar{5,6}=  271.2924;
inter.coupling.scalar{7,8}=  271.2924;
inter.coupling.scalar{5,7}=    0.5401;
inter.coupling.scalar{5,8}=  -25.9884;
inter.coupling.scalar{6,7}=    9.9625;
inter.coupling.scalar{6,8}=  -40.7675;       % Alpha, outside

inter.coupling.scalar{9,10}=   263.2659;
inter.coupling.scalar{11,12}=  263.2659;
inter.coupling.scalar{9,11}=     1.5172;
inter.coupling.scalar{9,12}=   -28.4799;
inter.coupling.scalar{10,11}=    3.6493;
inter.coupling.scalar{10,12}=  -22.9903;     % Beta,  inside

inter.coupling.scalar{13,14}=  263.2659;
inter.coupling.scalar{15,16}=  263.2659;
inter.coupling.scalar{13,15}=    1.5172;
inter.coupling.scalar{13,16}=  -28.4799;
inter.coupling.scalar{14,15}=    3.6493;
inter.coupling.scalar{14,16}=  -22.9903;     % Beta,  outside

% Cartesian coordinates
inter.coordinates={[-0.0551   -1.2087   -1.6523];
                   [-0.8604   -2.3200   -0.0624];
                   [-2.4464   -0.1125   -0.9776];
                   [-1.9914   -0.0836    1.0743];     % Alpha, inside
                   
                   [-0.0551   -1.2087   -1.6523];
                   [-0.8604   -2.3200   -0.0624];
                   [-2.4464   -0.1125   -0.9776];
                   [-1.9914   -0.0836    1.0743];     % Alpha, outside
                   
                   [ 0.3831   -0.4876    1.8363];
                   [ 1.3928   -1.9042    0.6634];
                   [ 2.3946    0.8398    0.6811];
                   [ 2.1450   -0.1131   -1.1770];     % Beta,  inside
                   
                   [ 0.3831   -0.4876    1.8363];
                   [ 1.3928   -1.9042    0.6634];
                   [ 2.3946    0.8398    0.6811];
                   [ 2.1450   -0.1131   -1.1770]};    % Beta,  outside
               
% Chemical subsystems
inter.chem.parts={[1  2  3  4],...  % Alpha, inside
                  [5  6  7  8],...  % Alpha, outside
                  [9  10 11 12],... % Beta,  inside
                  [13 14 15 16]};   % Beta,  outside

% Reaction rate matrix
inter.chem.rates=[-1.0045  1.7738  0       0;          
                   1.0045 -1.7738  0       0;         
                   0       0      -0.9304  1.4586; 
                   0       0       0.9304 -1.4586];

% Equilibrate translocation with alpha-beta imbalance as the start
inter.chem.concs=equilibrate(inter.chem.rates,[3.2258; 0; 3.1902; 0]);
              
% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={4.137e-9 ...  % Alpha, inside
             0.951e-9 ...  % Alpha, outside
             4.137e-9 ...  % Beta,  inside
             0.951e-9};    % Beta,  outside

% Do not draw colorbars
sys.disable={'colorbar'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.tmix=0.5;
parameters.offset=-49000;
parameters.sweep=[8000 8000];
parameters.npoints=[256 256];
parameters.zerofill=[1024 512];
parameters.spins={'19F'};
parameters.axis_units='ppm';
parameters.rho0=state(spin_system,'Lz','19F','chem');

% Simulation
fid=liquid(spin_system,@noesy,parameters,'nmr');

% Apodization
fid.cos=apodization(fid.cos,'sqcosbell-2d');
fid.sin=apodization(fid.sin,'sqcosbell-2d'); 

% F2 Fourier transform
f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1));
f1_sin=real(fftshift(fft(fid.sin,parameters.zerofill(2),1),1));

% States signal
f1_states=f1_cos-1i*f1_sin;

% F1 Fourier transform
spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2);

% Real part and sign flip
spectrum=-real(spectrum);

% Process the experimental spectrum
load('glucose_expt_a.mat','Expression1')
expt_spec=rot90(Expression1,2)/5;
expt_spec=keep_rank(expt_spec,50);

% Plotting - theory vs experiment
figure(); scale_figure([2.1 1.2]);
subplot(1,2,1);  
plot_2d(spin_system,spectrum,parameters,...
        20,[0.01 0.2 0.01 0.2],2,256,6,'positive');
ktitle('simulated spectrum');
subplot(1,2,2); 
plot_2d(spin_system,expt_spec,parameters,...
        20,[0.01 0.2 0.01 0.2],2,256,6,'positive');
ktitle('experimental spectrum');

% Plotting - deviation histogram
figure(); histogram(spectrum(:)-expt_spec(:));
kxlabel('deviation'); kylabel('point count');
ktitle('difference histogram');
kgrid; xlim([-300 300]);

end

