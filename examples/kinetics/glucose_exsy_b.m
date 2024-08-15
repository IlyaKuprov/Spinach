% 2D EXSY of transmembrane exchange of 3,3-difluoroglucose. See
% the fitting example set for the script that yielded the para-
% meters used below.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk
% dmitry.shishmarev@sydney.edu.au
% philip.kuchel@sydney.edu.au

function glucose_exsy_b()

% Magnet field
sys.magnet=9.3933;

% Isotopes
sys.isotopes={'19F','19F',...
              '19F','19F',...
              '19F','19F',...
              '19F','19F'};
          
% Chemical shifts
inter.zeeman.scalar={-113.6784 -129.6771 ...    % Alpha, inside
                     -113.8796 -129.8002 ...    % Alpha, outside
                     -116.2744 -134.2341 ...    % Beta,  inside
                     -116.4564 -134.3244};      % Beta,  outside
                 
% J-couplings
inter.coupling.scalar=cell(8,8);
inter.coupling.scalar{1,2}=238.0633;            % Alpha, inside
inter.coupling.scalar{3,4}=238.0633;            % Alpha, outside
inter.coupling.scalar{5,6}=239.2315;            % Beta,  inside
inter.coupling.scalar{7,8}=239.2315;            % Beta,  outside

% Cartesian coordinates
inter.coordinates={[-0.0551   -1.2087   -1.6523];      % Alpha, inside
                   [-0.8604   -2.3200   -0.0624];
                                      
                   [-0.0551   -1.2087   -1.6523];      % Alpha, outside
                   [-0.8604   -2.3200   -0.0624];
                   
                   [ 0.3831   -0.4876    1.8363];      % Beta,  inside
                   [ 1.3928   -1.9042    0.6634];
                   
                   [ 0.3831   -0.4876    1.8363];      % Beta,  outside
                   [ 1.3928   -1.9042    0.6634]};
               
% Chemical subsystems
inter.chem.parts={[1  2],...  % Alpha, inside
                  [3  4],...  % Alpha, outside
                  [5  6],...  % Beta,  inside
                  [7  8]};    % Beta,  outside

% Reaction rate matrix
inter.chem.rates=[-0.3438  0.1550  0       0;          
                   0.3438 -0.1550  0       0;         
                   0       0      -0.7995  0.4350; 
                   0       0       0.7995 -0.4350];

% Equilibrium concentrations with alpha-beta imbalance
inter.chem.concs=equilibrate(inter.chem.rates,[3.8034; 0; 14.2442; 0]);
              
% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory parameters
inter.relaxation={'redfield','t1_t2'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={0.9601e-9 ... % Alpha, inside
             0.5255e-9 ... % Alpha, outside
             0.9601e-9 ... % Beta,  inside
             0.5255e-9};   % Beta,  outside
inter.r1_rates={0 0 0 0 0 0 0 0};
inter.r2_rates=34.0359*{1 1 1 1 1 1 1 1};

% Do not draw colorbars
sys.disable={'colorbar'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.tmix=0.5;
parameters.offset=-46681;
parameters.sweep=[8650 8650];
parameters.npoints=[512 1024];
parameters.zerofill=[1024 1024];
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

% Real part
spectrum=-real(spectrum);

% Process the experimental spectrum
load('glucose_expt_b.mat','spec');
expt_spec=keep_rank(atranspose(spec),25);

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

