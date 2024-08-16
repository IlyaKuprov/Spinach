% Fitting of 2,2,3,3-tetrafluoroglucose NOESY with respect to the 
% reaction rates in a chemical exchange problem and the rotational
% correlation time within Redfield theory.
%
% Calculation time: hours
%
% i.kuprov@soton.ac.uk
% dmitry.shishmarev@sydney.edu.au
% philip.kuchel@sydney.edu.au

function glucose_exsy_a()

% Get a figure going
figure(); scale_figure([1.75 1.00]);

% Set the initial guess
guess=[0.9448    1.6667   -8.4603   -9.0902    1.0074    0.1052 ...
       0.1032    0.0787    0.1125    0.1076    0.1273    0.0636 ...
       0.1367    0.0952    0.1171    0.0818    0.1226    0.1223 ...
       0.1079    0.1548    0.0500   -5.8355  -17.7155   16.7992 ...
     -51.4794   52.2778   -4.8695  -38.2804   32.1033    2.8853 ...
     -38.9973    0.9992    0.8724    1.3650    0.6937   26.6705];

% Set optimiser options
options=optimset('Display','iter','MaxIter',100000,...
                 'MaxFunEvals',200000,'UseParallel',true);

% Run the optimisation
answer=fminsearch(@exsy_err,guess,options);

% Display the result
disp(answer);

% Save figure
try_to_savefig(gcf,'glucose_exsy_a.fig');

end

function err=exsy_err(params)

% Hush up Spinach
sys.output='hush';
sys.disable={'hygiene'};
sys.enable={'greedy'};

% Magnet field
sys.magnet=9.4;

% Isotopes
sys.isotopes={'19F','19F','19F','19F',...
              '19F','19F','19F','19F',...
              '19F','19F','19F','19F',...
              '19F','19F','19F','19F'};
          
% Chemical shifts
inter.zeeman.scalar={-120.5429 -133.9455 -129.3219 -129.5137...    % Alpha, inside
                     -120.8449 -134.3029 -129.6334 -129.7151...    % Alpha, outside
                     -136.4588 -139.5504 -131.8493 -132.1227...    % Beta,  inside
                     -136.8020 -139.7067 -132.2969 -132.1936}...   % Beta,  outside
                     +num2cell(params(6:21)-0.1*ones(1,16));
                 
% J-couplings
inter.coupling.scalar=cell(16,16);

inter.coupling.scalar{1,2}=-265+params(22);
inter.coupling.scalar{3,4}=-265+params(22);
inter.coupling.scalar{1,3}=1.3+params(23);
inter.coupling.scalar{1,4}=4.7+params(24);
inter.coupling.scalar{2,3}=13.1+params(25);
inter.coupling.scalar{2,4}=1.4+params(26);              % Alpha, inside

inter.coupling.scalar{5,6}=-265+params(22);
inter.coupling.scalar{7,8}=-265+params(22);
inter.coupling.scalar{5,7}=1.3+params(23);
inter.coupling.scalar{5,8}=4.7+params(24);
inter.coupling.scalar{6,7}=13.1+params(25);
inter.coupling.scalar{6,8}=1.4+params(26);              % Alpha, outside

inter.coupling.scalar{9,10}=-258+params(27);
inter.coupling.scalar{11,12}=-258+params(27);
inter.coupling.scalar{9,11}=1.0+params(28);
inter.coupling.scalar{9,12}=12.1+params(29);
inter.coupling.scalar{10,11}=3.7+params(30);
inter.coupling.scalar{10,12}=25+params(31);             % Beta,  inside

inter.coupling.scalar{13,14}=-258+params(27);
inter.coupling.scalar{15,16}=-258+params(27);
inter.coupling.scalar{13,15}=1.0+params(28);
inter.coupling.scalar{13,16}=12.1+params(29);
inter.coupling.scalar{14,15}=3.7+params(30);
inter.coupling.scalar{14,16}=25+params(31);             % Beta,  outside

inter.coordinates={[-0.0551   -1.2087   -1.6523];
                   [-0.8604   -2.3200   -0.0624];
                   [-2.4464   -0.1125   -0.9776];
                   [-1.9914   -0.0836    1.0743];
                   
                   [-0.0551   -1.2087   -1.6523];
                   [-0.8604   -2.3200   -0.0624];
                   [-2.4464   -0.1125   -0.9776];
                   [-1.9914   -0.0836    1.0743];
                   
                   [ 0.3831   -0.4876    1.8363];
                   [ 1.3928   -1.9042    0.6634];
                   [ 2.3946    0.8398    0.6811];
                   [ 2.1450   -0.1131   -1.1770];
                   
                   [ 0.3831   -0.4876    1.8363];
                   [ 1.3928   -1.9042    0.6634];
                   [ 2.3946    0.8398    0.6811];
                   [ 2.1450   -0.1131   -1.1770]};
               
% Chemical subsystems
inter.chem.parts={[1  2  3  4],...  % Alpha, inside
                  [5  6  7  8],...  % Alpha, outside
                  [9  10 11 12],... % Beta,  inside
                  [13 14 15 16]};   % Beta,  outside

% Reaction rate matrix
inter.chem.rates=[-params(1)  params(2)  0           0;          
                   params(1) -params(2)  0           0;         
                   0          0         -params(33)  params(34); 
                   0          0          params(33) -params(34)];

% Equilibrium concentrations with alpha-beta imbalance
inter.chem.concs=equilibrate(inter.chem.rates,[6-params(32);
                                               0; 
                                               4+params(32); 
                                               0]);
              
% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory parameters
inter.relaxation={'redfield','t1_t2'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={10^params(3) 10^params(4) 10^params(3) 10^params(4)};
inter.r1_rates=params(35)*num2cell(ones(1,16));  % Other mechanisms
inter.r2_rates=params(36)*num2cell(ones(1,16));  % covered here

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

% Real part
spectrum=-params(5)*real(spectrum);

% Process the experimental spectrum
load('glucose_expt_a.mat','Expression1')
expt_spec=rot90(Expression1,2)/5;

% Compute the least squares error and display parameters
err=1e-6*sum(sum((spectrum-expt_spec).^2)); disp(params);

% Plotting
if ~isworkernode

    % Cosmetic SVD denoising
    expt_spec=keep_rank(expt_spec,50);

    % Theory vs exepriment
    subplot(1,2,1); cla reset; 
    plot_2d(spin_system,spectrum,parameters,...
            20,[0.01 0.2 0.5 1.0],2,256,6,'both');
    ktitle('simulated spectrum');
    subplot(1,2,2); cla reset; 
    plot_2d(spin_system,expt_spec,parameters,...
            20,[0.01 0.2 0.5 1.0],2,256,6,'positive');
    ktitle('experimental spectrum'); drawnow();

end

end

