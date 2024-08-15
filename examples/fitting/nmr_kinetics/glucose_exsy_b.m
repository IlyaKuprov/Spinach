% Fitting of 3,3-difluoroglucose NOESY with respect to the reaction
% rates in a chemical exchange and the rotational correlation times
% within Redfield theory.
%
% Calculation time: hours
%
% i.kuprov@soton.ac.uk
% dmitry.shishmarev@sydney.edu.au
% philip.kuchel@sydney.edu.au

function glucose_exsy_b()

% Get a figure going
figure(); scale_figure([1.75 1.00]);

% Set the initial guess
guess=[ 0.2331    0.1114    6.8446    0.1210    0.1423    0.1205  ...
        0.1697    0.0755    0.1158    0.0936    0.1256   26.9142  ...
       18.7751   25.5476    0.7948    0.4331   -9.0277   -9.2828  ...
        1.1924   32.9179];

% Set optimiser options
options=optimset('Display','iter','MaxIter',100000,...
                 'MaxFunEvals',200000,'UseParallel',true);

% Run the optimisation
answer=fminsearch(@exsy_err,guess,options);

% Display the result
disp(answer);

% Save figure
savefig(gcf,'glucose_exsy_b.fig');

end

function err=exsy_err(params)

% Hush up Spinach
sys.output='hush';
sys.disable={'hygiene'};
sys.enable={'greedy'};

% Magnet field
sys.magnet=9.3933;

% Isotopes
sys.isotopes={'19F','19F',...
              '19F','19F',...
              '19F','19F',...
              '19F','19F'};
          
% Chemical shifts
inter.zeeman.scalar={-113.70 -129.72 ...    % Alpha, inside
                     -113.90 -129.87 ...    % Alpha, outside
                     -116.25 -134.25 ...    % Beta,  inside
                     -116.45 -134.35}...    % Beta,  outside
                     +num2cell(params(4:11)-0.1*ones(1,8));
                 
% J-couplings
inter.coupling.scalar=cell(8,8);

inter.coupling.scalar{1,2}=-265+params(12);            % Alpha, inside
inter.coupling.scalar{3,4}=-265+params(12);            % Alpha, outside
inter.coupling.scalar{5,6}=-258+params(13);            % Beta,  inside
inter.coupling.scalar{7,8}=-258+params(13);            % Beta,  outside

inter.coordinates={[-0.0551   -1.2087   -1.6523];
                   [-0.8604   -2.3200   -0.0624];
                                      
                   [-0.0551   -1.2087   -1.6523];
                   [-0.8604   -2.3200   -0.0624];
                   
                   [ 0.3831   -0.4876    1.8363];
                   [ 1.3928   -1.9042    0.6634];
                   
                   [ 0.3831   -0.4876    1.8363];
                   [ 1.3928   -1.9042    0.6634]};
               
% Chemical subsystems
inter.chem.parts={[1  2],...  % Alpha, inside
                  [3  4],...  % Alpha, outside
                  [5  6],...  % Beta,  inside
                  [7  8]};    % Beta,  outside

% Reaction rate matrix
inter.chem.rates=[-params(1)  params(2)  0           0;          
                   params(1) -params(2)  0           0;         
                   0          0         -params(15)  params(16); 
                   0          0          params(15) -params(16)];

% Equilibrium concentrations with alpha-beta imbalance
inter.chem.concs=equilibrate(inter.chem.rates,[params(3);
                                               0; 
                                               params(14); 
                                               0]);
              
% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory parameters
inter.relaxation={'redfield','t1_t2'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={10^params(17) 10^params(18) ...
             10^params(17) 10^params(18)};
inter.r1_rates=params(19)*num2cell(ones(1,8));  % Other mechanisms
inter.r2_rates=params(20)*num2cell(ones(1,8));  % covered here

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
expt_spec=atranspose(spec);

% Compute the least squares error and display parameters
err=1e-6*sum(sum((spectrum-expt_spec).^2)); disp(params);

% Plotting
if ~isworkernode

    % Cosmetic SVD denoising
    expt_spec=keep_rank(expt_spec,25);

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

