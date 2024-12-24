% Diffusion attenuation during soft pulses in a simplified model sequence
% of the Zangger-Sterk pure shift iDOSY with fitting using a modified ver-
% sion of the Stejskal Tanner equation, as described in:
%
%                https://doi.org/10.1016/j.jmr.2019.02.010
%
% Calculation time: seconds on NVidia Tesla A100, much longer on CPU
%
% m.g.concilio@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
% gareth.morris@manchester.ac.uk

function idosyzs_test_1()

% Magnetic field
sys.magnet=11.7426;

% Isotopes
sys.isotopes={'1H'};

% Chemical shift
inter.zeeman.scalar={4.6};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.disable={'pt'};
sys.enable={'greedy','gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sample geometry
parameters.dims=0.015;        % m
parameters.npts=4000;
parameters.deriv={'period',7};

% Diffusion coefficient
parameters.diff=18e-10;       % m^2/s

% Relaxation phantom
parameters.rlx_ph={};
parameters.rlx_op={};

% White margins on the initial condition
parameters.rho0_ph={[zeros(1000,1); ones(2000,1); zeros(1000,1)]};                                                                                            
parameters.rho0_st={state(spin_system,'Lz','1H','cheap')};
parameters.coil_ph={ones(parameters.npts,1)};
parameters.coil_st={state(spin_system,'L+','1H','cheap')};

% Sequence parameters
parameters.spins={'1H'};                % Working spins
parameters.offset=2500;                 % Transmitter offset
parameters.rf_phi=pi;                   % Phase of the inversion pulse
parameters.filename='gaussian_1000.pk'; % Soft pulse shape
parameters.pulse_npoints=100;           % Number of points in the soft pulse
parameters.rf_dur=0.045;                % Soft pulse duration
parameters.delta_sml=0.002;             % Gradient duration
parameters.delta_big=0.1;               % Diffusion delay duration   
parameters.sel_g_amp=0.0053;            % Amplitude of the ZS gradient
   
% Gradient amplitude list 
grad_amps=linspace(0.01,0.40,20); % T/m   
         
% Compute a series of spectra with different gradient amplitudes
parfor n=1:numel(grad_amps)
    
    % Set gradient amplitude   
    localpar=parameters; localpar.g_amp=grad_amps(n);

    % Run the simulation
    inten(n)=imaging(spin_system,@idosyzs,localpar);
              
end

% Normalise and gather intensities
inten=gather(inten/inten(1));

% Stejskal-Tanner fitting function
expfactor=(spin(parameters.spins{1})*parameters.delta_sml).^2*...
          (parameters.delta_big-(parameters.delta_sml)/3)*1e-10;
F=@(x,grad_amps)x(1)*exp(-x(2)*expfactor*(grad_amps-x(3)).^2);

% Run Stejskal-Tanner fitting
opts=optimoptions(@lsqcurvefit,'Display','iter');   
x=lsqcurvefit(F,[1 18 0.01],grad_amps,inten,[],[],opts);  

% Plot simulation data
figure(); plot(grad_amps,inten,'ro'); hold on;

% Plot fitting data
g_axis=linspace(min(grad_amps),max(grad_amps),100);
plot(g_axis,F(x,g_axis),'b-'); axis tight; kgrid;
xlabel('gradient amplitude, T/m'); 
legend({'simulation','fitting'});

% Report the parameters
disp(['Fitted D:       ' num2str(x(2)*1e-10)]);
disp(['Gradient shift: ' num2str(x(3))]);

end
