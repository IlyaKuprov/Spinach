% Powder NQR nutation curve for a system with a 
% single 127I nucleus.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function nutation_nqr_iodine()

% System specification
sys.magnet=0;
sys.isotopes={'127I'};
inter.coupling.matrix{1,1}=eeqq2nqi(560e6,0.01,5/2,[0 0 0]);

% Formalism and basis
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory
inter.relaxation={'damp'};
inter.damp_rate=1e5;
inter.rlx_keep='labframe';
inter.equilibrium='zero';
inter.temperature=298;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.spins={'127I'};
parameters.needs={'aniso_eq'};
parameters.sweep=[83.5e6 84.5e6];
parameters.npoints=512;
parameters.grid='rep_2ang_200pts_sph';
parameters.coil=state(spin_system,'L+','127I');
parameters.Lx=(operator(spin_system,'L+','127I')+...
               operator(spin_system,'L-','127I'))/2;
parameters.Ly=(operator(spin_system,'L+','127I')-...
               operator(spin_system,'L-','127I'))/2i;
parameters.rf_frq=84.0e6;
parameters.rf_pwr=2*pi*1e5;

% Get a figure started
figure(); scale_figure([2 1]); drawnow;

% Loop over the pulse durations
for n=1:10
    
    % Set pulse duration
    params=parameters; params.rf_dur=5e-7*n;

    % Run the simulation
    spectrum=powder(spin_system,@nqr_pa,params,'labframe');
    
    % Demodulate the transmitter offset
    spectrum=exp(-1i*(2*pi*parameters.rf_frq)*params.rf_dur)*spectrum;

    % Plotting
    frq_axis=linspace(parameters.sweep(1),...
                      parameters.sweep(2),...
                      parameters.npoints);
    subplot(1,10,n); plot(1e-6*frq_axis',imag(spectrum));
    ktitle([num2str(n/2) ' $\mu$s']); kxlabel('MHz'); 
    xlim tight; ylim([-3e-11 3e-11]); 
    set(gca,'YTick',[]); drawnow();
 
end

end

