% Magic angle is usually calibrated using KBr powder. When
% the angle is not correctly set, the spinning sideband pat-
% tern is blurred. This simulation demonstrates the effect.
%
% Calculation time: seconds
%
% guinevere.mathies@uni-konstanz.de

function magic_angle_calibration()

% Magnet field
sys.magnet=9.4;

% Isotopes
sys.isotopes={'79Br'};

% Quardupolar coupling tensor
inter.coupling.matrix{1,1}=eeqq2nqi(-92.4e3,-0.79,3/2,[0 0 0]);

% Chemical shift
inter.zeeman.scalar={60.0933};

% Formalism and basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.projections=+1;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment setup
parameters.spins={'79Br'};
parameters.rate=4000;
parameters.max_rank=25;
parameters.grid='rep_2ang_1600pts_sph';
parameters.sweep=1e5;
parameters.npoints=2048;
parameters.zerofill=4096;
parameters.offset=6000;
parameters.axis_units='Hz';
parameters.invert_axis=1;
parameters.rho0=state(spin_system,'L+','79Br');
parameters.coil=state(spin_system,'L+','79Br');

% Convert magic angle errors to radians
ma_errors=deg2rad([-1.00 -0.25  0.00  0.25  1.00]);

% Get a figure going and compute the time axis
figure(); scale_figure([2.5 2.5]);
time_axis=(0:1:parameters.npoints-1)/parameters.sweep;

% Loop over angle errors
for n=1:numel(ma_errors)
    
    % Tilt the spinnig axis away from the magic angle
    parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]*euler2dcm(0,ma_errors(n),0);

    % Run the MAS simulation
    fid=singlerot(spin_system,@acquire,parameters,'nmr');

    % Apodization
    fid=apodization(fid,'exp-1d',6);
    
    % Fourier transform
    spectrum=fftshift(fft(fid,parameters.zerofill));
    
    % Plot the time domain signal
    subplot(numel(ma_errors),2,2*n-1);
    plot(time_axis,real(fid)); kgrid; xlim tight;
    klegend(['error: ' num2str(rad2deg(ma_errors(n))),' deg']);
    kxlabel('time, seconds'); drawnow;
    
    % Plot the spectrum
    subplot(numel(ma_errors),2,2*n);
    plot_1d(spin_system,real(spectrum),parameters);
    klegend(['error: ' num2str(rad2deg(ma_errors(n))),' deg']);
    kgrid; xlim tight; kxlabel('frequency, Hz'); 
    ylim(get(gca,'YLim')/10); drawnow;
   
end

end

