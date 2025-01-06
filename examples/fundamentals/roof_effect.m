% Roof effect in a strongly J-coupled two-spin system.
%
% ilya.kuprov@weizmann.ac.il

function roof_effect()

% Isotopes
sys.isotopes={'1H','1H'};

% Magnetic induction
sys.magnet=5.9;

% Chemical shifts
inter.zeeman.scalar={0.95 1.45};

% Scalar couplings
inter.coupling.scalar=cell(2,2);
inter.coupling.scalar{1,2}=7.0; 

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=300;
parameters.sweep=300;
parameters.npoints=1024;
parameters.zerofill=4096;
parameters.axis_units='Hz';
parameters.invert_axis=1;

% Get the figure going
figure(); subp=0; scale_figure([1.5 1.0]);

% Loop over line positions
for ppm=[0.2 0.05 0.0125 0.00625]

    % Update the Zeeman frequencies
    spin_system.inter.zeeman.matrix{1}=-eye(3)*spin('1H')*(1+1e-6*(1.2+ppm))*sys.magnet;
    spin_system.inter.zeeman.matrix{2}=-eye(3)*spin('1H')*(1+1e-6*(1.2-ppm))*sys.magnet;

    % Run the simulation
    fid=liquid(spin_system,@acquire,parameters,'nmr');

    % Apodisation
    fid=apodisation(spin_system,fid,{{'exp',10}});

    % Fourier transform
    spectrum=fftshift(fft(fid,parameters.zerofill));

    % Plotting
    subp=subp+1; subplot(2,2,subp);
    plot_1d(spin_system,real(spectrum),parameters);
    kylabel('intensity, a.u.');

end

end

