% Example of a 2D-INADEQUATE spectrum of a 
% generic three-spin system.
%
% Calculation time: seconds.
%
% Theresa Hune
% Christian Griesinger

function inad_three_spin_2d()

% Magnetic field (700 MHz)
sys.magnet=16.44;

% Generic three-spin system
sys.isotopes={'13C','13C','13C'};
inter.zeeman.scalar={10 30 70};
inter.coupling.scalar{1,2}=20;
inter.coupling.scalar{1,3}=60;
inter.coupling.scalar{2,3}=0;
inter.coupling.scalar{3,3}=0;

% Formalism and basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Generate isotopomers
subsystems=dilute(spin_system,'13C',2);

% Sequence parameters
parameters.offset=17604.78;             % Offset at 100 ppm
parameters.J=50;
parameters.spins={'13C'};
parameters.sweep=[35213.086 34722.223]; % Sweep width of 200 ppm
parameters.npoints=[128 2048];
parameters.zerofill=[512 8192];
parameters.axis_units='ppm';

% Preallocate the answer
spectrum=zeros(parameters.zerofill(2),...
               parameters.zerofill(1));

% Loop over isotopomers
parfor n=1:numel(subsystems)
   
    % Build the basis
    subsystem=basis(subsystems{n},bas);
    
    % Simulation
    fid=liquid(subsystem,@inadequate_2d,parameters,'nmr');
    
    % Apodisation
    fid.cos=apodisation(spin_system,fid.cos,{{'cos'},{'cos'}});
    fid.sin=apodisation(spin_system,fid.sin,{{'cos'},{'cos'}});

    % F2 Fourier transform
    spec_cos=fftshift(fft(fid.cos,parameters.zerofill(2),1),1);
    spec_sin=fftshift(fft(fid.sin,parameters.zerofill(2),1),1);

    % Form States signal
    spec_states=real(spec_cos)+1i*real(spec_sin);

    % F1 Fourier transform
    spectrum=spectrum+fftshift(fft(spec_states,parameters.zerofill(1),2),2);

end

% Do the plotting
figure(); scale_figure([1.5 2.0]); 
plot_2d(spin_system,real(spectrum),parameters,...
        20,[0.01 0.5 0.01 0.5],2,256,6,'both');
kylabel('F1: DQ dimension / ppm');

end

