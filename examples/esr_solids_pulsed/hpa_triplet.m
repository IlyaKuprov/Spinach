% Hypothetical powder averaged X-band pulse-acquire ESR 
% spectrum of photogenerated pentacene triplet state.
%
% Calculation time: seconds.
%
% i.kuprov@soton.ac.uk
% guinevere.mathies@uni-konstanz.de

function hpa_triplet()

% Magnet field
sys.magnet=0.33;

% Triplet electron
sys.isotopes={'E3'};

% Zeeman tensor, assumed isotropic
inter.zeeman.matrix={diag([2.0 2.0 2.0])};

% ZFS, photo-excited pentacene triplet
D=1360.1*1e6; E=-47.2*1e6;    % [Hz]
inter.coupling.matrix={zfs2mat(D,E,0,0,0)};

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Disable trajectory-level SSR algorithms
sys.disable={'trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'E3'};
parameters.coil=state(spin_system,'L+','E3');
parameters.pulse_op=operator(spin_system,'Ly','E3');
parameters.pulse_angle=pi/4;
parameters.offset=0;
parameters.sweep=4e9;
parameters.npoints=128;
parameters.zerofill=512;
parameters.axis_units='GHz-labframe';
parameters.grid='rep_2ang_6400pts_sph';
parameters.derivative=0;
parameters.invert_axis=1;

% Zeeman tensor into Hz/Tesla
Z=-spin('E')*inter.zeeman.matrix{1,1}/(2*pi*2.00231930436256);

% Orientation-dependent intial condition
parameters.rho0=@(alp,bet,gam)zftrip(spin_system,euler2dcm(alp,bet,gam)*...
                                                 inter.coupling.matrix{1,1}*... 
                                                 euler2dcm(alp,bet,gam)',...
                                                 [0.56 0.31 0.13],...
                                                 euler2dcm(alp,bet,gam)*Z*...
                                                 euler2dcm(alp,bet,gam)',sys.magnet,1);

% Simulation of a pulse-acquire experiment
fid=powder(spin_system,@hp_acquire,parameters,'esr');

% Apodization
fid=apodization(fid,'crisp-1d');

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

