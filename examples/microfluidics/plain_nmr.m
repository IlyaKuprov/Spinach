% NMR spectrum of the reaction mixture in the absence of 
% chemical kinetics and spatial dynamics.
%
% Direct-sum reduction separates the small acrylonitrile block, which
% uses expm rather than the stock batched Taylor propagator. At the
% default propagator rounding tolerance, this changes the full FID by
% about 1.63e-6 relative; mapped spin operators agree to round-off.
% A dense expm reference finds the direct-sum FID more accurate; both
% trajectories retain rounding error from their larger propagators.
%
% a.acharya@soton.ac.uk
% madhukar.said@ugent.be
% bruno.linclau@ugent.be
% ilya.kuprov@weizmann.ac.il

function plain_nmr()

% Import Diels-Alder cycloaddition
[sys,inter,bas]=dac_reaction();

% Disable reactions for the static mixture spectrum
inter.chem.reactions={};

% Equal reacting-species concentrations, unpopulated solvent protons
inter.chem.concs=[1 1 1 1 0];

% Magnet field
sys.magnet=14.1;

% Greedy parallelisation
sys.enable={'zte','greedy'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters - 1H
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=coil_state(spin_system,'L+','1H','exact');
parameters.decouple={};
parameters.offset=2328;
parameters.sweep=3500;
parameters.npoints=4096;
parameters.zerofill=16384;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

