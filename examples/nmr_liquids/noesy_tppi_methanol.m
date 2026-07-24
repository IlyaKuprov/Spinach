% NOESY spectrum of 13C methanol with time-proportional phase
% incrementation (TPPI) acquisition and processing. J-couplings are
% from Pecul and Helgaker, and CSA tensors are from DFT.
%
% Calculation time: seconds
%
% tim.claridge@chem.ox.ac.uk
% ilya.kuprov@weizmann.ac.il

function noesy_tppi_methanol()

% Spin system properties (vacuum DFT calculation)
[sys,inter]=g2spinach(gparse('../standard_systems/methanol.log'),...
                     {{'H','1H'},{'C','13C'}},[31.8 182.4],[]);

% Remove the OH proton
sys.isotopes(end)=[];
inter.coordinates(end)=[];
inter.zeeman.matrix(end)=[];

% Put all chemical shifts on resonance
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,[1 2 3 4],[0 0 0 0]);

% Magnet field
sys.magnet=14.1;

% Assign J-couplings
inter.coupling.scalar=cell(4,4);
inter.coupling.scalar{1,2}=141;
inter.coupling.scalar{1,3}=141;
inter.coupling.scalar{1,4}=141;
inter.coupling.scalar{2,3}=-11;
inter.coupling.scalar{3,4}=-11;
inter.coupling.scalar{2,4}=-11;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='IME';
inter.temperature=298;
inter.rlx_keep='kite';
inter.tau_c={50e-12};

% Algorithmic options
sys.enable={'greedy'};
sys.tols.prox_cutoff=4.0;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.tmix=0.5;
parameters.offset=0;
parameters.sweep=[600 300];
parameters.npoints=[256 256];
parameters.zerofill=[1024 1024];
parameters.spins={'1H'};
parameters.axis_units='ppm';
parameters.needs={'rho_eq'};

% Simulation
fid=liquid(spin_system,@noesy,parameters,'nmr');

% Apply the phase increment to the indirect quadrature channels
phase_inc=(0:(parameters.npoints(1)-1))*pi/2;
fid_tppi=fid.cos.*cos(phase_inc)-fid.sin.*sin(phase_inc);

% Apodise the phase-incremented signal
fid_tppi=apodisation(spin_system,fid_tppi,{{'sqcos'},{'sqcos'}});

% F2 Fourier transform and absorptive component selection
f2_spec=real(fftshift(fft(fid_tppi,parameters.zerofill(2),1),1));

% Real F1 Fourier transform
tppi_spec=fftshift(fft(f2_spec,parameters.zerofill(1),2),2);

% Retain the non-redundant F1 spectral half
spectrum=tppi_spec(:,(parameters.zerofill(1)/2+1):end);

% Restore the physical F1 sweep width
parameters.sweep=[300 300];
parameters.zerofill(1)=parameters.zerofill(1)/2;

% Plotting
kfigure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,-real(spectrum),parameters,...
        20,[0.05 0.5 0.05 0.5],2,256,6,'both');

end


