% 1H-1H-15N NOESY-HSQC spectrum of 15N-labelled ubiquitin at 900
% MHz with 65 ms mixing time. It is assumed that the protein is 
% not 13C-labelled. 
%
% Calculation time: hours.
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk

function noesyhsqc_ubiquitin_prot()

% Protein data import
options.pdb_mol=1;
options.select='all';
options.noshift='delete';
[sys,inter]=protein('1D3Z.pdb','1D3Z.bmrb',options);

% Magnet field
sys.magnet=21.1356;

% Tolerances
sys.tols.inter_cutoff=2.0;
sys.tols.prox_cutoff=4.0;

% Relaxation theory
inter.relaxation={'redfield'};
inter.rlx_keep='kite';
inter.equilibrium='zero';
inter.tau_c={5e-9};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.connectivity='scalar_couplings';
bas.level=4; bas.space_level=3;

% Algorithmic options
sys.enable={'prop_cache','greedy'};
sys.disable={'asyredf'};

% Create the spin system structure
spin_system=create(sys,inter);

% Kill carbons
spin_system=kill_spin(spin_system,strcmp('13C',spin_system.comp.isotopes));

% Build the basis
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.tmix=0.065;
parameters.J=90.0;
parameters.npoints=[128 64 128];
parameters.zerofill=[512 256 512];
parameters.spins={'1H','15N','1H'};
parameters.offset=[4250 -10600 4250];
parameters.sweep=[10750 3000 10750];
parameters.axis_units='ppm';

% Simulation
fid=liquid(spin_system,@noesyhsqc,parameters,'nmr');

% Apodization
fid.pos_pos=apodization(fid.pos_pos,'sqcosbell-3d');
fid.pos_neg=apodization(fid.pos_neg,'sqcosbell-3d');
fid.neg_pos=apodization(fid.neg_pos,'sqcosbell-3d');
fid.neg_neg=apodization(fid.neg_neg,'sqcosbell-3d');

% F3 Fourier transform
f3_pos_pos=fftshift(fft(fid.pos_pos,parameters.zerofill(3),3),3);
f3_pos_neg=fftshift(fft(fid.pos_neg,parameters.zerofill(3),3),3);
f3_neg_pos=fftshift(fft(fid.neg_pos,parameters.zerofill(3),3),3);
f3_neg_neg=fftshift(fft(fid.neg_neg,parameters.zerofill(3),3),3);

% Absorption part of F3 signal
f3_pos=f3_pos_pos+conj(f3_neg_neg);
f3_neg=f3_neg_pos+conj(f3_pos_neg);

% F2 Fourier transform 
f3f2_pos=fftshift(fft(f3_pos,parameters.zerofill(2),2),2);
f3f2_neg=fftshift(fft(f3_neg,parameters.zerofill(2),2),2);

% Absorption part of F2 signal
f3f2=f3f2_pos+conj(f3f2_neg);

% F1 Fourier transform
spectrum=fftshift(fft(f3f2,parameters.zerofill(1),1),1);

% Plotting
figure(); plot_3d(spin_system,-real(spectrum),parameters,...
                  10,[0.01 0.1 0.01 0.1],2,'positive');

end

