% 3D HCCH COSY experiment on GB1 protein.
%
% Calculation time: hours.
%
% m.walker@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function hcch_cosy_gb1()

% Protein data import
options.pdb_mol=1;
options.noshift='delete';
options.select='all';
[sys,inter]=protein('2N9K.pdb','2N9K.bmrb',options);

% Magnet field
sys.magnet=14.1;

% Tolerances
sys.tols.inter_cutoff=20.0;
sys.tols.prox_cutoff=4.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.connectivity='scalar_couplings';
bas.level=4; bas.space_level=1;

% Algorithmic options
sys.enable={'greedy','prop_cache'};

% Sequence parameters
parameters.J_ch=140;
parameters.J_cc=35;
parameters.delta=1.1e-3;
parameters.sweep=[6000 13000 6000];
parameters.spins={'1H','13C','1H'};
parameters.offset=[2500 7000 2500];
parameters.npoints=[128 128 128];
parameters.zerofill=[256 256 256];
parameters.decouple_f3={'13C'};
parameters.axis_units='ppm';

% Create the spin system structure
spin_system=create(sys,inter);

% Kill nitrogens (not relevant)
spin_system=kill_spin(spin_system,strcmp('15N',spin_system.comp.isotopes));

% Build the basis
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@hcch_cosy,parameters,'nmr');

% Apodisation
fid.pos_pos=apodisation(spin_system,fid.pos_pos,{{'cos'},{'cos'},{'cos'}});
fid.pos_neg=apodisation(spin_system,fid.pos_neg,{{'cos'},{'cos'},{'cos'}});
fid.neg_pos=apodisation(spin_system,fid.neg_pos,{{'cos'},{'cos'},{'cos'}});
fid.neg_neg=apodisation(spin_system,fid.neg_neg,{{'cos'},{'cos'},{'cos'}});

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
kfigure(); plot_3d(spin_system,real(spectrum),parameters,...
                  10,[0.05 0.25 0.05 0.25],2,'positive');

end

