% Theoretical HN(CO)CA of human ubiquitin. It is assumed that 
% only the backbone is 13C,15N-labelled.
%
% Calculation time: minutes, faster with a Tesla A100 GPU.
%
% ledwards@cbs.mpg.de
% ilya.kuprov@weizmann.ac.il

function hncoca_ubiquitin()

% Protein data import
options.pdb_mol=1;
options.noshift='delete';
options.select='backbone-minimal';
[sys,inter]=protein('1D3Z.pdb','1D3Z.bmrb',options);

% Magnet field
sys.magnet=14.1;

% Tolerances
sys.tols.inter_cutoff=2.0;
sys.tols.prox_cutoff=4.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.connectivity='scalar_couplings';
bas.inter_level=4; bas.prox_level=1;

% Algorithmic options
sys.enable={'greedy'}; % 'gpu'
sys.disable={'krylov'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.tau=[2.25e-3, 2.75e-3, 8.00e-3, 7.00e-3];
parameters.spins={'15N','13C','1H'};
parameters.offset=[-7100 8450 4850];
parameters.sweep=[2500 4500 3000];
parameters.npoints=[64 64 64];
parameters.zerofill=[256 256 256];
parameters.axis_units='ppm';

% Simulation
fid=liquid(spin_system,@hncoca,parameters,'nmr');

% Apodisation
fid.pos_pos=apodisation(spin_system,fid.pos_pos,{{'sqcos'},{'sqcos'},{'sqcos'}});
fid.pos_neg=apodisation(spin_system,fid.pos_neg,{{'sqcos'},{'sqcos'},{'sqcos'}});
fid.neg_pos=apodisation(spin_system,fid.neg_pos,{{'sqcos'},{'sqcos'},{'sqcos'}});
fid.neg_neg=apodisation(spin_system,fid.neg_neg,{{'sqcos'},{'sqcos'},{'sqcos'}});

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
kfigure(); plot_3d(spin_system,-real(spectrum),parameters,...
                   10,[0.2 0.9 0.2 0.9],2,'positive');

end

