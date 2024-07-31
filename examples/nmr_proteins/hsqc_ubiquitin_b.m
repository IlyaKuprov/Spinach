% 1H-15N HSQC of human ubiquitin, without 1H decoupling in F1 and
% 15N decoupling in F2. Nitrogen-proton multiplicity is retained
% in both dimensions.
%
% Calculation time: hours, faster with a Tesla A100 GPU.
%
% Zenawi Welderufael
% Luke Edwards
% Ilya Kuprov

function hsqc_ubiquitin_b()

% Protein data import
options.pdb_mol=1;
options.noshift='delete';
options.select='backbone-hsqc';
[sys,inter]=protein('1D3Z.pdb','1D3Z.bmrb',options);

% Magnet field
sys.magnet=14.1;

% Tolerances
sys.tols.inter_cutoff=5.0;
sys.tols.prox_cutoff=4.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.connectivity='scalar_couplings';
bas.level=4; bas.space_level=1;

% Algorithmic options
sys.enable={'greedy','gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.J=90;
parameters.sweep=[2400 4800];
parameters.offset=[-7000 4500];
parameters.npoints=[256 256];
parameters.zerofill=[1024 1024];
parameters.spins={'15N','1H'};
parameters.decouple_f1={'13C'};
parameters.decouple_f2={'13C'};
parameters.axis_units='ppm';

% Simulation
fid=liquid(spin_system,@hsqc,parameters,'nmr');

% Apodization
fid.pos=apodization(fid.pos,'sqcosbell-2d');
fid.neg=apodization(fid.neg,'sqcosbell-2d');

% F2 Fourier transform
f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1);
f1_neg=fftshift(fft(fid.neg,parameters.zerofill(2),1),1);

% Form States signal
fid=f1_pos+conj(f1_neg);

% F1 Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill(1),2),2);

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,abs(spectrum),parameters,...
        20,[0.1 0.5 0.1 0.5],2,256,6,'positive');

end

