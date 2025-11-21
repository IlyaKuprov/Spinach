% 1H-13C HSQC spectrum of the example RNA molecule provided by
% the Wagner group.
%
% Calculation time: minutes
%
% Shunsuke Imai
% Scott Robson
% Gerhard Wagner
% Zenawi Welderufael
% Ilya Kuprov

function rna_hsqc_theo()

% Import RNA data
options.noshift='delete'; options.deut_list={};
[sys,inter]=nuclacid('example.pdb','example.txt',options);

% Magnet field
sys.magnet=17.62;

% Tolerances
sys.tols.inter_cutoff=5.0;
sys.tols.prox_cutoff=5.0;
sys.disable={'krylov','colorbar'};
sys.enable={'greedy'};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.connectivity='scalar_couplings';
bas.level=4; bas.space_level=1;

% Relaxation theory
inter.relaxation={'damp'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.damp_rate=5.0;

% Sequence parameters
parameters.J=90;
parameters.sweep=[2500 2000];
parameters.offset=[22000 6000];
parameters.npoints=[128 256];
parameters.zerofill=[1024 1024];
parameters.spins={'13C','1H'};
parameters.decouple_f1={'1H'};
parameters.decouple_f2={'13C'};
parameters.axis_units='ppm';

% Create the spin system structure
spin_system=create(sys,inter);

% Build the basis
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@hsqc,parameters,'nmr');

% Apodisation
fid.pos=apodisation(spin_system,fid.pos,{{'cos'},{'cos'}});
fid.neg=apodisation(spin_system,fid.neg,{{'cos'},{'cos'}});

% F2 Fourier transform
f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1);
f1_neg=fftshift(fft(fid.neg,parameters.zerofill(2),1),1);

% Form States signal
fid=f1_pos+conj(f1_neg);

% F1 Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill(1),2),2);

% Plotting
kfigure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,real(spectrum),parameters,...
        20,[0.1 0.5 0.1 0.5],2,256,6,'positive');

end

