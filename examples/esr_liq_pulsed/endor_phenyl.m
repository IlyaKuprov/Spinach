% Mims ENDOR on a phenyl radical in liquid state. The g-factor and the
% isotropic proton hyperfine couplings are specified explicitly from the
% experimental data of Kasai, Hedaya, and Whipple (J. Am. Chem. Soc.
% 1969, 91, 4364): a(ortho)=17.4 G, a(meta)=5.9 G, a(para)=1.9 G, and an
% isotropic g-factor of 2.0024. The two ortho and the two meta protons
% are magnetically equivalent, so the full symmetry treatment uses an
% S2 x S2 group direct product.
%
% Calculation time: seconds
%
% ledwards@cbs.mpg.de
% ilya.kuprov@weizmann.ac.il

function endor_phenyl()

% Magnet field
sys.magnet=0.33;

% Electron and the five ring protons (two ortho, two meta, one para)
sys.isotopes={'E','1H','1H','1H','1H','1H'};

% Phenyl radical isotropic g-factor
g_phenyl=2.0024;
inter.zeeman.scalar={g_phenyl 0 0 0 0 0};

% Isotropic proton hyperfine couplings, converted from milliTesla
% using the phenyl radical g-factor
inter.coupling.scalar=cell(6,6);
inter.coupling.scalar{2,1}=mt2hz(1.74,g_phenyl);
inter.coupling.scalar{3,1}=mt2hz(1.74,g_phenyl);
inter.coupling.scalar{4,1}=mt2hz(0.59,g_phenyl);
inter.coupling.scalar{5,1}=mt2hz(0.59,g_phenyl);
inter.coupling.scalar{6,1}=mt2hz(0.19,g_phenyl);

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Symmetry: equivalent ortho pair and equivalent meta pair
bas.sym_group={'S2','S2'};
bas.sym_spins={[2 3],[4 5]};

% Sequence parameters
parameters.offset=0;
parameters.npoints=512;
parameters.sweep=3e8;
parameters.tau=100e-9;
parameters.zerofill=4096;
parameters.spins={'E'};
parameters.axis_units='MHz';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@endor_mims,parameters,'esr');

% Crude apodisation
fid=apodisation(spin_system,fid-mean(fid),{{'kaiser',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,abs(spectrum),parameters);

end

