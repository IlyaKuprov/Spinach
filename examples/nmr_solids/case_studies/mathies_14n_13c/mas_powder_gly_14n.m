% 14N MAS spectrum of glycine powder (assuming decoupling of 1H
% and 13C), computed using the Fokker-Planck MAS formalism and
% a spherical grid. Numerical rotating frame transformation is
% used because 14N quadrupolar interaction is large.
%
% Calculation time: seconds.
%
% sanjay.vinod-kumar@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il

function mas_powder_gly_14n()

% Read CASTEP file
props=c2spinach('glycine.magres');

% Drop H, O, and C atoms
drop_mask=ismember(props.symbols,{'H','O','C'});
props.symbols(drop_mask)=[];
props.std_geom(drop_mask,:)=[];
props.cst(drop_mask)=[];
props.efg(drop_mask)=[];

% keep only 1 14N
sys.isotopes{1}='14N';

% Convert shielding tensor into shift 
inter.zeeman.matrix{1}=-props.cst{1};

% Set isotropic chemical shift to experimental value
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,1,110.0);

% Quadrupolar interaction from CASTEP
nqi=castep2nqi(props.efg{1},20.44e-3,1);
inter.coupling.matrix{1,1}=remtrace(nqi);

% Magnet field
sys.magnet=9.4;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment setup
parameters.sweep=3e6;
parameters.npoints=256;
parameters.zerofill=1024;
parameters.offset=0;
parameters.spins={'14N'};
parameters.grid='rep_2ang_12800pts_sph';
parameters.rho0=state(spin_system,'L+','14N');
parameters.coil=state(spin_system,'L+','14N');
parameters.axis_units='MHz';
parameters.verbose=0;

% Simulation
fid=powder(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); hold on;
plot_1d(spin_system,real(spectrum),parameters);

% Quadrupolar interaction measured by O'Dell PCCP 2009
dell_nqi=2*pi*eeqq2nqi(1.18e6,0.53,1,[0 0 0]);
spin_system.inter.coupling.matrix{1,1}=dell_nqi;

% Simulation
fid=powder(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting and legend
plot_1d(spin_system,real(spectrum),parameters);
klegend({'CASTEP','O''Dell PCCP 2009'});

end

