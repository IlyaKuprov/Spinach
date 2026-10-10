% Two-spin asymmetric magnetization flux problem.
%
% Calculation time: seconds.
%
% ilya.kuprov@weizmann.ac.il

function flux_asymmetric()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0,3};
inter.chem.parts={1,2};
inter.chem.concs=[2e3 5e2];
inter.chem.reactions={struct('reactants',[1 2],'products',[1 2],...
                            'matching',[1 2;2 1],'rate',1)};

% Basis specification
bas.formalism='sphten-liouv';
bas.approximation={'none','none'};

% Enable zero track elimination
sys.enable={'zte'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=coil_state(spin_system,'L+','1H','exact');
parameters.decouple={};
parameters.offset=900;
parameters.sweep=5000;
parameters.npoints=512;
parameters.zerofill=1024;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@flux_acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

% Freeze additive replacement at its invariant species concentrations
function fid=flux_acquire(spin_system,parameters,H,R,K)

% No concentration-changing reaction is present in this network
K=K(0,unit_state(spin_system));
fid=acquire(spin_system,parameters,H,R,K);

end

