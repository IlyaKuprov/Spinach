% Powder-averaged two-pulse out-of-phase ESEEM on the [P700+,A1-]
% spin-correlated electron pair in Photosystem I. Time-domain si-
% mulation in Liouville space with averaging over a powder grid.
% Set to reproduce Figure 3a in 
%
%    http://dx.doi.org/10.1021/bi048445d
%
% The eigenvalues of P700+ g-tensor (without angles) come from
%
%    http://dx.doi.org/10.1016/S0005-2728(01)00198-0
%
% The eigenvales of A- g-tensor (without angles) come from
%
%    http://dx.doi.org/10.1007/BF00019589
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function oop_eseem_photosystem_1()

% Magnet field
sys.magnet=0.3249;

% System specification
sys.isotopes={'E','E'};
inter.zeeman.eigs{1}=[2.00304 2.00262 2.00232];
inter.zeeman.euler{1}=[0 0 0];
inter.zeeman.eigs{2}=[2.00670 2.00560 2.00240];
inter.zeeman.euler{2}=[0 0 0];
inter.coupling.scalar=cell(2,2);
inter.coupling.scalar{1,2}=mt2hz(0.6e-3);
inter.coordinates={[ 0.00  0.00  0.00];
                   [25.35  0.00  0.00]};

% Relaxation theory
inter.relaxation={'damp'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.damp_rate=1.7e6;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Disable trajectory-level SSR algorithms
sys.disable={'trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set the sequence parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,{'Lz','Lz'},{1,2});
parameters.coil=state(spin_system,'L+','E');
parameters.pulse_op=(operator(spin_system,'L+','E')-...
                     operator(spin_system,'L-','E'))/2i;
parameters.offset=0;
parameters.npoints=200;
parameters.timestep=2e-8;
parameters.grid='rep_2ang_1600pts_sph';

% Simulation
fid=powder(spin_system,@oopeseem,parameters,'esr');

% Build the axis
time_axis=(0:(parameters.npoints-1))*parameters.timestep*1e6/2;

% Plot the results
kfigure(); plot(time_axis,-imag(fid)); 
kxlabel('time, $\mu$s'); xlim tight;
kylabel('echo intensity, a.u.'); kgrid;

end

