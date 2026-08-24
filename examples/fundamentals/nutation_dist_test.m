% Recovery of an RF field distribution from a nutation curve. A
% proton ensemble is driven on-resonance by a bimodal distribution
% of RF field amplitudes; the two magnetisation components that
% rotate under the RF field are combined into a complex nutation
% curve, noise is added, and nutation_dist.m is called to recover
% the nutation frequency distribution.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function nutation_dist_test()

% Isotopes
sys.isotopes={'1H'};

% Magnet field
sys.magnet=14.1;

% Chemical shift
inter.zeeman.scalar={0.0};

% No console output
sys.output='hush';

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Initial and observable states
rho=state(spin_system,'Lz','1H');
coil_z=state(spin_system,'Lz','1H');
coil_y=state(spin_system,'Ly','1H');

% RF field operator
Lx=operator(spin_system,'Lx','1H');

% Drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Bimodal RF field distribution in rad/s
b1_freq=2*pi*linspace(25e3,65e3,201).';
b1_dist=0.8*exp(-(b1_freq-2*pi*50e3).^2/(2*(2*pi*3.0e3)^2))+...
        0.2*exp(-(b1_freq-2*pi*38e3).^2/(2*(2*pi*2.5e3)^2));
b1_dist=b1_dist/trapz(b1_freq,b1_dist);

% Probability masses on the RF field grid
b1_mass=b1_dist*(b1_freq(2)-b1_freq(1));

% Nutation curve time grid
dt=2e-6; npts=256;

% Nutation curve accumulator
curve=zeros(npts,1);

% Loop over the RF field ensemble
for n=1:numel(b1_freq)

    % Magnetisation trajectories for the current RF field
    traj=evolution(spin_system,H+b1_freq(n)*Lx,[coil_z coil_y],...
                   rho,dt,npts-1,'multichannel');

    % Weighted accumulation of the complex nutation curve
    curve=curve+b1_mass(n)*(traj(1,:)+1i*traj(2,:)).';

end

% Complex measurement noise
rng(1); curve=curve+2e-3*(randn(npts,1)+1i*randn(npts,1));

% Nutation frequency distribution recovery
[freq,distr]=nutation_dist(curve,dt);

% Plot the source and recovered distributions
kfigure(); plot(b1_freq/(2*pi*1e3),2*pi*1e3*b1_dist,'b-');
hold on; plot(freq/(2*pi*1e3),2*pi*1e3*distr,'r.');
kxlabel('nutation frequency, kHz'); kgrid();
kylabel('probability density, kHz$^{-1}$');
klegend({'source distribution','recovered distribution'},...
        'Location','northwest');

end


