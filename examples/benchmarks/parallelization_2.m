% Parallelization test: multi-threaded evaluation of observables in Hilbert
% space time propagation for pyrene radical spin system at low field. For
% further information, see:
%
%                   http://dx.doi.org/10.1063/1.3679656
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk

function parallelization_2()

% Read the spin system properties (vacuum DFT calculation)
options.no_xyz=1;
[sys,inter]=g2spinach(gparse('../standard_systems/pyrene_cation.log'),...
                                {{'E','E'},{'H','1H'}},[0 0],options);
% Magnet field
sys.magnet=50e-6;

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Assumptions
spin_system=assume(spin_system,'labframe');

% Hamiltonian operator
[H,Q]=hamiltonian(spin_system);
H=H+orientation(Q,[pi/3,pi/4,pi/5]);

% Initial state
rho=operator(spin_system,'Lz','E');

% Parallel propagation benchmark, 200 steps
ncores=[2 4 8 16 32 64 128 256 512 1024];
ncores=ncores(ncores<=feature('numcores'));
for n=ncores
    delete(gcp('nocreate')); parpool(n); pause(10); tic;
    evolution(spin_system,H,rho,rho,5e-9,200,'observable');
    disp([num2str(n) ' cores, execution time ' num2str(toc) ' seconds']);
end

end

