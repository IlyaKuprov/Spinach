% Parallelization test: multi-threaded evaluation of observables in 
% Hilbert space time propagation.
%
% For further information, see: http://dx.doi.org/10.1063/1.3679656
% 
% Spin system of 3-phenylmethylene-1H,3H-naphtho-[1,8-c,d]-pyran-1-one.
% Source: Penchav, et al., Spec. Acta Part A, 78 (2011) 559-565.
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk

function parallelization_1()

% Magnetic induction
sys.magnet=14.095;

% Spin system
sys.isotopes={'1H','1H','1H','1H','1H','1H',...
              '1H','1H','1H','1H','1H','1H'};

% Chemical shifts
inter.zeeman.scalar={8.345,7.741,8.097,8.354,7.784,8.330,...
                     7.059,7.941,7.466,7.326,7.466,7.941};

% Scalar couplings
inter.coupling.scalar=cell(12,12);
inter.coupling.scalar{1,2}=7.8;
inter.coupling.scalar{1,3}=0.9;
inter.coupling.scalar{2,3}=7.8;
inter.coupling.scalar{4,5}=8.4;
inter.coupling.scalar{4,6}=1.2;
inter.coupling.scalar{5,6}=7.2;
inter.coupling.scalar{8,9}=7.8;
inter.coupling.scalar{8,10}=1.2;
inter.coupling.scalar{9,10}=7.8;
inter.coupling.scalar{10,11}=7.8;
inter.coupling.scalar{10,12}=1.2;
inter.coupling.scalar{11,12}=7.8;

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Assumptions
spin_system=assume(spin_system,'nmr');

% Hamiltonian operator
H=hamiltonian(spin_system);

% Initial state
rho=operator(spin_system,'Lx','all');

% Parallel propagation benchmark, 1000 steps
ncores=[2 4 8 16 32 64 128 256 512 1024];
ncores=ncores(ncores<=feature('numcores'));
for n=ncores
    delete(gcp('nocreate')); parpool(n); pause(10); tic;
    evolution(spin_system,H,rho,rho,1e-3,1000,'observable');
    disp([num2str(n) ' cores, execution time ' num2str(toc) ' seconds']);
end

end

