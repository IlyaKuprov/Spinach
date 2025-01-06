% A spin-locking experiment on a two-spin system.
%
% ilya.kuprov@weizmann.ac.il

function spin_lock()

% Isotopes
sys.isotopes={'1H','1H'};

% Magnetic induction
sys.magnet=5.9;

% Chemical shifts
inter.zeeman.scalar={1.0 1.5};

% Scalar couplings
inter.coupling.scalar=cell(2,2);
inter.coupling.scalar{1,2}=7.0; 

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Initial state
rho=4*state(spin_system,'Lz','all');

% Observable states
coil_x1=state(spin_system,{'Lx'},{1});
coil_x2=state(spin_system,{'Lx'},{2});
coil_y1=state(spin_system,{'Ly'},{1});
coil_y2=state(spin_system,{'Ly'},{2});
coil_z1=state(spin_system,{'Lz'},{1});
coil_z2=state(spin_system,{'Lz'},{2});

% Pulse operators
Lx=operator(spin_system,'Lx','all');
Ly=operator(spin_system,'Ly','all');

% Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Spin-locking field of 1.5 kHz along Y
H=H+2*pi*1.5e3*Ly;

% Initial 90-degree pulse
rho=step(spin_system,Lx,rho,pi/2);

% Time evolution
answer=evolution(spin_system,H,[coil_x1 coil_y1 coil_z1 ...
                                coil_x2 coil_y2 coil_z2],...
                 rho,1e-4,100,'multichannel');

% Bloch sphere plot
[X,Y,Z]=sphere; figure();
surf(X,Y,Z,'FaceAlpha',0.5,'EdgeAlpha',0.1);
answer=real(answer)'; hold on; colormap bone;
plot3(answer(:,1),answer(:,2),answer(:,3),'b-');
plot3(answer(:,4),answer(:,5),answer(:,6),'r-');
axis([-1 1 -1 1 -1 1]); axis square;

end

