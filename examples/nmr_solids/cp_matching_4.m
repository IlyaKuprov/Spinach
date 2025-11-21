% Hartmann-Hahn matching condition test for a cross-polarisation
% experiment between a proton and a 15N nucleus in the presence of
% conformational exchange between two geometries that differ by
% 90 degrees in the N-H vector direction.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il
% mariagrazia.concilio@sjtu.edu.cn

function cp_matching_4()

% System specification
sys.magnet=9.394;
sys.isotopes={'1H','15N','1H','15N'};
          
% Interactions
inter.zeeman.scalar={0.0 0.0 0.0 0.0};
inter.coordinates={[0 0 0]; [0 0 2];
                   [0 0 0]; [0 2 0]};
               
% Chemical exchange
inter.chem.parts={[1 2],[3 4]};
inter.chem.rates=[-5000 +5000;
                  +5000 -5000];
inter.chem.concs=[1 1];

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Relevant operators
Nx=operator(spin_system,'Lx','15N');
Ny=operator(spin_system,'Ly','15N');
Hx=operator(spin_system,'Lx','1H');
Hy=operator(spin_system,'Ly','1H');

% Power levels
powers=linspace(20e3,80e3,60);

% Experiment parameters
parameters.rate=10000;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.max_rank=3;
parameters.spins={'1H','15N'};
parameters.irr_opers={Hx Nx};
parameters.exc_opers={0*Hy 0*Ny};
parameters.rho0=state(spin_system,'Lx','1H');
parameters.coil=state(spin_system,'Lx','15N');
parameters.time_steps=4e-5*ones(1,10);
parameters.grid='rep_2ang_200pts_oct';

% Parallel loop over power levels
parfor n=1:60
    
    % MAS parameters
    localpar=parameters;
    localpar.irr_powers=[powers(n)*ones(1,10);
                              50e3*ones(1,10)];
    
    % Simulation
    fid=singlerot(spin_system,@cp_contact_hard,localpar,'nmr');
    cp(n)=real(fid(end));
    
end

% Plotting
kfigure(); plot(powers'/1e3,cp);
kylabel('$^{15}$N signal, a.u.'); kgrid;
kxlabel('$^{1}$H spin-lock RF power, kHz');

end

