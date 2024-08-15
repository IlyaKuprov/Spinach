% Cross-formalism test of state recovery towards the 
% thermodynamic equilibrium.
%
% i.kuprov@soton.ac.uk

function thermal_equilibrium_5()

% Magnet field
sys.magnet=9.4;

% Isotopes
sys.isotopes={'19F','19F','19F','19F'};

% Chemical shifts
inter.zeeman.scalar={-120.5380 -133.9429 -129.3169 -129.5320};

% J-couplings
inter.coupling.scalar=cell(4,4);
inter.coupling.scalar{1,2}=  271.2924;
inter.coupling.scalar{3,4}=  271.2924;
inter.coupling.scalar{1,3}=    0.5401;
inter.coupling.scalar{1,4}=  -25.9884;
inter.coupling.scalar{2,3}=    9.9625;
inter.coupling.scalar{2,4}=  -40.7675;       

% Relaxation theory parameters
inter.relaxation={'damp'};
inter.temperature=40;
inter.damp_rate=5.0;
inter.rlx_keep='diagonal';

% Formalisms and methods to test
formalisms={'sphten-liouv','zeeman-liouv'};
methods={'dibari','IME'}; traj=cell(2,2);

% Loop over formalisms
for n=1:numel(formalisms)

    % Loop over methods
    for k=1:numel(methods)

        % Thermalisation method
        inter.equilibrium=methods{k};

        % Basis set
        bas.formalism=formalisms{n};
        bas.approximation='none';

        % Spinach housekeeping
        spin_system=create(sys,inter);
        spin_system=basis(spin_system,bas);

        % Get thermal equilibrium state
        rho_eq=equilibrium(spin_system);

        % Flip the magnetisation over
        Fx=operator(spin_system,'Lx','19F');
        rho=step(spin_system,Fx,rho_eq,pi);

        % Get the evolution generator
        spin_system=assume(spin_system,'nmr');
        L=hamiltonian(spin_system)+1i*relaxation(spin_system);

        % Get recovery trajectories
        coil=state(spin_system,'Lz','19F');
        traj{n,k}=evolution(spin_system,L,coil,rho,1e-3,1000,'observable');

    end

end

% Differences between formalisms and methods
a=norm(traj{1,1}-traj{2,1},2)/norm(traj{1,1}+traj{2,1},2); 
b=norm(traj{1,2}-traj{2,2},2)/norm(traj{1,2}+traj{2,2},2);
c=norm(traj{1,1}-traj{1,2},2)/norm(traj{1,1}+traj{1,2},2); 
d=norm(traj{2,1}-traj{2,2},2)/norm(traj{2,1}+traj{2,2},2);

% Report the test result
if max([a b c d])>1e-3
    error('Thermalisation formalism-method cross-check FAILED.');
else
    disp('Thermalisation formalism-method cross-check PASSED.');
end

end

