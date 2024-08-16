% Test of the invariance of the thermal equilibrium state under the
% thermalised relaxation superoperator.
%
% i.kuprov@soton.ac.uk

function thermal_equilibrium_4()

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
inter.equilibrium='zero';        % Thermalised later
inter.temperature=40;
inter.damp_rate=5.0;
inter.rlx_keep='diagonal';

% Formalisms to test
formalisms={'sphten-liouv','zeeman-liouv'};

% Loop over formalisms
for n=1:numel(formalisms)

    % Basis set
    bas.formalism=formalisms{n};
    bas.approximation='none';

    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

    % Get thermal equilibrium state
    rho_eq=equilibrium(spin_system);

    % Get relaxation superoperator
    R=relaxation(spin_system);

    % Thermalise using Dibari-Levitt method
    HLSPS=hamiltonian(assume(spin_system,'labframe'),'left');
    Rt=thermalize(spin_system,R,HLSPS,inter.temperature,rho_eq,'dibari');

    % Check the norm of the action
    if norm(Rt*rho_eq,2)>1e-9, error('R*rho_eq action test FAILED.'); end

    % Thermalise using inhomogeneous master equation
    Rt=thermalize(spin_system,R,[],[],rho_eq,'IME');

    % Check the norm of the action
    if norm(Rt*rho_eq,2)>1e-9, error('R*rho_eq action test FAILED.'); end

end

% Report good test result
disp('R*rho_eq action test PASSED.');

end

