% Test of the thermal equilibrium functionality against the 
% textbook expressions for the Boltzmann populations.
%
% i.kuprov@soton.ac.uk

function thermal_equilibrium_3()

% X-band magnet
sys.magnet=0.34;

% Electron and two protons
sys.isotopes={'E','1H','1H'};

% Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5],[0 5 0]};
inter.zeeman.euler=(pi/180)*{[0 10 0],[0 0 10],[100 0 0]};

% Cartesian coordinates
inter.coordinates={[0.000 0.000 0.000];
                   [0.000 3.500 0.000];
                   [2.475 2.475 0.000]};

% Spin temperature
inter.temperature=80;

% Formalisms to test
formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'};

% Loop over formalisms
for n=1:numel(formalisms)

    % Formalism and basis set
    bas.formalism=formalisms{n};
    bas.approximation='none';

    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

    % Isotropic thermal equilibrium
    rho_eq=equilibrium(spin_system);

    % Detection states
    coil_a=state(spin_system,{'Lz'},{1});
    coil_b=state(spin_system,{'Lz'},{2});
    coil_c=state(spin_system,{'Lz'},{3});

    % Get the expectation values from Spinach
    expt_a_spinach=trace(coil_a'*rho_eq);
    expt_b_spinach=trace(coil_b'*rho_eq);
    expt_c_spinach=trace(coil_c'*rho_eq);

    % Get the expectation values from textbook
    [~,PE]=levelpop('E',sys.magnet,inter.temperature);
    [~,PH]=levelpop('1H',sys.magnet,inter.temperature);
    expt_a_textbook=0.5*PE(1)-0.5*PE(2);
    expt_b_textbook=0.5*PH(1)-0.5*PH(2);
    expt_c_textbook=0.5*PH(1)-0.5*PH(2);

    % Test for differences
    if (abs((expt_a_spinach-expt_a_textbook)/expt_a_textbook)>1e-3)||...
       (abs((expt_b_spinach-expt_b_textbook)/expt_b_textbook)>1e-3)||...
       (abs((expt_c_spinach-expt_c_textbook)/expt_c_textbook)>1e-3)
        error('Thermodynamic equilibrium test FAILED.');
    end

end

% Report success
disp('Thermodynamic equilibrium tests PASSED.');

end

