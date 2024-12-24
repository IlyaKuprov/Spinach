% Field dependence of the couping factor between 13C of CHCl3 and the
% electron spin of a nitroxide radical. Further particulars here:
%
%             https://doi.org/10.1016/j.jmro.2022.100040
%
% Experimental data from:
%
%               https://doi.org/10.1002/anie.201811892
%
% Calculation time: seconds
% 
% tomas.orlando@mpinat.mpg.de
% ilya.kuprov@weizmann.ac.il

function scalar_dnp()

% Spin system
sys.isotopes={'E','13C'};

% Coordinates for dipolar Redfield
inter.coordinates{1}=[0.0 0.0 0.0];   % Angstrom
inter.coordinates{2}=[0.0 0.0 3.1];   % Angstrom

% Zeeman interactions for CSA/g-aniso Redfield
inter.zeeman.eigs{1}=[2.0029 2.0065 2.0098]; % Bohr magneton units
inter.zeeman.eigs{2}=[100 120 120];          % ppm
inter.zeeman.euler{1}=[0 0 0];               % radians
inter.zeeman.euler{2}=[0 0 0];               % radians

% Static hyperfine coupling
inter.coupling.scalar=cell(2,2);
inter.coupling.scalar{1,2}=2e6;  % Hz

% Formalism and approximation
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theories
inter.relaxation={'SRFK','redfield','t1_t2'};

% Electron R1 and R2 for empirical T1/T2
inter.r1_rates{1}=2e6;  % Hz, T1 of 500 ns
inter.r2_rates{1}=5e6;  % Hz, T2 of 200 ns

% Nuclear R1 and R2 for empirical T1/T2
inter.r1_rates{2}=0.25;  % Hz, T1 of 4 s
inter.r2_rates{2}=3.00;  % Hz, T2 of 1/3 s

% Correlation time for rotational Redfield
inter.tau_c={30e-12}; % seconds

% Parameters of scalar collisional Redfield
inter.srfk_tau_c={[0.62, 30.0e-12],...  % [weight, tau_c]
                  [0.38, 0.80e-12]};    % [weight, tau_c]            
inter.srfk_mdepth=cell(2,2);
inter.srfk_mdepth{1,2}=3.6e6;           % Hz

% Temperature and thermalisation 
inter.temperature=298;
inter.equilibrium='dibari';

% Secular approximation
inter.rlx_keep='secular';

% Magnetic field axis
b_vector=[logspace(-2,1,30) 15 20 25 30];

% Prevent excessive output
sys.disable={'hygiene'}; sys.output='hush';

% Loop over magnet fields
Rx=zeros(1,34); R1n=zeros(1,34);
parfor n=1:numel(b_vector)

    % Localise and set the magnet field
    locsys=sys; locsys.magnet=b_vector(n);

    % Run Spinach housekeeping
    spin_system = create(locsys,inter);
    spin_system = basis(spin_system,bas);

    % Get the relaxation superoperator
    R=relaxation(spin_system);

    % Get and normalise the states of interest
    Nz=state(spin_system,'Lz','13C'); Nz = Nz/norm(Nz,2);
    Ez=state(spin_system,'Lz','E');   Ez = Ez/norm(Ez,2);

    % Extract the relaxation rates
    R1n(n)=real(Nz'*R*Nz); % nuclear relaxation rate
    Rx(n)=real(Nz'*R*Ez);  % cross-relaxation rate

end

% Compute DNP coupling factor
sigma_over_R1 = Rx./R1n;

% Experimental data
expt_data=[0.34, -0.84; 
           1.20, -0.48;
           3.40, -0.37;
           9.40, -0.16];

% Plotting
figure(); hold on;
plot(b_vector,sigma_over_R1,'b-'); 
plot(expt_data(:,1),expt_data(:,2),'ro');
set(gca,'XScale','log'); box on;
kxlabel('Magnetic field / Tesla');
kylabel('NOE coupling factor Rx/R1n'); 
kgrid; axis([1e-2 30 -1 0]);
legend({'simulation','experiment'},...
       'Location','NorthWest');

end

