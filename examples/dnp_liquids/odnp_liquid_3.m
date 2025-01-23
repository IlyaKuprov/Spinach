% Steady state nuclear magnetisation as a function of microwave frequency
% offset and the magnet field in a DNP experiment with an electron and a 
% nucleus connected by a hyperfine coupling. A g-hyperfine cross-correla-
% tion effect is visible at high field.
%
% The steady state is computed by setting the time derivative to zero in
% the inhomogeneous master equation, and solving the resulting algebraic
% equation for the steady state density matrix.
%
% Calculation time: minutes.
%
% mariagrazia.concilio@sjtu.edu.cn 
% ilya.kuprov@weizmann.ac.il

function odnp_liquid_3()

% Spin system
sys.isotopes={'1H','E'};

% Anisotropic Zeeman interactions
inter.zeeman.eigs{1}=[15 5 -20];                % ppm
inter.zeeman.eigs{2}=[2.00210 2.00250 2.00290]; % Bohr magneton units
inter.zeeman.euler{1}=[0.00 0.00 0.00];
inter.zeeman.euler{2}=[pi/3 pi/4 pi/5];

% Isotropic hyperfine coupling            
inter.coupling.scalar{1,2}=20e6;    % Hz
inter.coupling.scalar{2,2}=0;  

% Coordinates for dipolar coupling
inter.coordinates={[0.0 0.0 0.0]
                   [0.0 0.0 3.0]};  % Angstrom
               
% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory
inter.relaxation={'redfield'};
inter.equilibrium='zero';       % Required for steady state
inter.rlx_keep='secular';
inter.temperature=298;
inter.tau_c={10e-12};           % TEMPOL in water
sys.tols.rlx_integration=1e-10; % Needs to be this tight

% Sequence parameters
parameters.spins={'E'};
parameters.mw_pwr=2*pi*500e3;
parameters.method='lvn-backs';
parameters.needs={'rho_eq'};
parameters.g_ref=mean(inter.zeeman.eigs{2});

% Turn off excessive output
sys.disable={'hygiene'}; sys.output='hush';

% Microwave frequency offset grid
parameters.mw_frq=2*pi*linspace(-15,15,512)*1e6;

% Magnetic field grid, Tesla
field_grid=linspace(1,10,64);

% Preallocate result
answer=zeros(numel(parameters.mw_frq),numel(field_grid));

% Loop over magnet fields
parfor n=1:numel(field_grid)
    
    % Localise parameter arrays
    localsys=sys; locpar=parameters;

    % Set the magnet field
    localsys.magnet=field_grid(n);

    % Spinach housekeeping
    spin_system=create(localsys,inter);
    spin_system=basis(spin_system,bas);
    
    % Relevant operators and states
    locpar.coil=state(spin_system,'Lz','1H');
    locpar.mw_oper=operator(spin_system,'Lx','E');
    locpar.ez_oper=operator(spin_system,'Lz','E');

    % Steady state simulation
    answer(:,n)=liquid(spin_system,@dnp_freq_scan,locpar,'esr');
    
end

% Plotting
figure(); imagesc(field_grid,(parameters.mw_frq/(2*pi*1e6)),real(answer)); 
xlim([min(field_grid) max(field_grid)]); kxlabel('$B_0$, Tesla');
kylabel('MW freq offset from $g_{iso}$, MHz');
kcolourbar('$\langle H_{\rm{Z}} \rangle$ at the steady state');  

end

