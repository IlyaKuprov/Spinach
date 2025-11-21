% Cross-polarisation contact curve under magic angle spinning 
% in the presence of chemical exchange for H1, H4 and C19 in
% the unit cell of monohydrocalcite. Further details in:
%
%         https://doi.org/10.1038/s41467-023-44381-x
%
% Calculation time: hours, much faster on a GPU.
%
% guinevere.mathies@uni-konstanz.de

function cp_mas_powder_mhc_fplanck_exchange()

% 400 MHz NMR
sys.magnet=9.4;

% Read CASTEP file
props=c2spinach('mhc.magres');

% Drop O and Ca atoms
drop_mask=ismember(props.symbols,{'O','Ca'});
props.symbols(drop_mask)=[];
props.std_geom(drop_mask,:)=[];
props.cst(drop_mask)=[];

% Two chemical endpoints: H1, H4, and C19,
% with H1 and H4 under chemical exchange
sys.isotopes{1}='1H';  sys.isotopes{4}='1H';
sys.isotopes{2}='1H';  sys.isotopes{5}='1H';
sys.isotopes{3}='13C'; sys.isotopes{6}='13C';

% Convert shielding tensors into shift using the
% parametrisation of Huang et al. ACIE 2021
inter.zeeman.matrix{1}=29.25*eye(3)-props.cst{1};
inter.zeeman.matrix{2}=29.25*eye(3)-props.cst{4};
inter.zeeman.matrix{3}=169.86*eye(3)-props.cst{19};
inter.zeeman.matrix{4}=29.25*eye(3)-props.cst{4};
inter.zeeman.matrix{5}=29.25*eye(3)-props.cst{1};
inter.zeeman.matrix{6}=169.86*eye(3)-props.cst{19};

% Cartesian coordinates
inter.coordinates={props.std_geom(1,:);
                   props.std_geom(4,:);
                   props.std_geom(19,:);
                   props.std_geom(4,:);
                   props.std_geom(1,:);
                   props.std_geom(19,:)};

% Chemical kinetics endpoints
inter.chem.parts={[1 2 3],[4 5 6]};

% Initial concentrations (arb. units)
inter.chem.concs=[1 1];

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Disable start-up checks
sys.disable={'hygiene'};

% Enable GPU
% sys.enable={'gpu'};

% Exchange rate constant array
exch_rates=[1e1 1e2 1e3 1e4 1e5 1e6]; % Hz

% Experiment setup
parameters.spins={'1H','13C'};
parameters.rate=10000;
parameters.axis=[1 1 1];
parameters.max_rank=7;
parameters.grid='rep_2ang_800pts_sph';
parameters.offset=[2e3 1e4];      % 5 ppm 1H, 100 ppm 13C
parameters.hi_pwr=83e3;           % Hz 
parameters.cp_pwr=[60e3 50e3];    % Hz
parameters.nsteps=1000; 
parameters.timestep=1e-5;         % sec
parameters.needs={'iso_eq'};
parameters.verbose=1;

% Preallocate contact curve array
contact_curves=zeros(numel(exch_rates),...
                     parameters.nsteps+1);

% Loop over exchange rates
for n=1:numel(exch_rates)

    % Set the exchange rates
    inter.chem.rates=exch_rates(n)*[-1  1;
                                     1 -1];

    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
    
    % Detection state
    parameters.coil=state(spin_system,'L+','13C');

    % Simulation
    contact_curves(n,:)=singlerot(spin_system,@cp_contact_soft,parameters,'nmr');

end

% Generate the time axis
time_axis=linspace(0,parameters.timestep*parameters.nsteps,parameters.nsteps+1);

% Do the plotting
kfigure(); plot(1e3*time_axis,real(contact_curves));
kxlabel('contact time, ms'); grid on;
kylabel('$^{13}$C NMR signal, a.u.');
legend('10 Hz','100 Hz','1 kHz','10 kHz',...
       '100 kHz','1 MHz','Location','northeast');

end

