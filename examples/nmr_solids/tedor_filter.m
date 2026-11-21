% Compares the selection used by the Grenoble ZF-TEDOR schedules:
% the 13C-15N element-filter echo and the signed 13C-31P phase cycle.
% The small powder grid and shortened mixing sweep are illustrative.
%
% talos@spindynamics.org

function tedor_filter()

% Field and shared basis
sys.magnet=2*pi*400e6/spin('1H');
sys.parallel={'processes',4};
bas.formalism='sphten-liouv';
bas.approximation={'none'};

% Rotor sampling and pulse calibration
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.max_rank=8;
parameters.grid='leb_2ang_rank_5';
parameters.ncycles=0:12;
parameters.rf=150e3;
parameters.npoints=1;
parameters.offset=[0 0];

% Define the two source spin pairs and signed dipole constants
partners={'15N','31P'};
sequences={'zf-filter','zf-cycle'};
couplings=[35.07 -758.07];
rates=[7400 8000];
curves=zeros(2,numel(parameters.ncycles));
times=4e3*parameters.ncycles./rates.';
inter.coupling.euler={[] deg2rad([110.44 0.90 -110.25]); [] []};

% Run both selections through the same MAS context
for n=1:2

    % Specify the pair and its anisotropic interaction
    sys.isotopes={'13C',partners{n}};
    inter.coupling.eigs={[] couplings(n)*[-1 -1 2]; [] []};

    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

    % Prepare and detect the source carbon observable
    parameters.spins=sys.isotopes;
    parameters.sequence=sequences{n};
    parameters.rate=rates(n);
    parameters.rho0=state(spin_system,'Lx','13C');
    if n==1
        parameters.coil=coil_state(spin_system,'Lx','13C');
    else
        parameters.coil=coil_state(spin_system,'L+','13C');
    end
    curves(n,:)=singlerot(spin_system,@tedor,parameters,'nmr');

end

% Plot the signed real component without normalising the two signals
kfigure(); plot(times(1,:),real(curves(1,:)),'o-',times(2,:),real(curves(2,:)),'s-');
kxlabel('Total variable recoupling time, ms');
kylabel('Real echo amplitude'); klegend({'13C-15N element filter','13C-31P phase cycle'}); kgrid;

end


