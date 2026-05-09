% Example setup of a 19F Mims ENDOR powder simulation at 34 GHz. Syntax:
%
%      endor_kehl_mims_powder()
%
% Parameters:
%
%   None.
%
% Outputs:
%
%   none             - this example plots the simulated ENDOR spectrum.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=endor_kehl_mims_powder.m>

function endor_kehl_mims_powder()

    % Spinach spin-system constructor input
    T=0.066;
    D=6.5;

    sys.magnet=1.20521;
    sys.isotopes={'E','19F','1H','1H','14N'};
    sys.output='hush';
    sys.disable={'hygiene'};

    inter.zeeman.matrix=cell(1,5);
    inter.zeeman.matrix{1}=diag([2.00886,2.00610,2.00211]);
    inter.zeeman.matrix{2}=zeros(3,3);
    inter.zeeman.eigs=cell(1,5);
    inter.zeeman.euler=cell(1,5);
    inter.zeeman.eigs{2}=[205,310,358];
    inter.zeeman.euler{2}=[87,-14,-29]*pi/180;
    inter.zeeman.matrix{3}=zeros(3,3);
    inter.zeeman.matrix{4}=zeros(3,3);
    inter.zeeman.matrix{5}=zeros(3,3);

    inter.coupling.eigs=cell(5,5);
    inter.coupling.euler=cell(5,5);
    inter.coupling.eigs{1,2}=[2*T,-T,-T]*1e6;
    inter.coupling.euler{1,2}=[0,-1,-161]*pi/180;
    inter.coupling.eigs{1,5}=[15,11,95.8]*1e6;
    inter.coupling.euler{1,5}=[0,0,0];
    inter.coupling.eigs{5,5}=diag(remtrace(diag([1.2,0.54,-1.7]*1e6)))';
    inter.coupling.euler{5,5}=[0,0,0];
    inter.coupling.eigs{2,3}=[2*D,-D,-D]*1e3;
    inter.coupling.euler{2,3}=[0,-143,16.9]*pi/180;
    inter.coupling.eigs{2,4}=[2*D,-D,-D]*1e3;
    inter.coupling.euler{2,4}=[0,-34.3,47.9]*pi/180;

    inter.relaxation={'t1_t2'};
    inter.r1_rates={1/(4e-3) 0 0 0 0};
    inter.r2_rates={1/(5e-6) 1/(3e-3) 1/(3e-3) 1/(3e-3) 1/(3e-3)};
    inter.equilibrium='zero';
    inter.rlx_keep='diagonal';

    bas.formalism='zeeman-liouv';
    bas.approximation='none';

    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

    % ENDOR context metadata and simulation parameters
    parameters.mw_freq=2.122001324237840e11;
    parameters.static_field=1.20521;
    parameters.field_step=5e-5;
    parameters.endor_res=6283.185307179586;
    parameters.endor_range=1.884955592153876e6;
    parameters.pulse_times=[12e-9,2e-6,12e-9,1e-6,2e-4,1e-6,12e-9];
    parameters.rf_field_from_pulses=true;
    parameters.epr_freq_min=2.010619298297468e11;
    parameters.epr_freq_range=1.884955592153876e10;
    parameters.epr_freq_step=6.283185307179586e6;
    parameters.rf_flip_angle=2.844886680750757;
    parameters.pulse_file=fullfile(fileparts(mfilename('fullpath')),...
        'kehl_mlr09_12ns_pulse.txt');
    parameters.multipulses=false;
    parameters.freqDomain=true;
    parameters.powder=true;
    parameters.Nang=50;
    parameters.Bterm=false;
    parameters.endor_spins=[2,3,4];

    % Actual ENDOR calculation through Spinach-style context and experiment
    [endor_amp,~,x_coords,~]=endor_kehl_context(spin_system,@endor_kehl_mims,parameters,'labframe');

    % Plotting
    sim=endor_amp(:)-endor_amp(1);
    sim=sim/min(sim);
    plot(x_coords(:),sim);

end

