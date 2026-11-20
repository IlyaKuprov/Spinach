% Static powder-averaged cRW-OPT transfer from electron Lx to proton -Lz:
% on-resonance build-up and a 13-block electron offset profile. The
% waveform and model follow Goodwin et al., J. Magn. Reson. Open (2026),
% https://doi.org/10.1016/j.jmro.2026.100218. The powder and offset grids
% are reduced for this demonstration. Build-up is sampled at block
% boundaries across the original 5 microsecond window, without relaxation.
%
% Talos <talos@spindynamics.org>

function crw_opt_profile()

% Set the proton Larmor frequency and tilted axial hyperfine tensor
sys.isotopes={'E','1H'};
sys.magnet=2*pi*14.8e6/spin('1H');
axis=[1; 0; 1]/sqrt(2);
inter.coupling.matrix={[],1e6*(3*(axis*axis')-eye(3)); [],[]};

% Use the complete two-spin basis
bas.formalism='zeeman-liouv';
bas.approximation={'none'};
sys.parallel={'processes',4};
sys.output='hush';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set the original nutation frequencies and phases on 5 ns slices
nutation_hz=[23271500 31999700 31999700 31999800 27118800 31981000 28153700 24655300 ...
             30551100 31996600 31999900 7908330 31998400 31999900 31929500 31155000 ...
             20661300 29851800 24224600 31998500 32000000 31999700 18769000 31936100 ...
             23104800 31988900 513988 32000000 31988300 26226500];
phase_deg=[0 180 0 0 180 180 180 180 ...
           0 180 180 180 0 0 0 0 ...
           0 180 180 0 0 0 0 180 ...
           0 0 0 180 180 180];
parameters.slice_durs=5e-9*ones(size(nutation_hz));
omega=2*pi*nutation_hz;
parameters.amplitudes={omega.*cosd(phase_deg),omega.*sind(phase_deg)};
parameters.pulse_ops={operator(spin_system,'Lx','E'),...
                      operator(spin_system,'Ly','E')};

% Prepare the electron and detect negative proton longitudinal magnetisation
parameters.rho0=state(spin_system,'Lx','E');
parameters.coil=-coil_state(spin_system,'Lz','1H');
parameters.spins={'E'};
parameters.grid='rep_2ang_400pts_sph';

% Calculate the on-resonance build-up at successive block boundaries
parameters.offset=0;
parameters.nloops=33;
buildup=powder(spin_system,@dnp_train,parameters,'esr');
time_axis=(0:parameters.nloops)*sum(parameters.slice_durs)*1e6;

% Scan the signed electron detuning after 13 periods
offsets=(-60:2:60)*1e6;
parameters.nloops=13;
parameters.grid='rep_2ang_100pts_sph';
profile=zeros(size(offsets),'like',1i);
parfor n=1:numel(offsets)

    % Average each offset over the same powder grid
    localpar=parameters;
    localpar.offset=offsets(n);
    signal=powder(spin_system,@dnp_train,localpar,'esr');
    profile(n)=signal(end);

end

% Plot the signed transfer without peak normalisation
kfigure();
subplot(2,1,1);
plot(time_axis,real(buildup)); kgrid;
kxlabel('irradiation time, $\mu$s');
kylabel('$-\langle I_z\rangle$');

subplot(2,1,2);
plot(offsets/1e6,real(profile)); kgrid;
kxlabel('electron detuning, MHz');
kylabel('$-\langle I_z\rangle$');

end


