% Powder-averaged stimulated-echo ESEEM on a model electron-proton
% pair, comparing ideal and 25 ns rectangular 90-degree pulses after
% transverse preparation. The timing and interaction parameters follow
% Goodwin et al., https://doi.org/10.1016/j.jmro.2026.100218. Spinach
% uses its own spherical quadrature and a uniform electron offset
% distribution spanning 30 MHz; the receiver phase is zero (Lx).
%
% Calculation time: about one minute with four local workers
%
% Talos <talos@spindynamics.org>

function eseem_3p_powder()

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

% Specify preparation, receiver, and timing
parameters.spins={'E'};
parameters.rho0=state(spin_system,'Lx','E');
parameters.coil=coil_state(spin_system,'Lx','E');
parameters.pulse_op=operator(spin_system,'Lx','E');
parameters.tau=104e-9;
parameters.timestep=4e-9;
parameters.npoints=256;
parameters.grid='rep_2ang_100pts_sph';

% Average the complex echoes over a uniform electron offset distribution
offsets=-15e6:0.5e6:15e6;
durations=[0 25e-9];
fid=zeros(parameters.npoints,numel(durations));
for n=1:numel(durations)

    % Keep the flip angle fixed while changing the pulse duration
    parameters.pulse_duration=durations(n);
    for offset=offsets

        % Integrate over static orientations at this electron offset
        parameters.offset=offset;
        fid(:,n)=fid(:,n)+powder(spin_system,@eseem_3p,parameters,'esr');

    end

end
fid=fid/numel(offsets);

% Plot the absorptive quadrature on the storage-delay axis
time_axis=(0:(parameters.npoints-1))*parameters.timestep*1e6;
kfigure();
plot(time_axis,real(fid));
kxlabel('storage delay, $\mu$s');
kylabel('in-phase echo amplitude');
klegend({'Hard pulses','25 ns pulses'});
kgrid;

end


