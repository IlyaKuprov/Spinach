% Simulates toroid hard and shaped pulses under Skinner et al.
% Syntax:
%
%    [response,fig]=toroid_response_map(rf_hz,pulse_dt,...
%                                         offset_hz,radii_m)
%
% Parameters:
%
%    rf_hz      - two-row Cartesian waveform in Hz at the outer
%                 radius; [0;25e3] gives the Fig. 8A hard pulse
%
%    pulse_dt   - slice durations in seconds; 3.74e-6 for Fig. 8A
%
%    offset_hz  - resonance offsets in Hz; Fig. 8A uses +/-1.5 kHz
%
%    radii_m    - sample radii in metres spanning 1-6 mm
%
% Outputs:
%
%    response   - x magnetisation by radius and offset, plus
%                 the radial integral at each offset
%
%    fig        - offset-radius response and detected profile
%
% Source: Skinner et al., J. Magn. Reson. 209, 282-290 (2011).
% DOI: 10.1016/j.jmr.2011.01.026, Eqs. (1),(3), Fig. 8A.
% With [0;25e3] and 3.74e-6, this simulates the rectangular
% Fig. 8A benchmark, not the shaped pulse in Fig. 8B.
%
function [response,fig]=toroid_response_map(rf_hz,pulse_dt,...
                                               offset_hz,radii_m)

% Magnetic field and isotope
sys.magnet=14.1;
sys.isotopes={'1H'};

% Chemical shift, ppm
inter.zeeman.scalar={0};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation={'none'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Normalised initial and target states
rho_z=state(spin_system,'Lz',1);
rho_z=rho_z/norm(full(rho_z),2);
rho_x=state(spin_system,'Lx',1);
rho_x=rho_x/norm(full(rho_x),2);

% Control and offset operators
lx=operator(spin_system,'Lx',1);
ly=operator(spin_system,'Ly',1);
lz=operator(spin_system,'Lz',1);

% Drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Preallocate the radius-by-offset response
mag_x=zeros(numel(radii_m),numel(offset_hz));

% Inverse-radius RF scaling
rf_scale=6e-3./radii_m;

% Loop over resonance offsets
parfor n_offset=1:numel(offset_hz)

    % Offset Hamiltonian and response column
    drift=H+2*pi*offset_hz(n_offset)*lz;
    column=zeros(numel(radii_m),1);

    % Loop over sample radii
    for n_radius=1:numel(radii_m)

        % Scale the waveform at this radius
        rf_x=2*pi*rf_scale(n_radius)*rf_hz(1,:); %#ok<PFBNS>
        rf_y=2*pi*rf_scale(n_radius)*rf_hz(2,:);

        % Propagate the initial state
        final_state=shaped_pulse_xy(spin_system,drift,{lx,ly},...
                                    {rf_x,rf_y},pulse_dt,rho_z,'expv-pwc');

        % Detect transverse x magnetisation
        column(n_radius)=real(rho_x'*final_state);

    end

    % Store the response column
    mag_x(:,n_offset)=column;

end

% Package the radius-by-offset response
response.radii_m=radii_m;
response.offset_hz=offset_hz;
response.mag_x=mag_x;

% Apply the equal-radius detection integral
response.detected=trapz(radii_m,mag_x,1)/(radii_m(end)-radii_m(1));

% Start the response figure
fig=kfigure();
tiledlayout(1,2);

% Plot the spatial-offset surface
nexttile;
surf(offset_hz/1e3,rf_scale,mag_x,'EdgeColor','none');
kxlabel('offset (kHz)'); kylabel('relative $B_1$');
kzlabel('$M_x$'); ktitle('toroid radial response'); colorbar;
view(35,25);

% Plot the detected offset profile
nexttile;
plot(offset_hz/1e3,response.detected,'LineWidth',1.5);
kxlabel('offset (kHz)'); kylabel('detected $M_x$');
ktitle('radius-weighted signal'); ylim([-1 1]); kgrid;

end


