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

% Build the rotating-frame single-spin model
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
Sz=state(spin_system,'Lz',1);
Sz=Sz/norm(full(Sz),2);
Sx=state(spin_system,'Lx',1);
Sx=Sx/norm(full(Sx),2);
Lx=operator(spin_system,'Lx',1);
Ly=operator(spin_system,'Ly',1);
Lz=operator(spin_system,'Lz',1);
H=hamiltonian(assume(spin_system,'nmr'));

% Propagate each offset and physical radius independently
mag_x=zeros(numel(radii_m),numel(offset_hz));
rf_scale=6e-3./radii_m;
parfor n_offset=1:numel(offset_hz)
    drift=H+2*pi*offset_hz(n_offset)*Lz;
    column=zeros(numel(radii_m),1);
    for n_radius=1:numel(radii_m)
        rf_x=2*pi*rf_scale(n_radius)*rf_hz(1,:);
        rf_y=2*pi*rf_scale(n_radius)*rf_hz(2,:);
        final_state=shaped_pulse_xy(spin_system,drift,{Lx,Ly},...
                    {rf_x,rf_y},pulse_dt,Sz,'expv-pwc');
        column(n_radius)=real(Sx'*final_state);
    end
    mag_x(:,n_offset)=column;
end

% Apply the paper's equal-radius detection integral
response.radii_m=radii_m;
response.offset_hz=offset_hz;
response.mag_x=mag_x;
response.detected=trapz(radii_m,mag_x,1)/(radii_m(end)-radii_m(1));

% Display the spatial-offset surface and detected offset profile
fig=kfigure();
tiledlayout(1,2);
nexttile;
surf(offset_hz/1e3,rf_scale,mag_x,'EdgeColor','none');
xlabel('offset (kHz)'); ylabel('relative B_1');
zlabel('M_x'); title('toroid radial response'); colorbar;
view(35,25);
nexttile;
plot(offset_hz/1e3,response.detected,'LineWidth',1.5);
xlabel('offset (kHz)'); ylabel('detected M_x');
title('radius-weighted signal'); ylim([-1 1]); grid on;

end


