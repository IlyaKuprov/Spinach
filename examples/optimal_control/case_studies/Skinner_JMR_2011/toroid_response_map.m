% Simulates toroid hard and shaped pulses under Skinner et al.
% Syntax:
%
%    [response,fig]=toroid_response_map(rf_hz,pulse_dt,...
%                                         offset_hz,radii_m)
%
% Parameters:
%
%    rf_hz      - two-row Cartesian waveform in Hz at the outer
%                 radius; default is a 25 kHz y-phase hard pulse
%
%    pulse_dt   - slice durations in seconds; default is 3.74 us
%
%    offset_hz  - resonance offsets in Hz; default is +/-1.5 kHz
%
%    radii_m    - sample radii in metres; default is 1-6 mm
%
% Outputs:
%
%    response   - x magnetization and radial average vs offset
%
%    fig        - offset-radius response and detected profile
%
% Source: Skinner et al., J. Magn. Reson. 209, 282-290 (2011).
% DOI: 10.1016/j.jmr.2011.01.026, Eqs. (1),(3), Fig. 8A.
% This default is the paper's optimised rectangular benchmark,
% not its numerically optimised shaped pulse in Fig. 8B.
%
function [response,fig]=toroid_response_map(rf_hz,pulse_dt,...
                                               offset_hz,radii_m)

% Use the paper's probe dimensions and hard-pulse benchmark
if nargin<1
    rf_hz=[0;25e3]; pulse_label='hard-pulse';
else
    pulse_label='shaped-pulse';
end
if nargin<2, pulse_dt=3.74e-6; end
if nargin<3, offset_hz=linspace(-1.5e3,1.5e3,51); end
if nargin<4, radii_m=linspace(1e-3,6e-3,51); end

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
zlabel('M_x'); title([pulse_label ' toroid response']); colorbar;
view(35,25);
nexttile;
plot(offset_hz/1e3,response.detected,'LineWidth',1.5);
xlabel('offset (kHz)'); ylabel('detected M_x');
title([pulse_label ' detected signal']); ylim([-1 1]); grid on;

end


