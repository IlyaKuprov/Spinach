% Simulates an author-deposited robust universal 180 degree pulse
% Syntax:
%
%    [profiles,fig]=ur180_profiles(offset_hz,rf_scales)
%
% Parameters:
%
%    offset_hz - resonance offsets in Hz; default is 100 points
%                over a total 20 kHz span
%
%    rf_scales - relative RF amplitudes; default is five evenly
%                spaced values from 0.6 to 1.4
%
% Outputs:
%
%    profiles  - state mappings and propagator overlap scores
%
%    fig       - source waveform and robustness plots
%
% Source: Kobzar et al., J. Magn. Reson. 225, 142-160 (2012).
% DOI: 10.1016/j.jmr.2012.09.013. Shape UR180_730u comes from
% https://www.ioc.kit.edu/luy/311.php, archive
% UR180_BW20_0-40pmB1.zip, pm40B1/UR180_730u.
% Original numerical header states 730 us, 0.5 us slices,
% 20 kHz bandwidth, five B1 scales, and 10 kHz maximum RF.
%
function [profiles,fig]=ur180_profiles(offset_hz,rf_scales)

% Use the source waveform's optimised ensemble by default
if nargin<1, offset_hz=linspace(-10e3,10e3,100); end
if nargin<2, rf_scales=linspace(0.6,1.4,5); end

% Build the one-spin rotating-frame Spinach representation
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
Sx=state(spin_system,'Lx',1);
Sx=Sx/norm(full(Sx),2);
Sy=state(spin_system,'Ly',1);
Sy=Sy/norm(full(Sy),2);
Sz=state(spin_system,'Lz',1);
Sz=Sz/norm(full(Sz),2);
Lx=operator(spin_system,'Lx',1);
Ly=operator(spin_system,'Ly',1);
Lz=operator(spin_system,'Lz',1);
H=hamiltonian(assume(spin_system,'nmr'));
initial=[Sx Sy Sz];
target=[-Sx Sy -Sz];

% Read the deposited x/y RF frequencies and slice durations
source_dir=fileparts(mfilename('fullpath'));
waveform=readmatrix(fullfile(source_dir,'ur180_20khz_40pct.dat'));
assert(size(waveform,2)==3);
pulse_dt=waveform(:,3).';
xy_rad=2*pi*waveform(:,1:2).';

% Score all three Cartesian mappings of a y-axis pi rotation
fidelity=zeros(numel(offset_hz),numel(rf_scales));
parfor n_offset=1:numel(offset_hz)
    drift=H+2*pi*offset_hz(n_offset)*Lz;
    row=zeros(1,numel(rf_scales));
    for n_scale=1:numel(rf_scales)
        rf_x=rf_scales(n_scale)*xy_rad(1,:);
        rf_y=rf_scales(n_scale)*xy_rad(2,:);
        final_state=shaped_pulse_xy(spin_system,drift,{Lx,Ly},...
                    {rf_x,rf_y},pulse_dt,initial,'expv-pwc');
        row(n_scale)=real(trace(target'*final_state))/3;
    end
    fidelity(n_offset,:)=row;
end

% Convert the three-state SO(3) score to the source SU(2) overlap
quality=sqrt(max(0,(1+3*fidelity)/4));
profiles.offset_hz=offset_hz;
profiles.rf_scales=rf_scales;
profiles.fidelity=fidelity;
profiles.quality=quality;

% Display the source waveform and its simulated response
fig=kfigure();
tiledlayout(2,2);
time_us=1e6*cumsum(pulse_dt);
nexttile;
plot(time_us,xy_rad(1,:)/(2*pi*1e3),time_us,xy_rad(2,:)/(2*pi*1e3));
xlabel('time (microseconds)'); ylabel('RF quadratures (kHz)');
title('deposited UR180 pulse'); legend({'x','y'}); grid on;
nexttile;
plot(time_us,hypot(xy_rad(1,:),xy_rad(2,:))/(2*pi*1e3));
xlabel('time (microseconds)'); ylabel('RF amplitude (kHz)');
title('RF at nominal power'); grid on;
nexttile;
imagesc(offset_hz/1e3,rf_scales,quality.'); axis xy;
xlabel('offset (kHz)'); ylabel('relative B_1');
title('y-axis 180 degree overlap'); colorbar;
clim([min(quality(:)) 1]);
nexttile;
plot(offset_hz/1e3,quality);
xlabel('offset (kHz)'); ylabel('propagator overlap');
title('profiles by RF scale'); grid on;

end


