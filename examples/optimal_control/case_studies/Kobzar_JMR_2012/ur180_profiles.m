% Simulates an author-deposited robust universal 180 degree pulse
% Syntax:
%
%    [profiles,fig]=ur180_profiles(offset_hz,rf_scales)
%
% Parameters:
%
%    offset_hz - resonance offsets in Hz; 100 points over
%                20 kHz reproduce the source grid
%
%    rf_scales - relative RF amplitudes; five evenly spaced
%                values from 0.6 to 1.4 reproduce the source grid
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

% Normalised Cartesian states
rho_x=state(spin_system,'Lx',1);
rho_x=rho_x/norm(full(rho_x),2);
rho_y=state(spin_system,'Ly',1);
rho_y=rho_y/norm(full(rho_y),2);
rho_z=state(spin_system,'Lz',1);
rho_z=rho_z/norm(full(rho_z),2);

% Control and offset operators
lx=operator(spin_system,'Lx',1);
ly=operator(spin_system,'Ly',1);
lz=operator(spin_system,'Lz',1);

% Drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Initial and target mappings for a y-axis pi rotation
initial=[rho_x rho_y rho_z];
target=[-rho_x rho_y -rho_z];

% Read the deposited x/y RF frequencies and slice durations
source_dir=fileparts(mfilename('fullpath'));
waveform=readmatrix(fullfile(source_dir,'ur180_20khz_40pct.dat'));
assert(size(waveform,2)==3);

% Slice durations and Cartesian controls in rad/s
pulse_dt=waveform(:,3).';
xy_rad=2*pi*waveform(:,1:2).';

% Preallocate the offset-by-RF response
fidelity=zeros(numel(offset_hz),numel(rf_scales));

% Loop over resonance offsets
parfor n_offset=1:numel(offset_hz)

    % Offset Hamiltonian and response row
    drift=H+2*pi*offset_hz(n_offset)*lz;
    row=zeros(1,numel(rf_scales));

    % Loop over RF scalings
    for n_scale=1:numel(rf_scales)

        % Scale the deposited waveform
        rf_x=rf_scales(n_scale)*xy_rad(1,:); %#ok<PFBNS>
        rf_y=rf_scales(n_scale)*xy_rad(2,:);

        % Propagate all three Cartesian states
        final_state=shaped_pulse_xy(spin_system,drift,{lx,ly},...
                                    {rf_x,rf_y},pulse_dt,initial,'expv-pwc');

        % Score the three-state rotation
        row(n_scale)=real(trace(target'*final_state))/3;

    end

    % Store the response row
    fidelity(n_offset,:)=row;

end

% Convert the three-state SO(3) score to the source SU(2) overlap
quality=sqrt(max(0,(1+3*fidelity)/4));

% Package the grid and rotation scores
profiles.offset_hz=offset_hz;
profiles.rf_scales=rf_scales;
profiles.fidelity=fidelity;
profiles.quality=quality;

% Start the waveform and response figure
fig=kfigure();
tiledlayout(2,2);

% Plot the deposited Cartesian waveform
time_us=1e6*cumsum(pulse_dt);
nexttile;
plot(time_us,xy_rad(1,:)/(2*pi*1e3),...
     time_us,xy_rad(2,:)/(2*pi*1e3));
kxlabel('time (microseconds)'); kylabel('RF quadratures (kHz)');
ktitle('deposited UR180 pulse'); klegend({'x','y'}); kgrid;

% Plot the nominal RF amplitude
nexttile;
plot(time_us,hypot(xy_rad(1,:),xy_rad(2,:))/(2*pi*1e3));
kxlabel('time (microseconds)'); kylabel('RF amplitude (kHz)');
ktitle('RF at nominal power'); kgrid;

% Plot the offset-by-RF robustness map
nexttile;
imagesc(offset_hz/1e3,rf_scales,quality.'); axis xy;
kxlabel('offset (kHz)'); kylabel('relative $B_1$');
ktitle('y-axis 180 degree overlap'); colorbar;
clim([min(quality(:)) 1]);

% Plot the RF-scale response profiles
nexttile;
plot(offset_hz/1e3,quality);
kxlabel('offset (kHz)'); kylabel('propagator overlap');
ktitle('profiles by RF scale'); kgrid;

end


