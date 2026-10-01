% Simulates the BEBOP and BIBOP pulses published by Kobzar et al.
% Syntax:
%
%    [profiles,fig]=bebop_bibop_profiles(offset_hz,rf_scales)
%
% Parameters:
%
%    offset_hz - resonance offsets in Hz; 200 points across
%                20 kHz form the verification grid
%
%    rf_scales - relative RF amplitudes; five evenly spaced
%                values from 0.8 to 1.2 reproduce the source grid
%
% Outputs:
%
%    profiles  - excitation and inversion efficiencies on the
%                offset by RF-scale grid
%
%    fig       - waveform and robustness maps
%
% Source: Kobzar et al., J. Magn. Reson. 170, 236-243 (2004).
% DOI: 10.1016/j.jmr.2004.06.017. The numerical shapes are from
% https://www.ioc.kit.edu/luy/186.php and /luy/227.php.
% Both source files give Cartesian RF frequencies in Hz and slice
% durations in seconds; pulse675 and pulse615 are selected.
%
function [profiles,fig]=bebop_bibop_profiles(offset_hz,rf_scales)

% Magnetic field and isotope
sys.magnet=14.1;
sys.isotopes={'1H'};

% Chemical shift, ppm
inter.zeeman.scalar={0};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

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

% Source waveform files
source_dir=fileparts(mfilename('fullpath'));
wave_files={'bebop_20khz_20pct.dat','bibop_20khz_20pct.dat'};

% Target states and pulse labels
targets={rho_x,-rho_z};
names={'BEBOP excitation','BIBOP inversion'};

% Response grid
profiles.offset_hz=offset_hz;
profiles.rf_scales=rf_scales;

% Start the waveform and response figure
fig=kfigure();
tiledlayout(2,2);

% Simulate the published shapes over the offset and RF-scale grid
for pulse_index=1:2

    % Read the author-provided waveform
    waveform=readmatrix(fullfile(source_dir,wave_files{pulse_index}),...
                        'NumHeaderLines',5);
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

            % Propagate the initial state
            final_state=shaped_pulse_xy(spin_system,drift,{lx,ly},...
                                        {rf_x,rf_y},pulse_dt,rho_z,'expv-pwc');

            % Score the excitation or inversion transfer
            row(n_scale)=real(targets{pulse_index}'*final_state); %#ok<PFBNS>

        end

        % Store the response row
        fidelity(n_offset,:)=row;

    end

    % Store the pulse response
    if pulse_index==1
        profiles.bebop=fidelity;
    else
        profiles.bibop=fidelity;

    end

    % Plot the deposited waveform in physical units
    nexttile(2*pulse_index-1);
    time_us=1e6*cumsum(pulse_dt);
    plot(time_us,xy_rad(1,:)/(2*pi*1e3),...
         time_us,xy_rad(2,:)/(2*pi*1e3));
    kxlabel('time (microseconds)'); kylabel('RF quadratures (kHz)');
    ktitle(names{pulse_index}); klegend({'x','y'}); kgrid;

    % Plot the independently propagated offset and B1 response
    nexttile(2*pulse_index);
    imagesc(offset_hz/1e3,rf_scales,fidelity.'); axis xy;
    kxlabel('offset (kHz)'); kylabel('relative $B_1$');
    ktitle('transfer efficiency'); colorbar; clim([min(fidelity(:)) 1]);

end

end


