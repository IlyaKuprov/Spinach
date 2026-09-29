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

% Build one spin in the rotating frame
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Construct normalised Cartesian states and control operators
Sz=state(spin_system,'Lz',1);
Sz=Sz/norm(full(Sz),2);
Sx=state(spin_system,'Lx',1);
Sx=Sx/norm(full(Sx),2);
Lx=operator(spin_system,'Lx',1);
Ly=operator(spin_system,'Ly',1);
Lz=operator(spin_system,'Lz',1);
H=hamiltonian(assume(spin_system,'nmr'));

% Load the author-provided excitation and inversion waveforms
source_dir=fileparts(mfilename('fullpath'));
wave_files={'bebop_20khz_20pct.dat','bibop_20khz_20pct.dat'};
targets={Sx,-Sz};
names={'BEBOP excitation','BIBOP inversion'};
profiles.offset_hz=offset_hz;
profiles.rf_scales=rf_scales;
fig=kfigure();
tiledlayout(2,2);

% Simulate the published shapes over the offset and RF-scale grid
for pulse_index=1:2
    waveform=readmatrix(fullfile(source_dir,wave_files{pulse_index}),...
                        'NumHeaderLines',5);
    assert(size(waveform,2)==3);
    pulse_dt=waveform(:,3).';
    xy_rad=2*pi*waveform(:,1:2).';
    fidelity=zeros(numel(offset_hz),numel(rf_scales));
    parfor n_offset=1:numel(offset_hz)
        drift=H+2*pi*offset_hz(n_offset)*Lz;
        row=zeros(1,numel(rf_scales));
        for n_scale=1:numel(rf_scales)
            rf_x=rf_scales(n_scale)*xy_rad(1,:);
            rf_y=rf_scales(n_scale)*xy_rad(2,:);
            final_state=shaped_pulse_xy(spin_system,drift,{Lx,Ly},...
                        {rf_x,rf_y},pulse_dt,Sz,'expv-pwc');
            row(n_scale)=real(targets{pulse_index}'*final_state);
        end
        fidelity(n_offset,:)=row;
    end
    if pulse_index==1
        profiles.bebop=fidelity;
    else
        profiles.bibop=fidelity;
    end

    % Plot the deposited waveform in physical units
    nexttile(2*pulse_index-1);
    time_us=1e6*cumsum(pulse_dt);
    plot(time_us,xy_rad(1,:)/(2*pi*1e3),time_us,xy_rad(2,:)/(2*pi*1e3));
    xlabel('time (microseconds)'); ylabel('RF quadratures (kHz)');
    title(names{pulse_index}); legend({'x','y'}); grid on;

    % Plot the independently propagated offset and B1 response
    nexttile(2*pulse_index);
    imagesc(offset_hz/1e3,rf_scales,fidelity.'); axis xy;
    xlabel('offset (kHz)'); ylabel('relative B_1');
    title('transfer efficiency'); colorbar; clim([min(fidelity(:)) 1]);
end

end


