% Tests deterministic pulse waveform generators. Syntax:
%
%                    result=test_pulses_waveform_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks analytic waveform formulas, periodic phase tables,
% JCAMP pulse-file reading, chirp construction, and sech pulse coordinates.
%
% ilya.kuprov@weizmann.ac.il

function result=test_pulses_waveform_suite()

% State the waveform target of the test
result=new_test_result('kernel/pulses_waveform_suite',...
                       'Pulse waveform generators',...
                       'pulse waveform utilities must reproduce their analytic formulae and tabulated periodic sequences.');

% Check sawtooth and triangular wave formulae at simple fractions of the period
amp=2;
freq=1;
time_grid=[0 0.25 0.5 0.75 1.0];
saw_ref=[-2 -1 0 1 -2];
tri_ref=abs(saw_ref);
result=test_close(result,'sawtooth explicit samples',sawtooth(amp,freq,time_grid),saw_ref,1e-15,1e-15,...
                  'sawtooth() implements amplitude*(2*f*mod(t,1/f)-1)');
result=test_close(result,'triwave explicit samples',triwave(amp,freq,time_grid),tri_ref,1e-15,1e-15,...
                  'triwave() is the absolute value of the corresponding sawtooth');

% Check Uhrig delay formula for three pulses over a unit interval
T=1;
N=3;
pos=T*(sin(pi*(1:N)/(2*N+2)).^2-0.5);
delays=diff(pos);
chunk=(T-sum(delays))/2;
udd_ref=[chunk delays chunk];
result=test_close(result,'uhrig_times N=3',uhrig_times(T,N),udd_ref,1e-15,1e-15,...
                  'UDD pulse positions follow sin^2(pi*k/(2N+2)) and delays are adjacent differences');

% Check periodic phase tables and wrap-around indexing
result=test_close(result,'pmlg5 first phase',pmlg5(1),(pi/180)*339.22,1e-14,1e-14,...
                  'PMLG5 index one is the first tabulated phase');
result=test_close(result,'pmlg5 wrap phase',pmlg5(21),pmlg5(1),1e-14,1e-14,...
                  'PMLG5 has a 20-pulse periodic phase table');
result=test_close(result,'spinal first phase',spinal(1),(pi/180)*10,1e-14,1e-14,...
                  'SPINAL index one is the first tabulated phase');
result=test_close(result,'spinal wrap phase',spinal(65),spinal(1),1e-14,1e-14,...
                  'SPINAL has a 64-pulse periodic phase table');

% Check simple analytic pulse envelopes
result=test_close(result,'pulse_shape rectangular',pulse_shape('rectangular',4),ones(1,4),1e-15,1e-15,...
                  'the rectangular envelope is constant over all pulse slices');
result=test_close(result,'pulse_shape sinc3 nodes',pulse_shape('sinc3',3),[0 pi 0],1e-15,1e-15,...
                  'sinc3 sampled at -3, 0, and +3 gives 0, pi, and 0');

% Check reading of a distributed rectangular Bruker pulse file
[A,phi,Cx,Cy,scaling_factor]=read_wave('rectangular_1000.pk',4);
result=test_close(result,'read_wave rectangular amplitude',A,ones(1,4),1e-15,1e-15,...
                  'the distributed rectangular pulse has 100 percent amplitude at every point');
result=test_close(result,'read_wave rectangular phase',phi,zeros(1,4),1e-15,1e-15,...
                  'the distributed rectangular pulse has zero phase at every point');
result=test_close(result,'read_wave rectangular x',Cx,ones(1,4),1e-15,1e-15,...
                  'zero-phase polar coordinates convert into unit X control');
result=test_close(result,'read_wave rectangular y',Cy,zeros(1,4),1e-15,1e-15,...
                  'zero-phase polar coordinates have no Y control');
result=test_close(result,'read_wave scaling factor',scaling_factor,1,1e-15,1e-15,...
                  'the rectangular pulse file declares unit integral scaling');

% Check Veshtort-Griffin duration scaling without duplicating the coefficient table
vg_short=vg_pulse('E0A',7,1);
vg_long=vg_pulse('E0A',7,2);
result=test_close(result,'vg_pulse duration scaling',vg_long,vg_short/2,1e-14,1e-14,...
                  'VG pulse amplitudes are normalised as 2*pi*shape/duration');

% Check WURST chirp formula on a grid whose samples are simple fractions
[Cx,Cy,durs,ints,amps,phis,frqs]=chirp_pulse(5,1,4,2,'wurst');
time_grid=linspace(-0.5,0.5,5);
amps_ref=2*pi*sqrt(4)*(1-abs(sin(pi*time_grid).^2));
phis_ref=pi*4*(time_grid.^2);
frqs_ref=4*time_grid;
[x_ref,y_ref]=polar2cartesian(amps_ref,phis_ref);
result=test_close(result,'chirp_pulse wurst amplitudes',amps,amps_ref,1e-12,1e-12,...
                  'WURST amplitude is calibrated sqrt(bandwidth/duration)*(1-|sin(pi t)^p|)');
result=test_close(result,'chirp_pulse phases',phis,phis_ref,1e-12,1e-12,...
                  'linear chirp phase is pi*duration*bandwidth*t^2 on the normalised grid');
result=test_close(result,'chirp_pulse frequencies',frqs,frqs_ref,1e-12,1e-12,...
                  'instantaneous chirp frequency is bandwidth times normalised time');
result=test_close(result,'chirp_pulse durations',durs,ones(1,5)/5,1e-15,1e-15,...
                  'uniform piecewise-constant chirp slices divide the duration equally');
result=test_close(result,'chirp_pulse intervals',ints,diff(time_grid),1e-15,1e-15,...
                  'piecewise-linear intervals are differences of the waveform grid');
result=test_close(result,'chirp_pulse x control',Cx,x_ref,1e-12,1e-12,...
                  'chirp Cartesian X control is amplitude times cosine phase');
result=test_close(result,'chirp_pulse y control',Cy,y_ref,1e-12,1e-12,...
                  'chirp Cartesian Y control is amplitude times sine phase');

% Check hyperbolic secant pulse formulae
peak_amp=3;
freq_mod=2;
phase_mod=5;
dur=2;
npts=5;
[Cx,Cy,time_grid,amps,phis]=sech_pulse(peak_amp,freq_mod,phase_mod,dur,npts);
time_ref=linspace(-dur/2,dur/2,npts);
amps_ref=peak_amp*sech(freq_mod*time_ref);
phis_ref=phase_mod*log(cosh(freq_mod*time_ref));
[x_ref,y_ref]=polar2cartesian(amps_ref,phis_ref);
result=test_close(result,'sech_pulse time grid',time_grid,time_ref,1e-15,1e-15,...
                  'the sech pulse time grid is centred at zero');
result=test_close(result,'sech_pulse amplitudes',amps,amps_ref,1e-15,1e-15,...
                  'the sech pulse amplitude is peak_amp*sech(freq_mod*t)');
result=test_close(result,'sech_pulse phases',phis,phis_ref,1e-15,1e-15,...
                  'the sech pulse phase is phase_mod*log(cosh(freq_mod*t))');
result=test_close(result,'sech_pulse x control',Cx,x_ref,1e-15,1e-15,...
                  'sech pulse X control is the polar-to-Cartesian X component');
result=test_close(result,'sech_pulse y control',Cy,y_ref,1e-15,1e-15,...
                  'sech pulse Y control is the polar-to-Cartesian Y component');

end


