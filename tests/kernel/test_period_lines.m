% Test: kernel/period_lines | Complex periodic spectral lines
% Checks signed Floquet lines against native complex time acquisition.
% Syntax: result=test_period_lines()
%
% Outputs:
%
%   result - line reconstruction, coincident lines, and finite-window checks
%
% talos@spindynamics.org

function result=test_period_lines()

% Build one driven period with distinct, noncommuting Hamiltonians
result=new_test_result('kernel/period_lines','Complex periodic spectral lines',...
                       'Summed lines must reconstruct directly sampled complex FIDs.');
spin_system.bas.formalism='zeeman-hilb';
spin_system.tols.liouv_zero=1e-10;
parameters.rate=2500;
parameters.npoints=45;
parameters.start_tick=2;
parameters.tick_rf={};
parameters.dwell=1/(3*parameters.rate);
parameters.rho0=[0.2 0.3i;-0.3i 0.8];
parameters.coil=[0 1;0 0];
parameters.sum_up=false;
X=[0 1;1 0]/2;
Y=[0 -1i;1i 0]/2;
Z=diag([1 -1])/2;
H={2*pi*720*Z,2*pi*1030*X,2*pi*330*Y};
zero=zeros(2);

% Reconstruct every dwell from the unthresholded signed amplitudes
time=(0:parameters.npoints-1)*parameters.dwell;
fid=periodic_acquire(spin_system,parameters,H,zero,zero);
lines=periodic_lines(spin_system,parameters,H,zero,zero);
signal=lines(:,2).'*exp(-2i*pi*real(lines(:,1))*time);
result=test_close(result,'noncommuting complex reconstruction',signal,fid,1e-11,0,...
                  'Folded lines retain the phase and amplitude of the complex FID.');

% Check both a bin-centred tone and an off-bin tone in the same finite window
for freq=[500 537]

    % Use a signed one-spin frequency with a complex receiver
    H={2*pi*freq*Z};
    parameters.start_tick=1;
    parameters.dwell=1/parameters.rate;
    parameters.npoints=25;
    fid=periodic_acquire(spin_system,parameters,H,zero,zero);
    lines=periodic_lines(spin_system,parameters,H,zero,zero);
    time=(0:parameters.npoints-1)*parameters.dwell;
    bell=exp(-time/time(end));
    signal=lines(:,2).'*exp(-2i*pi*real(lines(:,1))*time);
    result=test_close(result,sprintf('finite spectrum at %g Hz',freq),...
                      fft(signal.*bell),fft(fid.*bell),1e-11,0,...
                      'Line synthesis and sampling use the same finite observation window.');

end

% Degenerate identity periods must sum all coincident amplitudes before output
Q=[1 1i;1i 1]/sqrt(2);
rho=parameters.rho0;
coil=parameters.coil;
[freq,amp,band]=period_lines({eye(2),eye(2)},parameters.rate,rho,coil);
result=test_close(result,'degenerate line',[freq amp band],[0 trace(coil'*rho) 0],1e-13,0,...
                  'All zero-frequency transitions contribute to the one physical line.');
[~,rotated,~]=period_lines({eye(2),eye(2)},parameters.rate,Q'*rho*Q,Q'*coil*Q);
result=test_close(result,'degenerate basis invariance',rotated,amp,1e-13,0,...
                  'A unitary change of basis in the degenerate subspace changes no amplitude.');

end


