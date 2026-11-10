% Detuned-probe radiation damping: line broadening and frequency
% pulling of a weakly tilted proton with transverse relaxation.
% Small-tip predictions are compared with the simulated FIDs.
%
% Ported from Spynach, copyright (c) 2026 Alexej Jerschow, MIT licence:
% https://github.com/ajerschow/spynach/blob/main/LICENSE
%
% Calculation time: seconds
%
% Alexej Jerschow (original demonstrations)
% talos@spindynamics.org

function rd_detuned_probe()

% Magnet field
sys.magnet=14.1;

% Isotopes
sys.isotopes={'1H'};

% Chemical shifts
inter.zeeman.scalar={0};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation={'none'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Rotation generators and detection states
Lx=operator(spin_system,'Lx','1H');
Ly=operator(spin_system,'Ly','1H');
Lz=operator(spin_system,'Lz','1H');
coil=[state(spin_system,'Lx','1H') ...
      state(spin_system,'Ly','1H') ...
      state(spin_system,'Lz','1H')];
rho_eq=coil(:,3); mz_eq=real(rho_eq'*rho_eq);

% Probe phases, initial tilt, and transverse relaxation
rate=40; phases=[0 0.6 1]; theta=0.05; r2=2; offset=40;
dt=5e-4; nsteps=1000; time=(0:nsteps)*dt;
H=2*pi*offset*Lz-1i*r2*(Lz*Lz);
fid=zeros(nsteps+1,numel(phases));

% Shared-coil feedback for each probe tuning
for k=1:numel(phases)
    gain=rate*cos(phases(k))*exp(1i*phases(k));
    field=@(rho)-1i*gain*(coil(:,1)'+1i*coil(:,2)')*rho/mz_eq;
    G=@(t,rho)H+real(field(rho))*Lx+imag(field(rho))*Ly;
    rho=step(spin_system,Ly,rho_eq,theta);
    fid(1,k)=(coil(:,1)'+1i*coil(:,2)')*rho/mz_eq;
    for n=1:nsteps
        rho=iserstep(spin_system,{G,time(n),'LG4'},rho,dt);
        fid(n+1,k)=(coil(:,1)'+1i*coil(:,2)')*rho/mz_eq;
    end

    % Small-tip damping and frequency-pulling predictions
    decay=r2+real(gain); frequency=offset-imag(gain)/(2*pi);
    envelope=polyfit(time,log(abs(fid(:,k)))',1);
    phase_fit=polyfit(time,unwrap(angle(fid(:,k)))',1);
    fprintf('Phase %.1f: decay %.6f (%.6f), frequency %.6f (%.6f)\n',...
            phases(k),-envelope(1),decay,phase_fit(1)/(2*pi),frequency);
end

% FID envelopes and broadened, shifted spectra
nfft=8192; frequency=(-nfft/2:nfft/2-1)/(nfft*dt);
spectrum=fftshift(fft(fid(1:end-1,:),nfft),1);
kfigure();
subplot(1,2,1); plot(time,abs(fid));
xlabel('time, seconds'); ylabel('|m_x+im_y| / equilibrium');
legend('phase 0','phase 0.6','phase 1','Location','best'); kgrid;
subplot(1,2,2); plot(frequency,abs(spectrum));
xlim([20 60]); set(gca,'XDir','reverse');
xlabel('frequency, Hz'); ylabel('|spectrum|, a.u.');
legend('phase 0','phase 0.6','phase 1','Location','best'); kgrid;

end

