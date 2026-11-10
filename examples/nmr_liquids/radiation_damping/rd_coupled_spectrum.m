% Radiation-damped AB proton spectrum and Cartesian magnetisation
% traces, compared with the same spin system without feedback.
% Exponential apodisation is applied only to the plotted spectra.
%
% Ported from the Spynach browser example try.html#ex-rd, copyright
% (c) 2026 Alexej Jerschow, MIT licence:
% https://github.com/ajerschow/spynach/blob/main/LICENSE
%
% Calculation time: seconds
%
% Alexej Jerschow (original demonstration)
% talos@spindynamics.org

function rd_coupled_spectrum()

% Magnet field
sys.magnet=14.1;

% Isotopes
sys.isotopes={'1H','1H'};

% Chemical shifts
inter.zeeman.scalar={0 0};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation={'none'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% AB Hamiltonian in the centre of the proton multiplet
J=10; dnu=12; rate=60; phase=0; theta=170*pi/180;
H=pi*dnu*(operator(spin_system,'Lz',1)-...
          operator(spin_system,'Lz',2))+...
  2*pi*J*(operator(spin_system,{'Lx','Lx'},{1,2})+...
          operator(spin_system,{'Ly','Ly'},{1,2})+...
          operator(spin_system,{'Lz','Lz'},{1,2}));

% Rotation generators and detection states
Lx=operator(spin_system,'Lx','1H');
Ly=operator(spin_system,'Ly','1H');
coil=[state(spin_system,'Lx','1H') ...
      state(spin_system,'Ly','1H') ...
      state(spin_system,'Lz','1H')];
rho_eq=coil(:,3); mz_eq=real(rho_eq'*rho_eq);
gain=rate*cos(phase)*exp(1i*phase);
field=@(rho)-1i*gain*(coil(:,1)'+1i*coil(:,2)')*rho/mz_eq;
G=@(t,rho)H+real(field(rho))*Lx+imag(field(rho))*Ly;

% Radiation-damped and uncoupled-from-probe trajectories
dt=5e-4; nsteps=800; time=(0:nsteps)*dt;
rho=step(spin_system,Ly,rho_eq,theta); rho_free=rho;
magn=zeros(3,nsteps+1,2);
magn(:,1,1)=real(coil'*rho)/mz_eq; magn(:,1,2)=magn(:,1,1);
for n=1:nsteps
    rho=iserstep(spin_system,{G,time(n),'LG4'},rho,dt);
    rho_free=step(spin_system,H,rho_free,dt);
    magn(:,n+1,1)=real(coil'*rho)/mz_eq;
    magn(:,n+1,2)=real(coil'*rho_free)/mz_eq;
end

% One-hertz exponential line broadening and Fourier transform
fid=squeeze(magn(1,1:end-1,:)+1i*magn(2,1:end-1,:));
fid=apodisation(spin_system,fid,{{'exp',pi*time(end-1)},{}});
nfft=8*nsteps; spectrum=fftshift(fft(fid,nfft),1);
frequency=(-nfft/2:nfft/2-1)/(nfft*dt);

% Spectra and Cartesian magnetisation traces
kfigure();
subplot(1,2,1); plot(frequency,abs(spectrum));
xlim([-60 60]); set(gca,'XDir','reverse');
xlabel('frequency, Hz'); ylabel('|spectrum|, a.u.');
ylim([0 1.1*max(abs(spectrum),[],'all')]);
legend('radiation damping','no feedback','Location','best'); kgrid;
subplot(1,2,2); plot(time,magn(:,:,1),time,magn(3,:,2),'--');
xlabel('time, seconds'); ylabel('magnetisation / equilibrium');
legend('m_x','m_y','m_z','m_z, no feedback','Location','best'); kgrid;

end

