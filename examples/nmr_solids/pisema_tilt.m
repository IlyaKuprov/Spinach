% Oriented PISEMA/FSLG for a static nitrogen-proton spin pair.
% The source two-channel preparation and signed offsets are retained.
% The spectra compare a fixed tilted dipole with its aligned counterpart.
%
% talos@spindynamics.org

function pisema_tilt()

% Nitrogen-proton pair with a 10 kHz axial dipole
sys.magnet=2*pi*700e6/spin('1H');
sys.isotopes={'15N','1H'};
sys.parallel={'processes',1};
inter.zeeman.scalar={0 0};
inter.coupling.matrix={[] 1e4*diag([-1 -1 2]); [] []};
bas.formalism='zeeman-hilb';
bas.approximation={'none'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Preparation, detection, and two-pulse FSLG dwell
parameters.spins=sys.isotopes;
parameters.offset=[0 0];
parameters.rf=50e3;
parameters.pulse_time=sqrt(2/3)/parameters.rf;
parameters.npoints=2048;
parameters.rho0=state(spin_system,'Lx','1H')-state(spin_system,'Ly','15N');
parameters.coil=coil_state(spin_system,'L+','1H');
parameters.nx=operator(spin_system,'Lx','15N');
parameters.ny=operator(spin_system,'Ly','15N');
parameters.nz=operator(spin_system,'Lz','15N');
parameters.hx=operator(spin_system,'Lx','1H');

% Change only the static crystal alignment
tilts=[pi/6 0];
fids=complex(zeros(2,parameters.npoints));
for n=1:2
    parameters.orientation=[0 tilts(n) 0];
    fids(n,:)=crystal(spin_system,@pisema_trace,parameters,'nmr');
end

% Check the tilted trace independently with explicit spin matrices
Ix=kron([0 1; 1 0]/2,eye(2));
Iy=kron([0 -1i; 1i 0]/2,eye(2));
Iz=kron(diag([1 -1])/2,eye(2));
Sx=kron(eye(2),[0 1; 1 0]/2);
Sp=kron(eye(2),[0 1; 0 0]);
Sz=kron(eye(2),diag([1 -1])/2);
H=2*pi*1e4*(3*cos(tilts(1))^2-1)*Iz*Sz;
rf=2*pi*parameters.rf;
theta=pi/2-atan(sqrt(2));
P=expm(1i*theta*Ix);
rho=P*(Sx-Iy)*P';
drive=rf*(Iy+sqrt(1.5)*Sx-Iz/sqrt(2));
P=expm(-1i*(H-drive)*parameters.pulse_time)*...
  expm(-1i*(H+drive)*parameters.pulse_time);
ref=complex(zeros(1,16));
for n=1:numel(ref)
    ref(n)=trace(Sp'*rho);
    rho=P*rho*P';
end
assert(norm(fids(1,1:numel(ref))-ref)<1e-10);
assert(norm(fids(1,:)-fids(2,:))>sqrt(eps));

% Apodise and display the proton-detected indirect spectra
dwell=2*parameters.pulse_time;
time=(0:parameters.npoints-1)*dwell;
nfft=4096;
spectra=fftshift(fft(fids.*exp(-pi*20*time),nfft,2),2);
freq=(-nfft/2:nfft/2-1)/(nfft*dwell);
kfigure(); plot(freq/1e3,real(spectra).');
kxlabel('Indirect frequency, kHz');
kylabel('Real spectral amplitude'); kgrid;
klegend({'Tilted by 30 degrees','Aligned'});

end

% Sample the proton receiver between finite FSLG cycles
function fid=pisema_trace(~,parameters,H,~,~)

% Prepare with an instantaneous nitrogen rotation about minus X
theta=pi/2-atan(sqrt(2));
U=expm(1i*theta*full(parameters.nx));
rho=U*parameters.rho0*U';

% Flip both RF phases and the nitrogen offset on the second segment
drive=2*pi*parameters.rf*(parameters.ny+sqrt(1.5)*parameters.hx-...
                          parameters.nz/sqrt(2));
U=expm(-1i*full(H-drive)*parameters.pulse_time)*...
  expm(-1i*full(H+drive)*parameters.pulse_time);

% Acquire with no intervening free evolution
fid=complex(zeros(1,parameters.npoints));
for n=1:numel(fid)
    fid(n)=trace(parameters.coil'*rho);
    rho=U*rho*U';
end

end


