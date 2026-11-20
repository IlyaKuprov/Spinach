% Periodic distant dipolar field, CRAZED, and diffusion tests. Syntax:
%
%                    result=test_ddf_imaging()
%
% Outputs:
%
%    result - signed analytical and coupled-pair Fourier reference checks

function result=test_ddf_imaging()

% Isolate the uniform and grating symbols for every field component
result=new_test_result('kernel/ddf_imaging','Periodic distant dipolar fields',...
                       'Analytical trajectories and a coupled Fourier-Galerkin reference.');
mu0=1.25663706127e-6;
wave=[1;-1]*[1 2 3];
magn=wave+[3 1 2];
field=ddf_field(magn,[0 0 1],0);
target=(wave-[3 1 2]/2).*[1 1 -2]/3;
result=test_close(result,'shape and grating',field/mu0,target,1e-13,0,...
                  'only the mean sees the axial needle boundary');
field=ddf_field(magn,[sqrt(2/3) 0 sqrt(1/3)],1/3);
result=test_close(result,'field null',field/mu0,zeros(2,3),1e-14,0,...
                  'magic-angle grating and spherical mean are both silent');

% Build a single-proton system and calibrate the raw transverse receiver
sys.magnet=1;
sys.isotopes={'1H'};
sys.output='hush';
sys.parallel={'processes',1};
sys.disable={'hygiene'};
bas.formalism='sphten-liouv';
bas.approximation={'none'};
s=basis(create(sys,struct()),bas);
parameters=phantom(s,24);
parameters.flip_angles=[pi/2 0];
pulse=imaging(s,@crazed_ddf,parameters);
gain=pulse(1)/1i;
result=test_close(result,'raw pulse phase',gain,parameters.coil_st{1}'*state(s,'Lx','1H'),...
                  1e-10,0,'a positive X pulse makes negative Y magnetisation');
time=(0:parameters.npoints-1)'/parameters.sweep;

% Compare the entire complex trajectory including precession direction
parameters.flip_angles=[pi/3 0];
for shape=[1/3 0 1]
    parameters.demag_factor=shape;
    signal=imaging(s,@crazed_ddf,parameters)/gain;
    rate=10*cos(pi/3)*(3*shape-1)/2;
    target=1i*sin(pi/3)*exp(-1i*rate*time);
    result=test_close(result,sprintf('tilt Nz=%g',shape),...
                      max(abs(signal-target)),0,1e-6,0,'signed uniform precession');
end

% Detect longitudinal magnetisation without changing the source calibration
parameters.coil_st={parameters.rho_ref};
signal=imaging(s,@crazed_ddf,parameters);
target=cos(pi/3)*(parameters.rho_ref'*parameters.rho_ref);
result=test_close(result,'conserved longitudinal component',...
                  max(abs(signal-target)),0,1e-7,0,'the tilted uniform state retains Mz');

% Isolate the imposed Hamiltonian offset with no dipolar feedback
parameters=phantom(s,24);
parameters.flip_angles=[pi/2 0];
parameters.offset=3;
parameters.magnetisation=0;
signal=imaging(s,@crazed_ddf,parameters)/gain;
result=test_close(result,'offset sign',max(abs(signal-1i*exp(-6i*pi*time))),...
                  0,1e-7,0,'positive rotating-frame offset rotates L+ detection negatively');

% Select the first three grating orders after two ninety-degree pulses
parameters=phantom(s,24);
area=2*pi/(spin('1H')*parameters.dims);
parameters.flip_angles=[pi/2 pi/2];
for zcos=[1 0 sqrt(1/3)]
    parameters.direction=[sqrt(1-zcos^2) 0 zcos];
    parameters.demag_factor=zcos^2;
    arg=10*(3*zcos^2-1)*time/2;
    for n=1:3
        parameters.gradient_areas=area*[1 n];
        target=-(1i)^n*(besselj(n-1,arg)+besselj(n+1,arg))/2;
        signal=imaging(s,@crazed_ddf,parameters)/gain;
        result=test_close(result,sprintf('Bessel order=%d cos=%g',n,zcos),...
                          max(abs(signal-target)),0,1e-5,0,'complex echo in single-pulse units');
    end
end

% Detect the encoded harmonic directly with a fixed phase-sensitive coil
parameters.magnetisation=0;
parameters.flip_angles=[pi/2 0];
parameters.gradient_areas=area*[2 0];
parameters.diff=2e-9;
phase=2*pi*((0:parameters.npts-1)'+0.5-parameters.npts/2)/parameters.npts;
parameters.coil_ph={exp(2i*phase)/parameters.npts};
signal=imaging(s,@crazed_ddf,parameters)/gain;
target=1i*exp(-parameters.diff*(4*pi/parameters.dims)^2*time);
result=test_close(result,'diffusion of second harmonic',...
                  max(abs(signal-target)),0,1e-7,0,'Fourier diffusion damps the encoded mode');

% Resolve unequal offsets and scalar coupling in the presence of diffusion
sys.isotopes={'1H','1H'};
shifts=[9 -13];
coupling=11;
inter.zeeman.scalar=num2cell(-2*pi*1e6*shifts/spin('1H'));
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=coupling;
s=basis(create(sys,inter),bas);
parameters=phantom(s,32);
parameters.diff=2e-9;
parameters.flip_angles=[1.1 0.9];
parameters.gradient_areas=area*[1 -2];
parameters.demag_factor=1/3;
parameters.sweep=400;
parameters.ddf_steps=8;
ref=pair_reference(parameters,shifts,coupling,12);
wide=pair_reference(parameters,shifts,coupling,18);
result=test_close(result,'Galerkin bandwidth',max(abs(ref-wide)),0,1e-9,0,...
                  'zero-padded convolution converges under independent mode refinement');

% Verify that the receiver distinguishes scalar coupling and transport
uncoupled=pair_reference(parameters,shifts,0,12);
static=parameters;
static.diff=0;
immobile=pair_reference(static,shifts,coupling,12);
result=test_true(result,'resolved coupling',max(abs(uncoupled-ref))>1e-5,...
                 'omitting scalar coupling changes the trace above its acceptance tolerance');
result=test_true(result,'resolved transport',max(abs(immobile-ref))>1e-5,...
                 'omitting diffusion changes the trace above its acceptance tolerance');

% Refine the time step and then the voxel grid at fixed physical length
coarse=imaging(s,@crazed_ddf,parameters);
parameters.ddf_steps=2*parameters.ddf_steps;
fine=imaging(s,@crazed_ddf,parameters);
parameters.npts=2*parameters.npts;
parameters.rho0_ph={ones(parameters.npts,1)};
parameters.coil_ph={ones(parameters.npts,1)/parameters.npts};
parameters.rlx_ph={ones(parameters.npts,1)};
finer=imaging(s,@crazed_ddf,parameters);
parameters.flip_angles=[pi/2 0];
parameters.gradient_areas=[0 0];
parameters.npoints=1;
pulse=imaging(s,@crazed_ddf,parameters);
gain=pulse(1)/1i;
result=test_close(result,'time refinement',max(abs((fine-coarse)/gain)),0,1e-6,0,...
                  'halving the integration step resolves the coupled signal');
result=test_close(result,'voxel refinement',max(abs((finer-fine)/gain)),0,1e-6,0,...
                  'doubling the number of voxels resolves the coupled signal');
result=test_close(result,'two-spin Fourier reference',max(abs(finer/gain-wide)),0,1e-6,0,...
                  'Pauli commutators and non-circular mode convolution reproduce complex reception');

end

% Specify a homogeneous periodic phantom with a volume-average receiver
function parameters=phantom(s,count)

% Sample a 0.5 mm pitch with a 100 ms dipolar time
parameters.npts=count;
parameters.dims=5e-4;
parameters.deriv={'fourier'};
parameters.diff=0;
parameters.spins={'1H'};
parameters.offset=0;
parameters.direction=[0 0 1];
parameters.demag_factor=1/3;
parameters.magnetisation=10/(spin('1H')*1.25663706127e-6);
parameters.rho_ref=state(s,'Lz','1H');
parameters.rho0_st={parameters.rho_ref};
parameters.rho0_ph={ones(count,1)};
parameters.coil_st={coil_state(s,'L+','1H')};
parameters.coil_ph={ones(count,1)/count};
parameters.rlx_ph={ones(count,1)};
parameters.rlx_op={sparse(numel(parameters.rho_ref),numel(parameters.rho_ref))};
parameters.gradient_areas=[0 0];
parameters.sweep=80;
parameters.npoints=9;
parameters.ddf_steps=16;

end

% Propagate a two-proton Hilbert density matrix for each spatial harmonic
function fid=pair_reference(parameters,shifts,coupling,halfwidth)

% Assemble scalar coupling and magnetic moments from Cartesian Pauli matrices
pauli={[0 1;1 0]/2,[0 -1i;1i 0]/2,diag([1 -1])/2};
mom=zeros(16,3);
rot=zeros(16,16,3);
H=zeros(4);
for k=1:3
    first=kron(pauli{k},eye(2));
    second=kron(eye(2),pauli{k});
    total=first+second;
    mom(:,k)=total(:);
    rot(:,:,k)=-1i*(kron(eye(4),total)-kron(total.',eye(4)));
    H=H+2*pi*coupling*first*second;
end
H=H+2*pi*(shifts(1)*first+shifts(2)*second);
linear=-1i*(kron(eye(4),H)-kron(H.',eye(4)));

% Rotate locally in spin space and shift each matrix element by coherence order
modes=-halfwidth:halfwidth;
rho=zeros(16,numel(modes));
rho(:,halfwidth+1)=mom(:,3);
orders=diag(total)-diag(total).';
for k=1:2
    rho=expm(parameters.flip_angles(k)*rot(:,:,1))*rho;
    pitch=round(parameters.gradient_areas(k)*spin('1H')*parameters.dims/(2*pi));
    for n=1:16
        shift=orders(n)*pitch;
        keep=abs(modes-shift)<=halfwidth;
        rho(n,keep)=rho(n,find(keep)-shift);
        rho(n,~keep)=0;
    end
end

% Differentiate the Fourier series by a non-circular coefficient convolution
weight=ones(size(modes))*(3*parameters.direction(3)^2-1)/2;
weight(halfwidth+1)=(3*parameters.demag_factor-1)/2;
field_scale=-spin('1H')*1.25663706127e-6*parameters.magnetisation*[1 1 -2]/(3*norm(mom(:,3))^2);
decay=parameters.diff*(2*pi*modes/parameters.dims).^2;
time=(0:parameters.npoints-1)'/parameters.sweep;
opts=odeset('AbsTol',1e-13,'RelTol',1e-11);
[~,traj]=ode45(@rhs,time,rho(:),opts);
coil=mom(:,1)+1i*mom(:,2);
fid=traj(:,16*halfwidth+(1:16))*conj(coil)/norm(mom(:,3))^2;

    % Couple each moment harmonic to the retained density harmonics
    function slope=rhs(~,vec)
        density=reshape(vec,16,[]);
        change=linear*density-density.*decay;
        for a=1:3
            coeff=field_scale(a)*weight.*(mom(:,a)'*density);
            change=change+conv2(rot(:,:,a)*density,coeff,'same');
        end
        slope=change(:);
    end

end


