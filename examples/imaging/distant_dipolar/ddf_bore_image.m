% Synthetic proton CRAZED projection of an 8-by-10 mm water cylinder
% with an off-centre 2 mm axial bore, a uniform coil, and no magnetic
% susceptibility contrast. Finite rectangular lobes wind one and two
% turns per mm, with hard X pulses before each lobe. The single-species
% Bloch reduction evolves diffusion, relaxation, and dipolar torque
% during both lobes and the 30 ms build-up. Instantaneous Cartesian
% encoding at the end gives a 64-by-64 projection over a 12 mm FOV.
% Water parameters follow https://doi.org/10.1063/1.3116107.
%
%             [kspace,img,profile,psf,signal]=ddf_bore_image(npts,steps)
%
% Parameters:
%
%    npts    - 1-by-3 voxel counts covering the 8-by-8-by-10 mm box;
%              e.g. [48 48 60], first index along X
%    steps   - 1-by-3 positive substep counts for the two 2 ms lobes
%              and subsequent 30 ms build-up, e.g. [32 64 30]
%
% Outputs:
%
%    kspace  - centred 64-by-64-by-2 complex raw Fourier integrals
%              of Mx-i*My in A m^2 (axial and oblique preparation)
%    img     - centred images in integrated single-pulse signal units
%              per pixel, including the rectangular voxel aperture
%    profile - complex X profiles at Y=0 in those same units
%    psf     - central continuous point-response cut along X, at
%              (-256:256)/8 image pixels, normalised to unit peak
%    signal  - complex transverse magnetisation in A/m at the end of
%              each shot, prod(npts)-by-2, before Fourier encoding

function [kspace,img,profile,psf,signal]=ddf_bore_image(npts,steps)

% Specify a centred cylinder and an empty bore at X=2 mm
validateattributes(npts,{'double'},{'size',[1 3],'integer','positive','finite'});
validateattributes(steps,{'double'},{'size',[1 3],'integer','positive','finite'});
spacing=[8 8 10]*1e-3./npts;
axes_m=cell(1,3);
for n=1:3
    axes_m{n}=((1:npts(n))-(npts(n)+1)/2)*spacing(n);
end
[x,y,z]=ndgrid(axes_m{:});
mask=(x.^2+y.^2<0.004^2)&((x-0.002).^2+y.^2>0.001^2);

% Specify water at equilibrium and reflecting diffusion
gamma=spin('1H');
m0=1/(gamma*1.25663706127e-6*0.098);
equilibrium=m0*mask(:);
rate=ddf_diffusion(mask,spacing,2.15e-9);
symbol=ddf_bounded(npts,spacing);
field=@(m) ddf_cartesian(m,npts,symbol);
rhs=@(m) gamma*cross(m,field(m),2)+rate*m+...
         [-m(:,1)/0.250 -m(:,2)/0.250 (equilibrium-m(:,3))/0.288];

% Integrate the piecewise-constant voxel model over the readout moments
freq=(-32:31)'/0.012;
encode_x=exp(-2i*pi*freq*axes_m{1}).*sinc(freq*spacing(1));
encode_y=exp(-2i*pi*freq*axes_m{2}).*sinc(freq*spacing(2));
kspace=zeros(64,64,2);
signal=zeros(prod(npts),2);
directions=[0 0 1; 1/sqrt(2) 0 1/sqrt(2)];
duration=[0.002 0.002 0.030];

% Compare axial and oblique preparations with identical readout
for n=1:2

    % Wind the grating along the selected preparation axis
    coord=directions(n,1)*x(:)+directions(n,3)*z(:);
    magn=[zeros(numel(mask),2) equilibrium];
    for k=1:3

        % Apply a hard positive X pulse before each finite lobe
        if k<=2
            magn(:,2:3)=[-magn(:,3) magn(:,2)];
            angle=-2*pi*k*coord/(1e-3*steps(k));
        else
            angle=zeros(size(coord));
        end
        phase=exp(1i*angle/2);
        dt=duration(k)/steps(k);

        % Split exact gradient rotation around the midpoint torque step
        for j=1:steps(k)
            transverse=(magn(:,1)+1i*magn(:,2)).*phase;
            magn(:,1:2)=[real(transverse) imag(transverse)];
            magn=magn+dt*rhs(magn+(dt/2)*rhs(magn));
            transverse=(magn(:,1)+1i*magn(:,2)).*phase;
            magn(:,1:2)=[real(transverse) imag(transverse)];
        end

    end

    % Fourier encode the projected transverse moment in physical units
    signal(:,n)=magn(:,1)-1i*magn(:,2);
    density=reshape(signal(:,n),npts);
    kspace(:,:,n)=prod(spacing)*encode_x*sum(density,3)*encode_y.';

end

% Reconstruct on centred pixels and refer both shots to one hard X pulse
reference=1i*sum(equilibrium)*prod(spacing);
img=fftshift(fftshift(ifft2(ifftshift(ifftshift(kspace,1),2)),1),2)/reference;
profile=reshape(img(:,33,:),64,2);
psf=mean(exp(2i*pi*(-256:256)'*(-32:31)/(8*64)),2);

end


