% Test: kernel/ddf_bore | Bore phantom imaging
% Compare the bore image to a rotated direct-space Bloch calculation.
% Separately refine the spatial grid and the integration time steps.
% The coarse oracle tests the discrete forward map, not the curved
% continuum phantom. Syntax:
%
%                      result=test_ddf_bore()
%
% Outputs:
%
%    result - pre-encoding magnetisation and complex image comparisons

function result=test_ddf_bore()

% Find the example and restore the caller path on return
root_dir=fileparts(fileparts(fileparts(mfilename('fullpath'))));
old_path=path;
restore_path=onCleanup(@()path(old_path));
addpath(fullfile(root_dir,'examples','imaging','distant_dipolar'));

% Acquire both preparations with the FFT field and split Bloch solver
result=new_test_result('kernel/ddf_bore','Bore phantom imaging',...
                       'Independent discrete oracle, padding, and separate refinements.');
npts=[6 6 7];
steps=[32 64 30];
[~,img,~,~,signal]=ddf_bore_image(npts,steps);
spacing=[8 8 10]*1e-3./npts;
[a,b,c]=ndgrid(1:6,1:6,1:7);
x=(a-3.5)*spacing(1);
y=(b-3.5)*spacing(2);
z=(c-4)*spacing(3);
mask=(x.^2+y.^2<0.004^2)&(x.^2+(y-0.002).^2>0.001^2);

% Integrate each distinct displacement by Gaussian face quadrature
[dx,dy,dz]=ndgrid((0:5)*spacing(1),(0:5)*spacing(2),(0:6)*spacing(3));
weights=zeros(npts);
[nodes,w]=gaussleg(-0.5,0.5,23);
for j=1:numel(nodes)
    for k=1:numel(nodes)
        radial=(dx-spacing(1)*nodes(j)).^2+(dy-spacing(2)*nodes(k)).^2;
        for h=[-1 1]
            axial=dz-h*spacing(3)/2;
            weights=weights+h*w(j)*w(k)*prod(spacing(1:2))*...
                            axial./(4*pi*(radial+axial.^2).^(3/2));
        end
    end
end
weights(1)=weights(1)+1/3;

% Compare the open-boundary field at all targets, including the bore
dx=abs(a(:)-a(mask).');
dy=abs(b(:)-b(mask).');
dz=abs(c(:)-c(mask).');
weights=weights(sub2ind(npts,dx+1,dy+1,dz+1));
magn=[zeros(numel(mask),2) double(mask(:))];
mu0=1.25663706127e-6;
field=ddf_cartesian(magn,npts,ddf_bounded(npts,spacing))/mu0;
result=test_close(result,'direct field at all targets',field(:,3),sum(weights,2),1e-10,0,...
                  'quadrature integrates source cells even for vacuum targets');

% Double the vacuum box without changing source locations or spacing
padded=zeros([2*npts 3]);
padded(1:6,1:6,1:7,3)=double(mask);
pad_field=ddf_cartesian(reshape(padded,[],3),2*npts,ddf_bounded(2*npts,spacing))/mu0;
pad_field=reshape(pad_field(:,3),2*npts);
result=test_close(result,'double vacuum padding',pad_field(1:6,1:6,1:7),reshape(field(:,3),npts),...
                  1e-12,0,'open convolution excludes fields from periodic copies');

% Form the dense field matrix and independent nearest-face exchange
coupling=weights(mask(:),:);
dx=dx(mask(:),:);
dy=dy(mask(:),:);
dz=dz(mask(:),:);
faces=(dx+dy+dz)==1;
hops=2.15e-9*sparse(faces.*(dx/spacing(1)^2+dy/spacing(2)^2+dz/spacing(3)^2));
count=nnz(mask);
rate=hops-spdiags(sum(hops,2),0,count,count);

% Evolve normalised Bloch vectors with the gradients inside the ODE
m0=1/(spin('1H')*1.25663706127e-6*0.098);
opts=odeset('RelTol',1e-8,'AbsTol',1e-10);
duration=[0.002 0.002 0.030];
for n=1:2

    % Rotate the original gradient with the cylinder around the Z axis
    coord=z(mask);
    if n==2
        coord=(y(mask)+z(mask))/sqrt(2);
    end
    magn=[zeros(count,2) ones(count,1)];
    for k=1:3
        omega=zeros(count,1);
        if k<=2
            magn(:,2:3)=[-magn(:,3) magn(:,2)];
            omega=2*pi*k*coord/(1e-3*duration(k));
        end
        fun=@(~,m) torque(m,omega,coupling,rate);
        [~,trajectory]=ode45(fun,[0 duration(k)],magn(:),opts);
        magn=reshape(trajectory(end,:),[],3);
    end

    % Undo the spatial rotation before comparing complex voxel signals
    voxels=zeros(npts);
    voxels(mask)=magn(:,1)-1i*magn(:,2);
    voxels=flip(permute(voxels,[2 1 3]),2);
    residual=norm(signal(:,n)/m0-voxels(:))/norm(voxels(:));
    result=test_close(result,'rotated voxel signal',residual,0,0.001,0,...
                      'rotation and direct integration preserve complex magnetisation');

    % Integrate each rectangular voxel directly at the same readout moments
    [kx,ky]=ndgrid((-32:31)/0.012);
    moments=zeros(64);
    for j=1:numel(voxels)
        moments=moments+voxels(j)*exp(-2i*pi*(kx*x(j)+ky*y(j)));
    end
    moments=moments.*sinc(kx*spacing(1)).*sinc(ky*spacing(2));
    oracle=fftshift(ifft2(ifftshift(moments)))/(1i*count);
    residual=norm(img(:,:,n)-oracle,'fro')/norm(oracle,'fro');
    result=test_close(result,'complex image',residual,0,0.01,0,...
                      'the independent reconstruction uses the same finite aperture');

end

% Compare successive spatial refinements at fixed time steps and readout bandwidth
steps=[16 32 30];
grids=[24 24 30; 32 32 40; 48 48 60];
coarse=ddf_bore_image(grids(1,:),steps);
middle=ddf_bore_image(grids(2,:),steps);
fine=ddf_bore_image(grids(3,:),steps);

% Halve each time step without changing the fine spatial grid
half_step=ddf_bore_image(grids(3,:),2*steps);
for n=1:2

    % Parseval equates the relative raw complex image and k-space differences
    scale=norm(fine(:,:,n),'fro');
    coarse_change=norm(coarse(:,:,n)-middle(:,:,n),'fro')/scale;
    fine_change=norm(middle(:,:,n)-fine(:,:,n),'fro')/scale;
    time_change=norm(fine(:,:,n)-half_step(:,:,n),'fro')/norm(half_step(:,:,n),'fro');
    result=test_true(result,sprintf('spatial refinement shot %d',n),fine_change<coarse_change,...
                     sprintf('successive relative changes %.6g, %.6g',coarse_change,fine_change));
    result=test_close(result,sprintf('time refinement shot %d',n),time_change,0,0.001,0,...
                      'step halving consumes at most a tenth of the one-percent image budget');

end

end

% Evaluate the unsplit direct-summation Bloch equations
function derivative=torque(state,omega,coupling,rate)

% Form the local angular frequency including the finite gradient
magn=reshape(state,[],3);
angular=(coupling*magn).*[-0.5 -0.5 1]/0.098;
angular(:,3)=angular(:,3)+omega;

% Add dipolar precession, reflecting transport, and equilibrium recovery
derivative=cross(magn,angular,2)+rate*magn+...
           [-magn(:,1)/0.250 -magn(:,2)/0.250 (1-magn(:,3))/0.288];
derivative=derivative(:);

end


