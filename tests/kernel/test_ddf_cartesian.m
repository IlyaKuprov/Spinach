% Cartesian dipolar field tests on anisotropic grids. Syntax:
%
%                     result=test_ddf_cartesian()
%
% Outputs:
%
%    result - oblique mode, direct-space, padding, and torque checks

function result=test_ddf_cartesian()

% Predict the field of an oblique Fourier mode in a rectangular box
result=new_test_result('kernel/ddf_cartesian','Cartesian dipolar fields',...
                       'Physical field checks for isolated and periodic boundaries.');
npts=[5 4 3];
spacing=[1.1 1.5 0.9]*1e-3;
[x,y,z]=ndgrid(0:4,0:3,0:2);
phase=2*pi*(x/5-y/4+z/3);
wave=cos(phase(:))*[1 2 3];
uniform=[3 2 1];
magn=wave+uniform;
mu0=1.25663706127e-6;
k=[1 -1 1]./(npts.*spacing);
symbol=ddf_periodic(npts,spacing,0.7);
field=ddf_cartesian(magn,npts,symbol);
target=((1/3-k(3)^2/sum(k.^2))*wave+(1/3-0.7)*uniform).*[-0.5 -0.5 1];
result=test_close(result,'oblique mode',field/mu0,target,1e-12,0,...
                  'the reciprocal axes scale inversely with the physical box');

% Compute all pairwise fields by quadrature of the two charged Z faces
dx=(x(:)-x(:).')*spacing(1);
dy=(y(:)-y(:).')*spacing(2);
dz=(z(:)-z(:).')*spacing(3);
[nodes,weights]=gaussleg(-0.5,0.5,31);
direct=eye(prod(npts))/3;
for a=1:numel(nodes)
    for b=1:numel(nodes)
        radial=(dx-nodes(a)*spacing(1)).^2+(dy-nodes(b)*spacing(2)).^2;
        for c=[-1 1]
            axial=dz-c*spacing(3)/2;
            direct=direct+c*weights(a)*weights(b)*prod(spacing(1:2))*...
                          axial./(4*pi*(radial+axial.^2).^(3/2));
        end
    end
end
symbol=ddf_bounded(npts,spacing);
field=ddf_cartesian(magn,npts,symbol);
target=(direct*magn).*[-0.5 -0.5 1];
result=test_close(result,'open-boundary sum',field/mu0,target,1e-11,0,...
                  'independent face quadrature includes the noncubic self cell');

% Add empty exterior voxels without changing the physical source
padded=zeros([2*npts 3]);
padded(1:npts(1),1:npts(2),1:npts(3),:)=reshape(magn,[npts 3]);
wide=ddf_bounded(2*npts,spacing);
padded=ddf_cartesian(reshape(padded,[],3),2*npts,wide);
padded=reshape(padded,[2*npts 3]);
inner=padded(1:npts(1),1:npts(2),1:npts(3),:);
result=test_close(result,'vacuum padding',reshape(inner,[],3)/mu0,target,1e-11,0,...
                  'vacuum does not create a periodic image');

% Evolve the same nonuniform initial magnetisation with either field
seed=magn;
reference=magn;
for n=1:20
    seed=seed+0.005*cross(seed,ddf_cartesian(seed,npts,symbol)/mu0,2);
    reference=reference+0.005*cross(reference,(direct*reference).*[-0.5 -0.5 1],2);
end
signal=seed(:,1)-1i*seed(:,2);
oracle=reference(:,1)-1i*reference(:,2);
result=test_close(result,'complex voxel signal',signal,oracle,1e-11,0,...
                  'the quadrature field produces the same nonlinear torque history');

end


