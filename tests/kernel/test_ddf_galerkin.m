% BlochDDFGalerkin midpoint-field and reflecting-sphere comparisons.
% Reference: Bouchard, J. Chem. Phys. 164, 144121 (2026). Syntax:
%
%                     result=test_ddf_galerkin()
%
% Outputs:
%
%    result - regression checks against direct fields, archived decay
%             traces, and the continuum Neumann eigenmode

function result=test_ddf_galerkin()

% Sample a trilinear magnetisation at the upstream centre quadrature points
result=new_test_result('kernel/ddf_galerkin','Bounded Galerkin comparisons',...
                       'Softened midpoint fields and archived sphere decays.');
spacing=[0.3 0.4 0.2];
npts=[3 2 4];
[x,y,z]=ndgrid(((1:3)-0.5)*spacing(1),((1:2)-0.5)*spacing(2),((1:4)-0.5)*spacing(3));
magn=[x(:)+y(:).*z(:) z(:)-x(:).*y(:) 1+x(:).*z(:)];

% Compare the corner vector with the upstream all-pairs evaluation
radii=[0.13 0.31];
terms={'omit','regularised'};
oracle=[ 0.2155816309073895   0.8175897386134489 -3.410611583222289;...
        -0.7129572857972075   0.4352501846762616  7.677235480956138;...
         0.03177422107908884  0.3620022140468250 -0.9763982913433374;...
        -0.03670283575015490  0.3338057788818423 -0.1587016715588386];
for n=1:2
    for k=1:2
        symbol=ddf_regularised(npts,spacing,radii(n),terms{k});
        field=ddf_cartesian(magn,npts,symbol)*(-4*pi/1.25663706127e-6);
        result=test_close(result,'bounded corner field',field(1,:),oracle(2*n-2+k,:),...
                          1e-12,0,'all three secular components match at the wall');
    end
end

% Reproduce the archived voxel decay errors and refine to the continuum
counts=[10 12];
radii=[0.25 0.30];
archived=[0.02588002956709021 0.008257106035406914];
for n=1:2

    % Match the archived time grid, sphere, diffusion, and relaxation
    [trace,coarse]=ddf_sphere_decay(counts(n),radii(n));
    decay_error=norm(trace(:,2)-trace(:,3))/norm(trace(:,3));
    result=test_close(result,'archived FD trace error',decay_error,archived(n),...
                      1e-12,0,'relative trace error matches the article archive');

    % Demand local as well as projected convergence on a refined mask
    [trace,fine]=ddf_sphere_decay(4*counts(n),radii(n));
    decay_error=norm(trace(:,2)-trace(:,3))/norm(trace(:,3));
    result=test_true(result,'reflecting wall refinement',fine<coarse&&fine<0.01&&decay_error<0.01,...
                     'relative local and projected errors both fall below one percent');

end

end


