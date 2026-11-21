% Test: kernel/ddf_sphere | Reflecting sphere
% Reflecting sphere diffusion and spherical cavity field checks. Syntax:
%
%                     result=test_ddf_sphere()
%
% Outputs:
%
%    result - regression checks for grid refinement and conservation

function result=test_ddf_sphere()

% Compare the radial Neumann mode with its continuum decay
result=new_test_result('kernel/ddf_sphere','Reflecting sphere',...
                       'Diffusion and cavity fields on refined voxel spheres.');
root=fzero(@(q) sin(q)-q*cos(q),[pi 1.5*pi]);
errors=zeros(2,1);
for n=1:2

    % Sample a sphere of unit radius on successively finer cubic grids
    count=12*n;
    pitch=2/count;
    axis_vals=((1:count)-(count+1)/2)*pitch;
    [x,y,z]=ndgrid(axis_vals);
    radius=sqrt(x.^2+y.^2+z.^2);
    mask=radius<1;
    mode=zeros(size(mask));
    mode(mask)=sin(root*radius(mask))./(root*radius(mask));

    % Propagate the eigenmode with stable explicit conservative fluxes
    rate=ddf_diffusion(mask,pitch*[1 1 1],1);
    duration=1/root^2;
    steps=ceil(12*duration/pitch^2);
    density=mode(:);
    for k=1:steps
        density=density+(duration/steps)*(rate*density);
    end
    errors(n)=norm(density-exp(-1)*mode(:))/norm(mode(:));
    result=test_close(result,'mass balance',sum(density),sum(mode(:)),1e-11,0,...
                      'there is no flux through the voxelised sphere boundary');
    result=test_close(result,'no leakage',norm(density(~mask)),0,0,0,...
                      'exterior and interior cells share no open faces');

    % Measure the interior field away from the staircased wall
    magn=[zeros(numel(mask),2) double(mask(:))];
    symbol=ddf_bounded([count count count],pitch*[1 1 1]);
    field=ddf_cartesian(magn,[count count count],symbol);
    interior=field(radius(:)<0.5,3)/1.25663706127e-6;
    result=test_close(result,'sphere bulk field',max(abs(interior)),0,0.01,0,...
                      'a uniform continuum sphere has zero cavity field');

end

% Demand improvement of the full mode, not just its spatial integral
result=test_true(result,'sphere refinement',errors(2)<errors(1)&&errors(2)<0.04,...
                 sprintf('relative mode errors: %.4g, %.4g',errors));

end


