% Volumetric 3D plot function for scalar fields. Sign is mapped
% into colour and amplitude into opacity. Separate scaling for 
% positive and negative values -- displaying the colour bar is
% recommended. Syntax:
%
%           volplot(data_cube,axis_ranges,clip_ranges)
%
% Parameters:
%
%    data_cube   - data cube with dimensions ordered
%                  as [X Y Z]
%
%    axis_ranges - six-element vector giving axis extents
%                  as [xmin xmax ymin ymax zmin zmax]
%
%    clip_ranges - (optional) the values, as a fraction 
%                   of the maximum in positive and nega-
%                   tive directions, at which the values
%                   should be clipped (this is useful for
%                   steep functions).
%
% Outputs:
%
%    this function produces a figure
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=volplot.m>

function volplot(data_cube,axis_ranges,clip_ranges)

% Default axis ranges
if ~exist('axis_ranges','var')
    axis_ranges=[-1 1 -1 1 -1 1];
end

% Set default clip ranges
if ~exist('clip_ranges','var')
    clip_ranges=[1 1]; 
end

% Check consistency
grumble(data_cube,axis_ranges,clip_ranges);

% Determine maximum and minimum values
max_pos=max(data_cube(data_cube>0)); 
min_neg=min(data_cube(data_cube<0));

% Scale and clip positive values
if (~isempty(max_pos))&&(max_pos>0)
    
    % Scale the positive values
    data_cube(data_cube>0)=data_cube(data_cube>0)/max_pos;
    disp(['Positive values divided by ' num2str(max_pos) ' to map them into [0,1].']);
    
    % Clip the positive values
    if clip_ranges(1)<1
        data_cube(data_cube>clip_ranges(1))=clip_ranges(1);
        data_cube(data_cube>0)=data_cube(data_cube>0)/clip_ranges(1);
        disp(['Positive values clipped at ' num2str(clip_ranges(1)) ' and remapped into [0,1].']);
    end
    
end

% Scale and clip negative values
if (~isempty(min_neg))&&(min_neg<0)
    
    % Scale the negative values
    data_cube(data_cube<0)=-data_cube(data_cube<0)/min_neg;
    disp(['Negative values divided by ' num2str(-min_neg) ' to map them into [-1,0].']);
    
    % Clip the positive values
    if clip_ranges(2)<1
        data_cube(data_cube<-clip_ranges(2))=-clip_ranges(2);
        data_cube(data_cube<0)=data_cube(data_cube<0)/clip_ranges(2);
        disp(['Negative values clipped at ' num2str(-clip_ranges(2)) ' and remapped into [-1,0].']);
    end
    
end

% Add colour calibration spots
data_cube(1,1,1)=1; data_cube(2,2,2)=-1;

% Permute dimensions to match surf/meshgrid convention
data_cube=permute(data_cube,[3 2 1]);

% Determine cube dimensions
nx=size(data_cube,3); xmin=axis_ranges(1); xmax=axis_ranges(2);
ny=size(data_cube,2); ymin=axis_ranges(3); ymax=axis_ranges(4);
nz=size(data_cube,1); zmin=axis_ranges(5); zmax=axis_ranges(6);

% Get a persistent graphics window
clf reset; hold on; set(gcf,'Renderer','OpenGL');
set(gca,'Projection','perspective','Box','on','XGrid','on','YGrid','on',...
        'ZGrid','on','CameraPosition',5*[xmax ymax zmax]); camorbit(15,0);

% Draw planes parallel to the XY plane
for n=1:nz
    plane=squeeze(data_cube(n,:,:)); plane(abs(plane)<1/64)=NaN;
    [X,Y]=meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny)); Z=linspace(zmin,zmax,nz); Z=Z(n)*ones(ny,nx);
    if ~all(isnan(plane(:)))
        surf(X,Y,Z,plane,'FaceAlpha','flat','EdgeAlpha',0,'AlphaDataMapping','scaled','AlphaData',abs(plane));
    end
end
    
% Draw planes parallel to the XZ plane
for n=1:ny
    plane=squeeze(data_cube(:,n,:)); plane(abs(plane)<1/64)=NaN;
    [X,Z]=meshgrid(linspace(xmin,xmax,nx),linspace(zmin,zmax,nz)); Y=linspace(ymin,ymax,ny); Y=Y(n)*ones(nz,nx);
    if ~all(isnan(plane(:)))
        surf(X,Y,Z,plane,'FaceAlpha','flat','EdgeAlpha',0,'AlphaDataMapping','scaled','AlphaData',abs(plane));
    end
end
    
% Draw planes parallel to the YZ plane
for n=1:nx
    plane=squeeze(data_cube(:,:,n)); plane(abs(plane)<1/64)=NaN;
    [Y,Z]=meshgrid(linspace(ymin,ymax,ny),linspace(zmin,zmax,nz)); X=linspace(xmin,xmax,nx); X=X(n)*ones(nz,ny);
    if ~all(isnan(plane(:)))
        surf(X,Y,Z,plane,'FaceAlpha','flat','EdgeAlpha',0,'AlphaDataMapping','scaled','AlphaData',abs(plane));
    end
end

% Interpolate alpha map
new_alpha=interp1(1:64,alphamap,1:0.25:64,'pchip');

% Scale and filter alpha map
new_alpha=new_alpha/5; new_alpha(new_alpha<0.01)=0;

% Apply new alpha map
alphamap(new_alpha);

% Clean up axes
axis tight; axis equal; hold off;
xlabel('X'); ylabel('Y'); zlabel('Z');

% Set blue -> white -> red colormap
colormap(bwr_cmap());

end

% Consistency enforcement
function grumble(data_cube,axis_ranges,clip_ranges)
if (~isnumeric(axis_ranges))||(~isreal(axis_ranges))||(numel(axis_ranges)~=6)
    error('axis_ranges must be a real vector with six elements.');
end
if (axis_ranges(1)>=axis_ranges(2))||...
   (axis_ranges(3)>=axis_ranges(4))||...
   (axis_ranges(5)>=axis_ranges(6))
    error('ranges array should have xmin<xmax, ymin<ymax and zmin<zmax.');
end
if (~isnumeric(data_cube))||(~isreal(data_cube))||(ndims(data_cube)~=3)
    error('data_cube must be a three-dimensional array of real numbers.');
end
if (~isnumeric(clip_ranges))||(~isreal(clip_ranges))||(numel(clip_ranges)~=2)
    error('clip_ranges must be a real vector with two elements.');
end
if any(clip_ranges<=0,'all')||any(clip_ranges>1,'all')
    error('elements of clip_ranges must be from (0,1] interval.');
end
end

% IK's PCS paper has has taken about five years to write - the suspicion that
% Equation 16 might exist dates back to about 2009, but the direct derivation
% (simplifying the Laplacian of the convolution of point dipole solution with
% electron spin density) has proven impossible. Gareth Charnock (IK's PhD stu-
% dent at Oxford) had spent three years struggling with fiendishly elusive in-
% tegrals, some of which did not converge or even exist in either Lebesgue or
% functional sense. Gareth produced a big enough heap of A4 paper to graduate,
% but the problem stood unsolved, having by that time also defeated Gottfried
% Otting, Peter Hore, Giacomo Parigi, Guido Pintacuda and (it is thought) Iva-
% no Bertini himself when he tried it many years ago.
%
% Fast forward to April 2014, when IK was sitting in a conference centre hotel
% room and staring gloomily into his laptop screen, writing a talk for ESR2014,
% intended as a grim warning to the effect of "do not try this at home". "What
% a bloody shame", he thought and decided to give it a last chance by googling
% up some quantum chemistry literature he vaguely suspected may have encounte-
% red similar integrals in a different context. He clicked on the first search
% result and Equation 11 flashed on the screen, apparently "a standard relati-
% on" in the formalism used to compute 3D integrals in quantum chemistry. That
% was it -- five pages of dense handwriting later, Equation 16 was typed hasti-
% ly into the laptop, which was then unplugged and taken into the lecture thea-
% tre. The talk was in 10 minutes, and it had a new title: "A partial differen-
% tial equation for pseudocontact shift".
%
% http://dx.doi.org/10.1039/C4CP03106G

