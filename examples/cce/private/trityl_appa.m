% High-field analytical pair product approximation (APPA), Eqs. (8)-(11)
% of doi:10.5194/mr-6-15-2025, on the nine-ring Ci quadrature used there.
% This is the projected CCE-2 limit with secular hyperfine coupling and
% an infinite-temperature proton bath; every one-proton factor is unity.
%
% Syntax:
%
%   [answer,diagnostics]=trityl_appa(coords,hyperfine,tunnel,times)
%
% Parameters:
%
%   coords     - N-by-3 proton coordinates, Angstrom
%   hyperfine  - 3-by-3-by-N electron-proton tensors, rad/s
%   tunnel     - N-by-N symmetric matrix of tunnel splittings, rad/s
%   times      - column of nonnegative total Hahn echo times, seconds
%
% Outputs:
%
%   answer      - normalised Hahn echo column
%   diagnostics - times, directions (one row per node), weights (column),
%                 and echo traces (one column per node)
%

function [answer,diagnostics]=trityl_appa(coords,hyperfine,tunnel,times)

% Check the physical input dimensions and units
grumble(coords,hyperfine,tunnel,times);

% Build the nine-ring hemisphere quadrature, integrating each annular cell
knots=9;
step=pi/(2*(knots-1));
directions=[0 0 1];
weights=1-cos(step/2);
for n=1:knots-1

    % Use four additional azimuths for each step away from the pole
    theta=n*step;
    azimuth=2*pi*(0:4*n-1).'/(4*n);
    directions=[directions; sin(theta)*cos(azimuth),...
                            sin(theta)*sin(azimuth),...
                            cos(theta)*ones(4*n,1)]; %#ok<AGROW>
    lower=theta-step/2;
    upper=min(theta+step/2,pi/2);
    weights=[weights; repmat((cos(lower)-cos(upper))/(4*n),4*n,1)]; %#ok<AGROW>

end

% Compute the dipolar tensor for each proton pair
pairs=nchoosek(1:size(coords,1),2);
dipolar=zeros(3,3,size(pairs,1));
for n=1:size(pairs,1)

    % Spinach returns Cartesian dipolar tensors in angular frequency units
    [~,~,~,~,dipolar(:,:,n)]=xyz2dd(coords(pairs(n,1),:),...
                                   coords(pairs(n,2),:),'1H','1H');

end
exchange=4*tunnel(sub2ind(size(tunnel),pairs(:,1),pairs(:,2)))/3;

% Average completed cluster products, not the individual pair signals
values=zeros(numel(times),numel(weights));
for n=1:numel(weights)

    % Contract each tensor with the laboratory field direction
    field=directions(n,:);
    projection=kron(field,field);
    secular_hf=projection*reshape(hyperfine,9,[]);
    pair_dd=(projection*reshape(dipolar,9,[])).'+exchange;
    difference=secular_hf(pairs(:,1)).'-secular_hf(pairs(:,2)).';
    frequency=hypot(difference,pair_dd);

    % A pair with zero zero-quantum splitting contributes exactly one
    active=frequency>0;
    depth=4*((difference(active)./frequency(active)).*...
             (pair_dd(active)./frequency(active))).^2;
    values(:,n)=prod(1-depth.*sin((frequency(active)/8)*times.').^4,1).';

end
answer=values*weights;
diagnostics=struct('times',times,'directions',directions,...
                   'weights',weights,'values',values);

end

% Check consistency
function grumble(coords,hyperfine,tunnel,times)
n=size(coords,1);
if ~isnumeric(coords)||~isreal(coords)||~ismatrix(coords)||...
   size(coords,2)~=3||n<2||any(~isfinite(coords),'all')
    error('coords must be a finite N-by-3 real matrix with at least two rows.');
end
if ~isnumeric(hyperfine)||~isreal(hyperfine)||...
   ~isequal(size(hyperfine),[3 3 n])||any(~isfinite(hyperfine),'all')
    error('hyperfine must be a finite real 3-by-3-by-N array.');
end
if ~isnumeric(tunnel)||~isreal(tunnel)||...
   ~isequal(size(tunnel),[n n])||any(~isfinite(tunnel),'all')||~issymmetric(tunnel)
    error('tunnel must be a finite real symmetric N-by-N matrix.');
end
if ~isnumeric(times)||~isreal(times)||~iscolumn(times)||...
   isempty(times)||any(~isfinite(times)|times<0)
    error('times must be a nonempty column of finite nonnegative values.');
end
end


