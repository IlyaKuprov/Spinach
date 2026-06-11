% Interpolates measured gradient coil field maps at spin coordinates.
% Syntax:
%
%          fields=gradient_fields_at_spins(field_maps,coordinates,currents)
%          fields=gradient_fields_at_spins(field_maps,coordinates,currents,map_names)
%
% Parameters:
%
%    field_maps - either a structure containing coil field maps, or a
%                 path to a MAT file containing them. Each map must be
%                 a finite numeric array with columns x y z Bx By Bz,
%                 with coordinates in metres and fields in tesla per
%                 ampere
%
%  coordinates - spin coordinates, either an N-by-3 array in metres or
%                a cell array of Spinach coordinates in Angstrom
%
%     currents - coil currents, amperes
%
%    map_names - optional cell array of field names in field_maps; the
%                default is {'B_X_Coils','B_Y_Coils','B_Z_Coils'}
%
% Outputs:
%
%       fields - N-by-3 array of interpolated Bx By Bz fields, tesla
%
% Note: the interpolant is linear inside the measured map and nearest
%       neighbour outside it.
%
% a.arnab@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=gradient_fields_at_spins.m>

function fields=gradient_fields_at_spins(field_maps,coordinates,currents,map_names)

% Default map names
if nargin<4
    map_names={'B_X_Coils','B_Y_Coils','B_Z_Coils'};
end
% Check consistency
grumble(field_maps,coordinates,currents,map_names);

% Load maps if necessary
if ischar(field_maps)||isstring(field_maps)
    field_maps=load(field_maps,map_names{:});
end

% Convert Spinach coordinate cells from Angstrom to metres
if iscell(coordinates)
    coordinates=cell2mat(coordinates)*1e-10;
end

% Accumulate coil contributions
fields=zeros(size(coordinates));
for n=1:numel(map_names)
    map=field_maps.(map_names{n});
    validate_field_map(map,map_names{n});
    fields=fields+currents(n)*interpolate_field_map(map,coordinates);
end

end

% Field map interpolation
function field=interpolate_field_map(map,positions)

% Preallocate output
field=zeros(size(positions));

% Interpolate Cartesian field components
for component=1:3
    interpolant=scatteredInterpolant(map(:,1),map(:,2),map(:,3), ...
        map(:,3+component),'linear','nearest');
    field(:,component)=interpolant(positions(:,1),positions(:,2), ...
                                   positions(:,3));
end

end

% Field map validation
function validate_field_map(map,label)

if size(map,2)~=6 || ~isnumeric(map) || any(~isfinite(map(:)))
    error([label ' must be a finite numeric array with columns x y z Bx By Bz.']);
end

if numel(unique(map(:,1)))<2 || numel(unique(map(:,2)))<2 || ...
   numel(unique(map(:,3)))<2
    error([label ' must contain a three-dimensional coordinate grid.']);
end

end

% Consistency enforcement
function grumble(field_maps,coordinates,currents,map_names)

if ~(isstruct(field_maps)||ischar(field_maps)||isstring(field_maps))
    error('field_maps must be a structure or a MAT file path.');
end
if isstring(field_maps)&&(~isscalar(field_maps))
    error('field_maps string must be scalar.');
end
if ~(iscell(coordinates)||(isnumeric(coordinates)&&isreal(coordinates)))
    error('coordinates must be a numeric N-by-3 array or a Spinach coordinate cell array.');
end
if iscell(coordinates)
    try
        coordinates=cell2mat(coordinates);
    catch
        error('coordinate cell array must contain numeric three-element vectors.');
    end
end
if size(coordinates,2)~=3 || any(~isfinite(coordinates(:)))
    error('coordinates must contain one finite x y z row per spin.');
end
if (~isnumeric(currents))||(~isreal(currents))||...
   (numel(currents)~=numel(map_names))||any(~isfinite(currents(:)))
    error('currents must be a finite real vector with one element per field map.');
end
if ~iscell(map_names)||isempty(map_names)||...
   any(~cellfun(@(x)ischar(x)&&isrow(x)&&(~isempty(x)),map_names))
    error('map_names must be a non-empty cell array of character vectors.');
end
if isstruct(field_maps)
    for n=1:numel(map_names)
        if ~isfield(field_maps,map_names{n})
            error(['missing field map ' map_names{n} '.']);
        end
    end
end

end
