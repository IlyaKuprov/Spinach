% Interpolates measured magnetic field maps at spin coordinates.
% Syntax:
%
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
%    map_names - cell array of field names in field_maps
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

% Check consistency
grumble(field_maps,coordinates,currents,map_names);

% Load maps if necessary
if ischar(field_maps)||isstring(field_maps)
    field_maps=load(field_maps,map_names{:});
end

% Convert Spinach coordinate cells from Angstrom to metres
if iscell(coordinates)
    coordinates=cell2mat(coordinates(:))*1e-10;
end

% Preallocate the field array
fields=zeros(size(coordinates));

% Accumulate coil contributions
for n=1:numel(map_names)

    % Extract the current field map
    field_map=field_maps.(map_names{n});

    % Check field map consistency
    if (~isnumeric(field_map))||(~isreal(field_map))||...
       (size(field_map,2)~=6)||any(~isfinite(field_map(:)))
        error([map_names{n} ' must be a finite real array with columns x y z Bx By Bz.']);
    end

    % Check field map dimensionality
    if (numel(unique(field_map(:,1)))<2)||...
       (numel(unique(field_map(:,2)))<2)||...
       (numel(unique(field_map(:,3)))<2)
        error([map_names{n} ' must contain a three-dimensional coordinate grid.']);
    end

    % Accumulate interpolated Cartesian field components
    for k=1:3

        % Build a scattered field interpolant
        field_int=scatteredInterpolant(field_map(:,1),field_map(:,2),...
                                       field_map(:,3),field_map(:,3+k),...
                                       'linear','nearest');

        % Add the weighted field component
        fields(:,k)=fields(:,k)+currents(n)*...
                    field_int(coordinates(:,1),coordinates(:,2),...
                              coordinates(:,3));
    end
end

end

% Consistency enforcement
function grumble(field_maps,coordinates,currents,map_names)

if (~iscell(map_names))||isempty(map_names)||...
   any(~cellfun(@(x)ischar(x)&&isrow(x)&&(~isempty(x)),map_names))
    error('map_names must be a non-empty cell array of character vectors.');
end
if ~(isstruct(field_maps)||ischar(field_maps)||isstring(field_maps))
    error('field_maps must be a structure or a MAT file path.');
end
if ischar(field_maps)&&(~isrow(field_maps))
    error('field_maps character array must be a row.');
end
if isstring(field_maps)&&(~isscalar(field_maps))
    error('field_maps string must be scalar.');
end
if (ischar(field_maps)||isstring(field_maps))&&(~isfile(field_maps))
    error('field_maps must be an existing MAT file path.');
end
if isstruct(field_maps)
    for n=1:numel(map_names)
        if ~isfield(field_maps,map_names{n})
            error(['missing field map ' map_names{n} '.']);
        end
    end
end
if ~(iscell(coordinates)||(isnumeric(coordinates)&&isreal(coordinates)))
    error('coordinates must be a numeric N-by-3 array or a Spinach coordinate cell array.');
end
if iscell(coordinates)
    for n=1:numel(coordinates)
        if (~isnumeric(coordinates{n}))||(~isreal(coordinates{n}))||...
           (~isequal(size(coordinates{n}),[1 3]))
            error('coordinates cell array must contain real 1x3 vectors.');
        end
    end
    coordinates=cell2mat(coordinates(:));
end
if (size(coordinates,2)~=3)||any(~isfinite(coordinates(:)))
    error('coordinates must contain one finite x y z row per spin.');
end
if (~isnumeric(currents))||(~isreal(currents))||...
   (numel(currents)~=numel(map_names))||any(~isfinite(currents(:)))
    error('currents must be a finite real vector with one element per field map.');
end

end

