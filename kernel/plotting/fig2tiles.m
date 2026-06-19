% Combines Matlab figure files into a single tiled figure. Syntax:
%
%              [fig_obj,tile_obj]=fig2tiles(fig_files)
%
% Parameters:
%
%    fig_files - cell array of character strings containing
%                Matlab *.fig file names
%
% Outputs:
%
%    fig_obj  - handle of the new Matlab figure
%
%    tile_obj - handle of the tiled layout object
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=fig2tiles.m>

function [fig_obj,tile_obj]=fig2tiles(fig_files)

% Check consistency
grumble(fig_files);

% Create the new figure
fig_obj=kfigure();

% Create a loose tiled layout
tile_obj=tiledlayout(fig_obj,'flow','TileSpacing','loose',...
                    'Padding','loose');

% Import saved figures
for n=1:numel(fig_files)

    % Open source figure without displaying it
    src_fig=openfig(fig_files{n},'new','invisible');

    try

        % Collect source figure children
        src_obj=get(src_fig,'Children');

        % Reject empty figures
        if isempty(src_obj)
            error('source figure contains no graphics objects.');
        end

        % Copy source graphics into the layout
        copy_obj=copyobj(src_obj,tile_obj);

        % Put copied tile-capable graphics into the next tile
        tile_found=false;
        for k=1:numel(copy_obj)
            if isprop(copy_obj(k),'Type')
                obj_type=get(copy_obj(k),'Type');
            else
                obj_type='';
            end
            if isprop(copy_obj(k),'Layout')&&...
               (~strcmp(obj_type,'legend'))&&...
               (~strcmp(obj_type,'colorbar'))
                layout_obj=copy_obj(k).Layout;
                if (~isempty(layout_obj))&&isprop(layout_obj,'Tile')
                    copy_obj(k).Layout.Tile=n;
                    tile_found=true;
                end
            end
        end

        % Reject unsupported source figures
        if ~tile_found
            error('source figure contains no tile-capable graphics objects.');
        end

        % Preserve axes colour maps
        src_axes=findobj(src_obj,'Type','axes');
        copy_axes=findobj(copy_obj,'Type','axes');
        if numel(src_axes)==numel(copy_axes)
            for k=1:numel(copy_axes)
                copy_axes(k).Colormap=src_axes(k).Colormap;
            end
        end

    catch err

        % Delete the invisible source figure
        delete(src_fig);

        % Rethrow the error
        rethrow(err);

    end

    % Delete the invisible source figure
    delete(src_fig);

end

% Update layout geometry
drawnow;

end

% Consistency enforcement
function grumble(fig_files)
if (~iscell(fig_files))||isempty(fig_files)||(~isvector(fig_files))
    error('fig_files must be a non-empty cell array of file names.');
end
for n=1:numel(fig_files)
    if (~ischar(fig_files{n}))||(~isrow(fig_files{n}))||isempty(fig_files{n})
        error('elements of fig_files must be character strings.');
    end
    if (numel(fig_files{n})<5)||(~strcmpi(fig_files{n}(end-3:end),'.fig'))
        error('elements of fig_files must be *.fig file names.');
    end
    if exist(fig_files{n},'file')~=2
        error('all files listed in fig_files must exist.');
    end
end
end


