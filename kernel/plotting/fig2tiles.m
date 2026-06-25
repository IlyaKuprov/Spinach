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

% Get the tile grid size
tile_rows=size(fig_files,1);
tile_cols=size(fig_files,2);

% Map cell array entries to tiled layout positions
tile_nums=zeros(numel(fig_files),1);
for n=1:numel(fig_files)
    [row_num,col_num]=ind2sub(size(fig_files),n);
    tile_nums(n)=(row_num-1)*tile_cols+col_num;
end

% Create the new figure
fig_obj=kfigure();

% Create a loose tiled layout
tile_obj=tiledlayout(fig_obj,tile_rows,tile_cols,'TileSpacing','loose',...
                    'Padding','loose');

% Allocate all outer tiles
tile_axes=gobjects(numel(fig_files),1);
tile_pos=zeros(numel(fig_files),4);
for n=1:numel(fig_files)
    tile_axes(n)=nexttile(tile_obj,tile_nums(n));
    tile_axes(n).Visible='off';
end



% Record outer tile positions
drawnow;
for n=1:numel(fig_files)
    tile_pos(n,:)=tile_axes(n).OuterPosition;
end

% Initialise overlay panel tracking
panel_axes=gobjects(numel(fig_files),1);
panel_objs=gobjects(numel(fig_files),1);
panel_count=0;

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

        % Count source tiled layouts
        tile_count=0;
        for k=1:numel(src_obj)
            if strcmp(get(src_obj(k),'Type'),'tiledlayout')
                tile_count=tile_count+1;
            end
        end

        % Collect source tiled layouts
        src_tiles=gobjects(tile_count,1);
        tile_count=0;
        for k=1:numel(src_obj)
            if strcmp(get(src_obj(k),'Type'),'tiledlayout')
                tile_count=tile_count+1;
                src_tiles(tile_count)=src_obj(k);
            end
        end

        % Copy tiled figures directly into the outer layout
        if isscalar(src_tiles)
            delete(tile_axes(n));
            copy_obj=copyobj(src_obj,tile_obj);
            tile_idx=find(src_obj==src_tiles,1);
            copy_obj(tile_idx).Layout.Tile=tile_nums(n);
            for k=1:numel(copy_obj)
                if isprop(copy_obj(k),'Type')
                    obj_type=get(copy_obj(k),'Type');
                else
                    obj_type='';
                end
                if (k~=tile_idx)&&isprop(copy_obj(k),'Layout')&&...
                   (~strcmp(obj_type,'legend'))&&...
                   (~strcmp(obj_type,'colorbar'))
                    layout_obj=copy_obj(k).Layout;
                    if (~isempty(layout_obj))&&isprop(layout_obj,'Tile')
                        copy_obj(k).Layout.Tile=tile_nums(n);
                    end
                end
            end
            tile_found=true;
        else
            tile_found=false;
        end

        % Process figures without a source tiled layout
        if ~tile_found

            % Count top-level source axes
            axis_count=0;
            for k=1:numel(src_obj)
                if strcmp(get(src_obj(k),'Type'),'axes')
                    axis_count=axis_count+1;
                end
            end

            % Collect top-level source axes
            src_axes=gobjects(axis_count,1);
            axis_count=0;
            for k=1:numel(src_obj)
                if strcmp(get(src_obj(k),'Type'),'axes')
                    axis_count=axis_count+1;
                    src_axes(axis_count)=src_obj(k);
                end
            end

            % Copy non-axes chart objects directly
            if isempty(src_axes)
                delete(tile_axes(n));
                copy_obj=copyobj(src_obj,tile_obj);
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
                            copy_obj(k).Layout.Tile=tile_nums(n);
                            tile_found=true;
                        end
                    end
                end
            else

                % Measure source axes positions
                axis_pos=zeros(numel(src_axes),4);
                plot_pos=zeros(numel(src_axes),4);
                for k=1:numel(src_axes)
                    old_units=src_axes(k).Units;
                    src_axes(k).Units='normalized';
                    axis_pos(k,:)=src_axes(k).OuterPosition;
                    plot_pos(k,:)=src_axes(k).Position;
                    src_axes(k).Units=old_units;
                end

                % Detect overlaid source axes
                pos_tol=sqrt(eps);
                overlap_found=false;
                for k=1:numel(src_axes)
                    for m=k+1:numel(src_axes)
                        hor_ovlp=min(plot_pos(k,1)+plot_pos(k,3),...
                                     plot_pos(m,1)+plot_pos(m,3))-...
                                 max(plot_pos(k,1),plot_pos(m,1));
                        ver_ovlp=min(plot_pos(k,2)+plot_pos(k,4),...
                                     plot_pos(m,2)+plot_pos(m,4))-...
                                 max(plot_pos(k,2),plot_pos(m,2));
                        if (hor_ovlp>pos_tol)&&(ver_ovlp>pos_tol)
                            overlap_found=true;
                        end
                    end
                end

                % Preserve overlaid axes in a tile-sized panel
                if overlap_found
                    panel_obj=uipanel(fig_obj,'BorderType','none',...
                                      'Units','normalized',...
                                      'Position',tile_pos(n,:));
                    copy_obj=copyobj(src_obj,panel_obj);
                    panel_count=panel_count+1;
                    panel_axes(panel_count)=tile_axes(n);
                    panel_objs(panel_count)=panel_obj;
                    tile_found=true;
                else

                    % Locate source axes centres
                    axis_centre=[axis_pos(:,1)+axis_pos(:,3)/2 ...
                                 axis_pos(:,2)+axis_pos(:,4)/2];

                    % Seed row and column grids
                    col_grid=uniquetol(axis_centre(:,1));
                    row_grid=uniquetol(axis_centre(:,2));

                    % Refine columns from the densest source row
                    row_nums=zeros(numel(src_axes),1);
                    for k=1:numel(src_axes)
                        [~,row_nums(k)]=min(abs(row_grid-axis_centre(k,2)));
                    end
                    row_count=zeros(numel(row_grid),1);
                    for k=1:numel(row_grid)
                        row_count(k)=sum(row_nums==k);
                    end
                    [max_cols,best_row]=max(row_count);
                    if max_cols>1
                        col_grid=sort(axis_centre(row_nums==best_row,1));
                    end

                    % Refine rows from the densest source column
                    col_nums=zeros(numel(src_axes),1);
                    for k=1:numel(src_axes)
                        [~,col_nums(k)]=min(abs(col_grid-axis_centre(k,1)));
                    end
                    col_count=zeros(numel(col_grid),1);
                    for k=1:numel(col_grid)
                        col_count(k)=sum(col_nums==k);
                    end
                    [max_rows,best_col]=max(col_count);
                    if max_rows>1
                        row_grid=sort(axis_centre(col_nums==best_col,2));
                    end

                    % Create a nested layout for this source figure
                    delete(tile_axes(n));
                    nest_obj=tiledlayout(tile_obj,numel(row_grid),numel(col_grid),...
                                         'TileSpacing','loose','Padding','loose');
                    nest_obj.Layout.Tile=tile_nums(n);

                    % Copy source graphics into the nested layout
                    copy_obj=copyobj(src_obj,nest_obj);

                    % Put copied axes into matching nested tiles
                    axis_tile=zeros(numel(src_axes),1);
                    for k=1:numel(src_axes)
                        col_hits=find((col_grid>=axis_pos(k,1)-pos_tol)&...
                                      (col_grid<=axis_pos(k,1)+axis_pos(k,3)+pos_tol));
                        row_hits=find((row_grid>=axis_pos(k,2)-pos_tol)&...
                                      (row_grid<=axis_pos(k,2)+axis_pos(k,4)+pos_tol));
                        if isempty(col_hits)
                            [~,col_hits]=min(abs(col_grid-axis_centre(k,1)));
                        end
                        if isempty(row_hits)
                            [~,row_hits]=min(abs(row_grid-axis_centre(k,2)));
                        end
                        col_idx=min(col_hits);
                        row_idx=numel(row_grid)-max(row_hits)+1;
                        col_span=max(col_hits)-col_idx+1;
                        row_span=max(row_hits)-min(row_hits)+1;
                        axis_tile(k)=(row_idx-1)*numel(col_grid)+col_idx;
                        axis_idx=find(src_obj==src_axes(k),1);
                        if isempty(axis_idx)
                            error('source axes must be top-level figure children.');
                        end
                        copy_obj(axis_idx).Layout.Tile=axis_tile(k);
                        copy_obj(axis_idx).Layout.TileSpan=[row_span col_span];
                    end

                    % Confirm that source axes were placed
                    tile_found=true;

                end

            end

        end

        % Reject unsupported source figures
        if ~tile_found
            error('source figure contains no tile-capable graphics objects.');
        end

        % Preserve axes colour maps
        src_map_axes=findobj(src_obj,'Type','axes');
        copy_axes=findobj(copy_obj,'Type','axes');
        if numel(src_map_axes)==numel(copy_axes)
            for k=1:numel(copy_axes)
                copy_axes(k).Colormap=src_map_axes(k).Colormap;
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
if panel_count>0
    move_panels(fig_obj,[],panel_axes(1:panel_count),...
                panel_objs(1:panel_count));
    fig_obj.SizeChangedFcn={@move_panels,panel_axes(1:panel_count),...
                            panel_objs(1:panel_count)};
end

end

% Synchronises overlay panels with their tiled placeholders
function move_panels(~,~,panel_axes,panel_objs)
for n=1:numel(panel_objs)
    if isvalid(panel_axes(n))&&isvalid(panel_objs(n))
        panel_objs(n).Position=panel_axes(n).OuterPosition;
    end
end
end

% Consistency enforcement
function grumble(fig_files)
if (~iscell(fig_files))||isempty(fig_files)||(~ismatrix(fig_files))
    error('fig_files must be a non-empty two-dimensional cell array of file names.');
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

