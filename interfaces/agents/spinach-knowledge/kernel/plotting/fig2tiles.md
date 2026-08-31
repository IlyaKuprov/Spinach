# kernel/plotting/fig2tiles.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/fig2tiles.m`
- Signature: `[fig_obj,tile_obj]=fig2tiles(fig_files,fig_size)`
- Total lines: 396

## Purpose

Combines Matlab figure files into a single tiled figure. Syntax: [fig_obj,tile_obj]=fig2tiles(fig_files,fig_size)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `move_panels()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 36-37: Check consistency; implemented by `grumble(fig_files,fig_size)`.
- Lines 39-40: Get the tile grid size; implemented by `tile_rows=size(fig_files,1)`.
- Lines 43-44: Map cell array entries to tiled layout positions; implemented by `tile_nums=zeros(numel(fig_files),1)`.
- Lines 50-51: Note the figure visibility the caller has set as the default; implemented by `def_visible=get(groot,'defaultFigureVisible')`.
- Lines 53-55: Create the new figure off the screen at the size requested; implemented by `fig_obj=kfigure('Units','pixels','Position',[0 0 fig_size], 'Visible','off')`.
- Lines 57-59: Create a loose tiled layout; implemented by `tile_obj=tiledlayout(fig_obj,tile_rows,tile_cols,'TileSpacing','loose', 'Padding','loose')`.
- Lines 61-62: Allocate all outer tiles; implemented by `tile_axes=gobjects(numel(fig_files),1)`.
- Lines 69-70: Record outer tile positions; implemented by `drawnow`.
- Lines 75-76: Initialise overlay panel tracking; implemented by `panel_axes=gobjects(numel(fig_files),1)`.
- Lines 80-81: Import saved figures; implemented by `for n=1:numel(fig_files)`.
- Lines 83-84: Open source figure without displaying it; implemented by `src_fig=openfig(fig_files{n},'new','invisible')`.
- Lines 88-89: Collect source figure children; implemented by `src_obj=get(src_fig,'Children')`.
- Lines 91-92: Reject empty figures; implemented by `if isempty(src_obj)`.
- Lines 96-97: Count source tiled layouts; implemented by `tile_count=0`.
- Lines 104-105: Collect source tiled layouts; implemented by `src_tiles=gobjects(tile_count,1)`.
- Lines 114-115: Copy tiled figures directly into the outer layout; implemented by `if isscalar(src_tiles)`.
- Lines 140-141: Process figures without a source tiled layout; implemented by `if ~tile_found`.
- Lines 143-144: Count top-level source axes; implemented by `axis_count=0`.

### Control flow inferred from the code

- Line 45: `for` loop over `n=1:numel(fig_files)`.
- Line 64: `for` loop over `n=1:numel(fig_files)`.
- Line 71: `for` loop over `n=1:numel(fig_files)`.
- Line 81: `for` loop over `n=1:numel(fig_files)`.
- Line 92: conditional branch on `isempty(src_obj)`.
- Line 98: `for` loop over `k=1:numel(src_obj)`.
- Line 99: conditional branch on `strcmp(get(src_obj(k),'Type'),'tiledlayout')`.
- Line 107: `for` loop over `k=1:numel(src_obj)`.
- Line 108: conditional branch on `strcmp(get(src_obj(k),'Type'),'tiledlayout')`.
- Line 115: conditional branch on `isscalar(src_tiles)`.
- Line 120: `for` loop over `k=1:numel(copy_obj)`.
- Line 121: conditional branch on `isprop(copy_obj(k),'Type')`.
- Line 126: conditional branch on `(k~=tile_idx)&&isprop(copy_obj(k),'Layout')&&`.
- Line 130: conditional branch on `(~isempty(layout_obj))&&isprop(layout_obj,'Tile')`.

### Key state/data transformations

- Lines 40: computes `tile_rows` using `tile_rows=size(fig_files,1)`.
- Lines 41: computes `tile_cols` using `tile_cols=size(fig_files,2)`.
- Lines 44: computes `tile_nums` using `tile_nums=zeros(numel(fig_files),1)`.
- Lines 46: computes `[row_num,col_num]` using `[row_num,col_num]=ind2sub(size(fig_files),n)`.
- Lines 47: computes `tile_nums(n)` using `tile_nums(n)=(row_num-1)*tile_cols+col_num`.
- Lines 51: computes `def_visible` using `def_visible=get(groot,'defaultFigureVisible')`.
- Lines 54-55: computes `fig_obj` using `fig_obj=kfigure('Units','pixels','Position',[0 0 fig_size], 'Visible','off')`.
- Lines 58-59: computes `tile_obj` using `tile_obj=tiledlayout(fig_obj,tile_rows,tile_cols,'TileSpacing','loose', 'Padding','loose')`.
- Lines 62: computes `tile_axes` using `tile_axes=gobjects(numel(fig_files),1)`.
- Lines 63: computes `tile_pos` using `tile_pos=zeros(numel(fig_files),4)`.
- Lines 65: computes `tile_axes(n)` using `tile_axes(n)=nexttile(tile_obj,tile_nums(n))`.
- Lines 66: computes `tile_axes(n).Visible` using `tile_axes(n).Visible='off'`.
- Lines 72: computes `tile_pos(n,:)` using `tile_pos(n,:)=tile_axes(n).OuterPosition`.
- Lines 76: computes `panel_axes` using `panel_axes=gobjects(numel(fig_files),1)`.
- Lines 77: computes `panel_objs` using `panel_objs=gobjects(numel(fig_files),1)`.
- Lines 78: computes `panel_count` using `panel_count=0`.
- Lines 84: computes `src_fig` using `src_fig=openfig(fig_files{n},'new','invisible')`.
- Lines 89: computes `src_obj` using `src_obj=get(src_fig,'Children')`.

### Local helper functions

- Line 362: `move_panels()` — `function move_panels(~,~,panel_axes,panel_objs)`. Consistency enforcement
  - Representative operation: `for n=1:numel(panel_objs)`.
  - Representative operation: `if isvalid(panel_axes(n))&&isvalid(panel_objs(n))`.
- Line 371: `grumble()` — `function grumble(fig_files,fig_size)`.
  - Representative operation: `if (~isnumeric(fig_size))||(~isreal(fig_size))||(~isrow(fig_size))|| (numel(fig_size)~=2)||any(fig_size<1)||any(mod(fig_size,1)~=0)`.
  - Representative operation: `(numel(fig_size)~=2)||any(fig_size<1)||any(mod(fig_size,1)~=0)`.

## Parameters / inputs

- fig_files -cell array of character strings containing
- Matlab *.fig file names
- fig_size -width and height of the merged figure in scr-
- een pixels, a two-element row vector

## Outputs

- fig_obj -handle of the new Matlab figure
- tile_obj -handle of the tiled layout object
- Note: tile geometry is measured at the moment the layout is cre-
- ated, and so the figure must already have its final size at
- that point. Matlab shrinks a visible figure to fit the disp-
- lay; the merge therefore runs off the screen, and the figure
- is only shown at the end, at the visibility the caller has
- set as the figure default, when its outer extent fits.
- Figures bigger than the screen stay invisible and must be
- written out with exportgraphics.m or print.m; if they are
- reopened later with openfig.m, Matlab refits them to the
- screen and the size requested here is lost.

## Implementation structure

- Combines Matlab figure files into a single tiled figure. Syntax:
- [fig_obj,tile_obj]=fig2tiles(fig_files,fig_size)
- fig_files -cell array of character strings containing
- Matlab *.fig file names
- fig_size -width and height of the merged figure in scr-
- een pixels, a two-element row vector
- fig_obj -handle of the new Matlab figure
- tile_obj -handle of the tiled layout object
- Note: tile geometry is measured at the moment the layout is cre-
- ated, and so the figure must already have its final size at
- that point. Matlab shrinks a visible figure to fit the disp-
- lay; the merge therefore runs off the screen, and the figure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ind2sub()`, `tile_nums()`, `get()`, `kfigure()`, `tiledlayout()`, `gobjects()`, `tile_axes()`, `nexttile()`, `tile_pos()`, `openfig()`, `strcmp()`, `src_obj()`, `src_tiles()`, `isscalar()`, `delete()`.
