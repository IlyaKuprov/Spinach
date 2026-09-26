# kernel/plotting/fig2tiles.m

- Signature: `[fig_obj,tile_obj]=fig2tiles(fig_files,fig_size)`

## Purpose

Combines Matlab figure files into a single tiled figure. Syntax: [fig_obj,tile_obj]=fig2tiles(fig_files,fig_size)

## Physical / mathematical content

## Numerical / algorithmic content

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

- After the layout geometry is updated, text objects tagged `kletter` (panel letters drawn by `kletter.m` in the source figures) are deleted and re-applied on their retiled axes, so their offsets inside the plot box are exact in the merged figure.

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
