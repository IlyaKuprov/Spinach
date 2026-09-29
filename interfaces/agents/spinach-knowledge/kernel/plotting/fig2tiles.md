# kernel/plotting/fig2tiles.m

- Signature: `[fig_obj,tile_obj]=fig2tiles(fig_files,fig_size)`

Combines saved MATLAB figure files into a new tiled figure. It returns the new figure handle and its outer tiled-layout handle.

## Inputs

- `fig_files`: nonempty 2-D cell array of existing `.fig` file-name character vectors. Its row and column dimensions determine the outer tile grid; entries are placed in that same row-major order.
- `fig_size`: two-element row vector of positive whole numbers, interpreted as the merged figure width and height in screen pixels.

## Layout and graphical content

The helper creates the output figure at pixel position `[0 0 fig_size]`, then makes a loose-spacing, loose-padding tiled layout with one outer tile per input file. It measures tile geometry after the layout is created and rendered. Each source figure is opened invisibly while its graphics are copied and then closed. A source tiled layout is copied into its outer tile; for figures without one, the source axes positions are used to arrange their contents in a nested layout. Matching axes colormaps are carried over when available.

The routine also retile-copies associated graphics such as legends and overlays. Panel objects are aligned to their axes' outer positions and follow later figure resizes. Text objects tagged `kletter` are removed and redrawn on their retiled axes so their offsets are recalculated for the merged layout.

The new figure is shown using MATLAB's default root visibility when its outer dimensions fit the screen. If it is larger than the screen it remains invisible and a warning points to `savefig` and `exportgraphics` for retaining or exporting it. Because tile positions are measured when the layout is created, choose the final `fig_size` at construction rather than resizing the figure later.

## Existing syntax

`[fig_obj,tile_obj]=fig2tiles(fig_files,fig_size)`

## References

- [Source: `kernel/plotting/fig2tiles.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/fig2tiles.m)
- [Spinach Wiki: `fig2tiles.m`](https://spindynamics.org/wiki/index.php?title=fig2tiles.m)
