# kernel/plotting/kvertshade.m

- MATLAB implementation: [kernel/plotting/kvertshade.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/kvertshade.m)

## Purpose and use

Call kvertshade(a,b) to shade a vertical band between x=a and x=b in the current axes (gca). It is the shaded-strip counterpart of Matlab's xline and the companion of khorzshade. It is a plotting helper, not a numerical or physical operation, and returns no value.

## Behaviour

The band is a Matlab xregion object with FaceColor [0 0 0], FaceAlpha 0.1, and no edge line; on a white background it renders as 10% grey, and overlapping bands darken multiplicatively (two layers give 19% grey). Being anchored to the axes, it spans their full height, follows later changes of the axis limits, and is drawn on the bottom layer, underneath lines and markers. Existing plot content is kept without hold on.

The arguments a and b are in x-axis data units and may be given in either order; each must be a finite real numeric scalar, otherwise the function stops with an error.

Source documentation: [kvertshade.m](https://spindynamics.org/wiki/index.php?title=kvertshade.m).
