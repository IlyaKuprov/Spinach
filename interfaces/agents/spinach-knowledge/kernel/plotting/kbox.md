# kernel/plotting/kbox.m

Source: [kernel/plotting/kbox.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/kbox.m) · [Wiki](https://spindynamics.org/wiki/index.php?title=kbox.m)

- Signature: `kbox()`

## Behaviour and coordinates

`kbox` operates on `gca` and draws a tickless outline with a tagged line object in the axes' data coordinates. Its bounds are the current `XLim`, `YLim`, and (in 3-D) `ZLim`; it does not choose numeric axis limits. In the standard 2-D view, where `abs(View*[0;1]-90) <= sqrt(eps)` and no child has nonempty `ZData`, it connects the four corners in the `x-y` plane at `z=0`. Otherwise it connects the 12 edges of the current limit box. The coordinate data are row vectors; each edge is represented by its two endpoints followed by `NaN` to break the line between edges. The 2-D line has empty `ZData`.

This is a rendering helper: the coordinates describe the axes frame, not a propagated signal or physical quantity.

## Rendering and updates

The function turns the axes `Box` off and `Layer` to `top`. The overlay uses the axes `XColor` and `LineWidth`, is unclipped, hidden from handle discovery and picking, and is excluded from the axis-limit calculations where those line properties are available. It is moved to the front of the axes' child order. Axes appdata holds its state and listeners; updates follow changes to observable limit, style, and child properties and the axes' `MarkedClean` event. A re-entry guard and a cached signature avoid recursive or unchanged updates. Calling `kbox` again removes the previous overlay and its listeners before creating a replacement; it also removes tagged orphan overlay axes left by an earlier implementation.

There are no arguments or returned values. The update helper returns early if the axes is no longer valid or an update is already in progress; there is no separate validation of axis bounds.
