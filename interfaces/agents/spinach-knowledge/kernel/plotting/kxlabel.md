# kernel/plotting/kxlabel.m

- MATLAB implementation: [kernel/plotting/kxlabel.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/kxlabel.m)

Set an axes' x-axis label with LaTeX rendering and apply the Spinach tick-label style.

## Call

`kxlabel(varargin)`

Supply the same arguments accepted by MATLAB's `xlabel` function (for example, label text and its property-value options). The wrapper forwards those arguments to `xlabel` and appends `'Interpreter','latex'`; the label is therefore requested in LaTeX, rather than the default interpreter. There is no separate unit convention imposed by this plotting helper.

## Effect and scope

After setting the x label, the function sets `TickLabelInterpreter` to `'latex'` and `FontSize` to `12` on `gca` (the current axes). These values concern the axes' tick labels, not a numerical result or data transformation. No output is returned.

The label call receives the caller's arguments, but the follow-up style call explicitly uses `gca`; if an axes target is supplied to `xlabel`, the tick styling still targets the current axes. There is no physical or numerical calculation here.

## Source link

[`kxlabel.m` on the Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=kxlabel.m)
