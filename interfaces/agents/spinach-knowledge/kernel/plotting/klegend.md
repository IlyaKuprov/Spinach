# kernel/plotting/klegend.m

- MATLAB implementation: [kernel/plotting/klegend.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/klegend.m)

## Purpose and use

leg_obj=klegend(varargin) is a styled wrapper around MATLAB legend. Pass the normal legend arguments as varargin; the wrapper forwards them to legend and returns the resulting legend object as leg_obj. Input conventions and units are those of MATLAB's legend function.

## Styling applied by the wrapper

The call appends the Interpreter='latex' and IconColumnWidth=10 settings to the forwarded arguments. It then sets the legend box face to ColorType='truecoloralpha' with ColorData=uint8([200 200 200 64]'), giving the box a translucent grey fill. The function does not otherwise define legend contents or labels; supply those through the standard legend arguments.

Source documentation: [klegend.m](https://spindynamics.org/wiki/index.php?title=klegend.m).
