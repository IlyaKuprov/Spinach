# kernel/plotting/write_movie.m

- Source: [kernel/plotting/write_movie.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/write_movie.m)
- Wiki: [write_movie.m](https://spindynamics.org/wiki/index.php?title=write_movie.m)

## Purpose

Capture a rotating view of the current 3D plot into a video file. The documented callers include `volplot()`, `cst_display()`, and `hfc_display()`.

## Input and export

`file_name` must be a character array and is passed to `VideoWriter(file_name,'MPEG-4')`; the writer quality is set to 100. The function sets the current axes to `vis3d`, opens the writer, and repeats 359 times: capture `getframe(gcf)` and pass it to the active `writeVideo` call, pause for `0.05` seconds, then call `camorbit(1,0)`. It closes the writer after the loop. The export profile is MPEG-4 and the frames come from the current figure.

## Figure and file effects

This writes the specified video file and leaves the current axes in `vis3d` mode with the camera orbited; it does not restore the prior view. Source inspection confirms that `writeVideo(writerObj,getframe(gcf))` is executable, not commented out, so each iteration captures and writes a frame. The function returns no MATLAB output.
