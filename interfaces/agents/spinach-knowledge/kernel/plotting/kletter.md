# kernel/plotting/kletter.m

Source: [kernel/plotting/kletter.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/kletter.m) · [Wiki](https://spindynamics.org/wiki/index.php?title=kletter.m)

- Signature: `kletter(letter_label)`

## Position and units

The function labels the current axes returned by `gca`. It temporarily changes the axes `Units` to `points`, obtains `tightPosition(ax_obj)`, and restores the previous units. For the measured plot-box width and height, `w` and `h`, in points, it computes the normalised text position as `x = 10/w` and `y = 1 - 10/h`. Thus the left edge and the text cap line are offset 10 points from the plot box's left and top edges. The position passed to `text` uses normalised units.

The created text is left-aligned, bold, 16-point, and tagged `kletter`; vertical alignment is `cap`. Its position is calculated when called, not maintained by a resize listener. The source notes that `fig2tiles.m` reapplies the tagged label when retiled axes are made for a merged figure; the knowledge page is [fig2tiles.md](fig2tiles.md).

## Input and output

`letter_label` must be a one-element character array; other inputs raise an error. The function creates a text object in the current axes and returns no output. It has no plot-data or physical-coordinate calculation, and no explicit guard for the measured plot-box dimensions.
