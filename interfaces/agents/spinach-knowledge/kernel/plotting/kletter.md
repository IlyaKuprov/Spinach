# kernel/plotting/kletter.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/kletter.m`
- Signature: `kletter(letter_label)`
- Total lines: 87

## Purpose

Draws an academic journal style letter label in the top left corner of the current axis set. The label is placed inside the plot box, with its left edge and its cap line 10 points away from the left and the top edge of the box, whatever the proportions of the figure or of the axes.

## Physical / mathematical content

- None: `kletter` is a figure annotation utility; it draws a panel letter and carries no physical or mathematical content.

## Numerical / algorithmic content

- The offsets are absolute (10 points) inside the rendered plot box when the function is called, so the figure should already have its final size, and a tiled layout all of its tiles, at that point; the position is then stored as a fraction of the plot box (normalised text units), so the label follows the axes if they are resized later; the text object carries the tag `kletter`, and `fig2tiles` deletes and re-applies tagged labels on the retiled axes of a merged figure so their offsets are exact there.
- The label is a `text` object that belongs to the current axes, so it is carried along when the axes are copied into a tiled figure by `fig2tiles`.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(letter_label)`.
- Lines 37-38: Offset from the box edges, points; `edge_offset=10`.
- Lines 40-42: Rendered plot box of the current axes in points; the axes units are switched to points, `tightPosition` (the rendered plot box, which is smaller than `Position` under a constrained aspect ratio such as `axis square`) is read, and the units are restored.
- Lines 44-46: Label position as a fraction of the plot box; the left offset divided by the box width, and one minus the top offset divided by the box height.
- Lines 48-51: Place the tagged label with its cap line at the top offset; `text(...,'Units','normalized','HorizontalAlignment','left','VerticalAlignment','cap','FontWeight','bold','FontSize',16,'Tag','kletter')`.

## Implementation structure

- Check consistency
- Offset from the box edges, points
- Rendered plot box of the current axes in points
- Label position as a fraction of the plot box
- Place the tagged label with its cap line at the top offset
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gca()`, `tightPosition()`, `text()`.
- The grumbler requires `letter_label` to be a one-element character string.
