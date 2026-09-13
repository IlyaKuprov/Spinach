# kernel/plotting/kletter.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/kletter.m`
- Signature: `kletter(letter_label)`
- Total lines: 86

## Purpose

Draws an academic journal style letter label in the top left corner of the current axis set. The label is placed inside the outer box of the axes (the whole figure when there is one axis set), with its left edge and its cap line 10 points away from the left and the top edge of that box, whatever the proportions of the figure.

## Physical / mathematical content

- None: `kletter` is a figure annotation utility; it draws a panel letter and carries no physical or mathematical content.

## Numerical / algorithmic content

- The offsets are absolute (10 points) when the function is called, so the figure should already have its final size at that point; the position is then stored as a fraction of the plot box (normalised text units), so the label follows the axes if they are resized or retiled later.
- The label is a `text` object that belongs to the current axes, so it is carried along when the axes are copied into a tiled figure by `fig2tiles`.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consistency; implemented by `grumble(letter_label)`.
- Lines 35-36: Offset from the box edges, points; `edge_offset=10`.
- Lines 38-41: Rendered plot box and outer box of the current axes in points; the axes units are switched to points, `tightPosition` (the rendered plot box, which is smaller than `Position` under a constrained aspect ratio such as `axis square`) and `OuterPosition` are read, and the units are restored.
- Lines 43-45: Label position as a fraction of the plot box; the left offset is added to the outer box left edge and the top offset is subtracted from the outer box top edge, both measured from the plot box corner and divided by the plot box width and height.
- Lines 47-50: Place the label with its cap line at the top offset; `text(...,'Units','normalized','HorizontalAlignment','left','VerticalAlignment','cap','FontWeight','bold','FontSize',16)`.

## Implementation structure

- Check consistency
- Offset from the box edges, points
- Rendered plot box and outer box of the current axes in points
- Label position as a fraction of the plot box
- Place the label with its cap line at the top offset
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gca()`, `tightPosition()`, `text()`.
- The grumbler requires `letter_label` to be a one-element character string.
