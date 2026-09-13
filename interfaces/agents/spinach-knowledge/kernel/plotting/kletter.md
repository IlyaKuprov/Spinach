# kernel/plotting/kletter.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/kletter.m`
- Signature: `kletter(letter_label)`
- Total lines: 83

## Purpose

Draws an academic journal style letter label in the top left corner of the current axis set. The label is placed inside the outer box of the axes (the whole figure when there is one axis set), with its left edge and its cap line 10 points away from the left and the top edge of that box, whatever the proportions of the figure.

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The offsets are absolute (points), so they do not scale with the axes width or height; they are measured when the function is called, and the figure must have its final size at that point.
- The label is a `text` object that belongs to the current axes, so it is carried along when the axes are copied into a tiled figure by `fig2tiles`.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(letter_label)`.
- Lines 32-33: Offset from the box edges, points; `edge_offset=10`.
- Lines 35-38: Plot box and outer box of the current axes in points; the axes units are switched to points, `Position` and `OuterPosition` are read, and the units are restored.
- Lines 40-42: Label position relative to the plot box corner; the left offset is added to the outer box left edge and the top offset is subtracted from the outer box top edge, both expressed from the plot box corner because text positions are measured from there.
- Lines 44-47: Place the label with its cap line at the top offset; `text(...,'Units','points','HorizontalAlignment','left','VerticalAlignment','cap','FontWeight','bold','FontSize',16)`.

## Implementation structure

- Check consistency
- Offset from the box edges, points
- Plot box and outer box of the current axes in points
- Label position relative to the plot box corner
- Place the label with its cap line at the top offset
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gca()`, `text()`.
- The grumbler requires `letter_label` to be a one-element character string.
