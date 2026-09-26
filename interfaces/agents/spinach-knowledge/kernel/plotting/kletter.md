# kernel/plotting/kletter.m

- Signature: `kletter(letter_label)`

## Purpose

Draws an academic journal style letter label in the top left corner of the current axis set. The label is placed inside the plot box, with its left edge and its cap line 10 points away from the left and the top edge of the box, whatever the proportions of the figure or of the axes.

## Physical / mathematical content

- None: `kletter` is a figure annotation utility; it draws a panel letter and carries no physical or mathematical content.

## Numerical / algorithmic content

- The offsets are absolute (10 points) inside the rendered plot box when the function is called, so the figure should already have its final size, and a tiled layout all of its tiles, at that point; the position is then stored as a fraction of the plot box (normalised text units), so the label follows the axes if they are resized later; the text object carries the tag `kletter`, and `fig2tiles` deletes and re-applies tagged labels on the retiled axes of a merged figure so their offsets are exact there.
- The label is a `text` object that belongs to the current axes, so it is carried along when the axes are copied into a tiled figure by `fig2tiles`.

## Implementation structure

- Check consistency
- Offset from the box edges, points
- Rendered plot box of the current axes in points
- Label position as a fraction of the plot box
- Place the tagged label with its cap line at the top offset
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- The grumbler requires `letter_label` to be a one-element character string.
