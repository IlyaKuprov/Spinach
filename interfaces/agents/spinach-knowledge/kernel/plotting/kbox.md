# kernel/plotting/kbox.m

- Signature: `kbox() % #NGRUM`

## Purpose

Creates a tickless frame around the current axes using a line overlay in data coordinates. It draws a rectangle in 2D and the 12 edges of the axis-limit box in 3D.

## Physical / mathematical content

## Numerical / algorithmic content

- The overlay follows the axes limits, line width, and colour; listeners update it when relevant axes properties or children change. It is excluded from autoscaling and placed above the other axes children.

## Outputs

- Creates or updates one box overlay in the current axes; returns no output.
