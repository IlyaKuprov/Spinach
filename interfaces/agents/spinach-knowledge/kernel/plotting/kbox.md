# kernel/plotting/kbox.m

- Signature: `kbox() % #NGRUM`

## Purpose

Creates a tickless boxed frame around the current axes using ordinary line objects that live in the same data space as the plot. This is needed for spectrograms be- cause Matlab has dumb plotting defaults. Syntax: kbox()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- creates or updates a tickless axis box
- in the current axes

## Implementation structure

- Creates a tickless boxed frame around the current axes
- using ordinary line objects that live in the same data
- space as the plot. This is needed for spectrograms be-
- cause Matlab has dumb plotting defaults. Syntax:
- kbox()
- creates or updates a tickless axis box
- in the current axes
- Get the current axis object
- Remove any previous Spinach box overlay
- Delete orphaned overlay axes from the rejected implementation
- Create the box line in the plot axes
- Exclude the box line from autoscaling
