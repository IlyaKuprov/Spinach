# kernel/plotting/kfigure.m

- Signature: `handle=kfigure(varargin)`

## Purpose

Sets four MATLAB root figure defaults to pre-R2025a values, creates a figure with the supplied arguments, and returns its handle.

## Physical / mathematical content

## Numerical / algorithmic content

- Sets `DefaultFigurePosition` to `[680 458 560 420]`, `DefaultFigureWindowStyle` to `normal`, `DefaultFigureMenuBar` to `figure`, and `DefaultFigureToolbar` to `figure` on `groot` before calling `figure`.

## Parameters / inputs

- `varargin` - arguments forwarded to MATLAB's `figure` function

## Outputs

- `handle` - handle returned by `figure(varargin{:})`
