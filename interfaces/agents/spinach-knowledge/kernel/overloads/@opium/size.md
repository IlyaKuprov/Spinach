# kernel/overloads/@opium/size.m

- Signature: `varargout=size(op,dim)`

## Purpose

The size of the matrix represented by the OPIUM. Syntax: answer=size(op,dim)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- op -an opium object
- dim -optional, dimension whose
- size is required

## Outputs

- answer -a vector with one or two elements

## Implementation structure

- Check that the requested dimension, if given, is 1 or 2
- Return the OPIUM object's dimensions in the requested output form
