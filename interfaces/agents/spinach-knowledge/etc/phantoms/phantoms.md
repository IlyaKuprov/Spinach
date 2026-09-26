# etc/phantoms/phantoms.m

- Signature: `[R1Ph,R2Ph,PDPh,dims,npts]=phantoms(ph_name)`

## Purpose

MRI phantom library. Syntax: [R1Ph,R2Ph,PDPh,dims,npts]=phantoms(ph_name)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- ph_name -character string giving the name of the
- phantom (see the function text)

## Outputs

- R1Ph -a cube of R1 values
- R2Ph -a cube of R2 values
- PDPh -a cube of PD values
- dims -row vector of three cube dimensions, m
- npts -row vector of three cube dimensions, points

## Implementation structure

- MRI phantom library. Syntax:
- [R1Ph,R2Ph,PDPh,dims,npts]=phantoms(ph_name)
- ph_name -character string giving the name of the
- phantom (see the function text)
- R1Ph -a cube of R1 values
- R2Ph -a cube of R2 values
- PDPh -a cube of PD values
- dims -row vector of three cube dimensions, m
- npts -row vector of three cube dimensions, points
- Check consistency
- Get own location
- Load the breast phantom
