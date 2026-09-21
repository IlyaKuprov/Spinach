# kernel/overloads/@polyadic/gpuArray.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/gpuArray.m`
- Signature: `p=gpuArray(p)`
- Total lines: 59

## Purpose

Uploads all components of a polyadic object to the GPU. The object still looks like a polyadic to Matlab, but all of its constituent matrices become gpuArrays. Syntax: p=gpuArray(p)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- p -polyadic object

## Outputs

- p -polyadic object with all cores,
- prefixes and suffixes uploaded
- to the current GPU
- Note: GPUs are not good at permuting array dimensions. If you find
- yourself using polyadics, do check that the CPU isn't faster.

## Implementation structure

- Uploads all components of a polyadic object to the GPU. The object
- still looks like a polyadic to Matlab, but all of its constituent
- matrices become gpuArrays. Syntax:
- p=gpuArray(p)
- p -polyadic object
- p -polyadic object with all cores,
- prefixes and suffixes uploaded
- to the current GPU
- Note: GPUs are not good at permuting array dimensions. If you find
- yourself using polyadics, do check that the CPU isn't faster.
- Check consistency
- Upload cores

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`.
