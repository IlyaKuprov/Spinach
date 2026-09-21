# kernel/utilities/dfpt.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/dfpt.m`
- Signature: `subgraphs=dfpt(conmatrix,max_sg_size)`
- Total lines: 110

## Purpose

Graph partitioning module. Analyzes the system connectivity graph and creates a list of all connected subgraphs of up to the user-specified size by crawling the graph in all available directions. Syntax: subgraphs=dfpt(conmatrix,max_sg_size)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- conmatrix -[nspins x nspins] matrix with 1 for
- connected spins and 0 elsewhere
- max_sg_size -maximum connected subgraph size

## Outputs

- subgraphs -[n_subgraphs x nspins] matrix; each
- row contains 1 for spins that belong
- to the subgraph and 0 for spins that
- do not

## Implementation structure

- Graph partitioning module. Analyzes the system connectivity graph and
- creates a list of all connected subgraphs of up to the user-specified
- size by crawling the graph in all available directions. Syntax:
- subgraphs=dfpt(conmatrix,max_sg_size)
- conmatrix -[nspins x nspins] matrix with 1 for
- connected spins and 0 elsewhere
- max_sg_size -maximum connected subgraph size
- subgraphs -[n_subgraphs x nspins] matrix; each
- row contains 1 for spins that belong
- to the subgraph and 0 for spins that
- do not
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `uint32()`, `subgraphs()`, `any()`, `conmatrix()`, `neighbours()`, `false()`, `subgraph()`, `cell2mat()`, `row_index()`, `logical()`, `islogical()`, `isscalar()`.
