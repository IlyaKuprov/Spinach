# kernel/overloads/save_anyway.m

- Signature: `save_anyway(file_name,variable)`

## Purpose

A wrapper intended to trick SPMD blocks into saving data. Can only save one variable at a time, its name in the mat file is "variable". Syntax: save_anyway(file_name,variable)

## Physical / mathematical content

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- file_name -a character string specifying the
- file name
- variable -the variable to be saved

## Implementation structure

- A wrapper intended to trick SPMD blocks into saving data. Can
- only save one variable at a time, its name in the mat file is
- "variable". Syntax:
- save_anyway(file_name,variable)
- file_name -a character string specifying the
- file name
- variable -the variable to be saved
- Check consistency
- Just call save
- Consistncy enforcement
- Life struggles to survive here, and while some clings
- to a tenacious existence, it is anemic and sickly.
