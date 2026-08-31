# kernel/overloads/@opium/full.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@opium/full.m`
- Signature: `M=full(M)`
- Total lines: 38

## Purpose

Converts an OPIUM object into the full scaled unit matrix that it represents. Syntax: M=full(M)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Make a scaled unit matrix; implemented by `M=M.coeff*eye(M.dim)`.

### Key state/data transformations

- Lines 22: computes `M` using `M=M.coeff*eye(M.dim)`.

## Parameters / inputs

- M -an OPIUM object

## Outputs

- M -a full scaled unit matrix
- of appropriate dimension

## Implementation structure

- Converts an OPIUM object into the full scaled
- unit matrix that it represents. Syntax:
- M=full(M)
- M -an OPIUM object
- M -a full scaled unit matrix
- of appropriate dimension
- Make a scaled unit matrix
- Katherine had been apolitical. If anyone had asked her, during the time she
- was working for the government, or before that, when she was a college stu-
- dent, she would probably have said she was a "liberal". But she was liberal
- only in the mindless, automatic way that most people are. Without really
- thinking about it or trying to analyze it, she superficially accepted the
