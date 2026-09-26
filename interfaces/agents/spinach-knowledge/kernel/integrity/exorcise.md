# kernel/integrity/exorcise.m

- Signature: `exorcise(mode)`

## Purpose

Searches Spinach distribution folders for any functions that do not conform to the house style. Opens the first one and complains to the console. Syntax: exorcise(mode)

## Physical / mathematical content

- Integrity-control utilities. These files check distribution state, path collisions, style conformance, sniffer databases, and other safeguards that protect Spinach reproducibility.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- mode -'online' checks the documentation
- Wiki for the corresponding page; 'offline'
- skips the Wiki check
- Any user contribution that this function has something to
- say about will either be brought under the house style, or
- rejected back to the user, depending on the amount of work
- involved. Always run this function before a commit if you
- have write access to Spinach repository.

## Implementation structure

- Searches Spinach distribution folders for any functions that
- do not conform to the house style. Opens the first one and
- complains to the console. Syntax:
- exorcise(mode)
- mode -'online' checks the documentation
- Wiki for the corresponding page; 'offline'
- skips the Wiki check
- Any user contribution that this function has something to
- say about will either be brought under the house style, or
- rejected back to the user, depending on the amount of work
- involved. Always run this function before a commit if you
- have write access to Spinach repository.
