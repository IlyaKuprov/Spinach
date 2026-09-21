# interfaces/spinxml/parsexml.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/spinxml/parsexml.m`
- Signature: `xml=parsexml(filename)`
- Total lines: 138

## Purpose

Converts an XML file into a Matlab structure. Syntax: xml=parsexml(filename) This function is called by x2spinach() during import of SpinXML files. Direct calls are discouraged.

## Physical / mathematical content

- This file belongs to the `interfaces` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `parseChildNodes()`, `makeStructFromNode()`, `parseAttributes()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- filename -a string with the XML file name

## Outputs

- xml -Matlab structure containing the
- information from the XML file

## Implementation structure

- Converts an XML file into a Matlab structure. Syntax:
- xml=parsexml(filename)
- This function is called by x2spinach() during import
- of SpinXML files. Direct calls are discouraged.
- filename -a string with the XML file name
- xml -Matlab structure containing the
- information from the XML file
- Check consistency
- Use Matlab's native XML engine, avoiding Java
- Read the file
- Parse child nodes
- Child node parsing

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `xmlread()`, `parseChildNodes()`, `isprop()`, `childNodes()`, `children()`, `makeStructFromNode()`, `char()`, `parseAttributes()`, `attributes()`, `ischar()`, `exist()`.
