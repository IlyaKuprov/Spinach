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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Check consistency; implemented by `grumble(filename)`.
- Lines 26-29: Use Matlab's native XML engine, avoiding Java; implemented by `try`.
- Lines 28-29: Read the file; implemented by `try`.
- Lines 35-36: Parse child nodes; implemented by `try`.

### Key state/data transformations

- Lines 30: computes `tree` using `tree=xmlread(filename,'XMLEngine','maxp')`.
- Lines 37: computes `xml` using `xml=parseChildNodes(tree)`.

### Local helper functions

- Line 45: `parseChildNodes()` — `function children=parseChildNodes(theNode)`.
  - Representative operation: `children=[]`.
  - Representative operation: `if isprop(theNode,'Children')&&(~isempty(theNode.Children))`.
- Line 61: `makeStructFromNode()` — `function nodeStruct=makeStructFromNode(theNode)`.
  - Representative operation: `if isprop(theNode,'TagName')`.
  - Representative operation: `nodeName=char(theNode.TagName)`.
- Line 87: `parseAttributes()` — `function attributes=parseAttributes(theNode)`.
  - Representative operation: `attributes=[]`.
  - Representative operation: `if isprop(theNode,'HasAttributes')&&theNode.HasAttributes`.
- Line 103: `grumble()` — `function grumble(filename)`. Malcolm Levitt's email to IK, 08 March 2013, responding to IK's insistence that the description of MRI in terms
  - Representative operation: `if (~ischar(filename))||isempty(filename)`.
  - Representative operation: `error('filename must be a non-empty character string.')`.

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
