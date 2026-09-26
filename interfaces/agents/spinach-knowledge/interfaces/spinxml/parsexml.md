# interfaces/spinxml/parsexml.m

- Signature: `xml=parsexml(filename)`

## Purpose

Converts an XML file into a Matlab structure. Syntax: xml=parsexml(filename) This function is called by x2spinach() during import of SpinXML files. Direct calls are discouraged.

## Physical / mathematical content

## Numerical / algorithmic content

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
