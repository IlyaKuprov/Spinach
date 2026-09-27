# interfaces/spinxml/parsexml.m

- Signature: `xml=parsexml(filename)`

## Purpose

Reads an XML file and converts its document nodes to a MATLAB structure. x2spinach() uses it when importing SpinXML files; direct calls are discouraged.

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- filename — a non-empty character string naming an existing XML file

## Outputs

- xml — a structure with node names, attributes, data, and child nodes

## Implementation structure

- Checks that filename is a non-empty character string and names an existing file.
- Reads the XML with MATLAB's xmlread using the XMLEngine setting `maxp`.
- Recursively converts child nodes to structures with `name`, `attributes`, `data`, and `children` fields; attributes are represented by name/value pairs.
- Uses element tag names, and labels text and comment nodes `#text` and `#comment`, respectively. Parse or read failures raise an error.
