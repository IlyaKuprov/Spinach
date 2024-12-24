% Converts an XML file into a Matlab structure. Syntax:
%
%                 xml=parsexml(filename)
%
% This function is called by x2spinach() during import
% of SpinXML files. Direct calls are discouraged.
%
% Parameters:
%
%   filename - a string with the XML file name
%
% Outputs:
%
%   xml      - Matlab structure containing the 
%              information from the XML file
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=parsexml.m>

function xml=parsexml(filename)

% Check consistency
grumble(filename);

% TODO: validation against the schema,
%       best just call Java here

% Read the file
try
   tree=xmlread(filename);
catch
   error('Could not read XML file %s.',filename);
end

% Parse child nodes
try
   xml=parseChildNodes(tree);
catch
   error('Could not parse XML file %s.',filename);
end

end

% Child node parsing
function children=parseChildNodes(theNode)
children=[];
if theNode.hasChildNodes
    childNodes=theNode.getChildNodes;
    numChildNodes=childNodes.getLength;
    allocCell=cell(1,numChildNodes);
    children=struct('name',allocCell,'attributes',allocCell,...
                    'data',allocCell,'children',allocCell);
    for n=1:numChildNodes
        theChild=childNodes.item(n-1);
        children(n)=makeStructFromNode(theChild);
    end
end
end

% Structure creation
function nodeStruct=makeStructFromNode(theNode)
nodeStruct=struct('name',char(theNode.getNodeName),...
                  'attributes',parseAttributes(theNode),...
                  'data','','children',parseChildNodes(theNode));
if ismember('getData',methods(theNode))
   nodeStruct.data=char(theNode.getData); 
else
   nodeStruct.data='';
end
end

% Attribute parsing
function attributes=parseAttributes(theNode)
attributes=[];
if theNode.hasAttributes
   theAttributes=theNode.getAttributes;
   numAttributes=theAttributes.getLength;
   allocCell=cell(1,numAttributes);
   attributes=struct('name',allocCell,'value',allocCell);
   for n=1:numAttributes
      attrib=theAttributes.item(n-1);
      attributes(n).name=char(attrib.getName);
      attributes(n).value=char(attrib.getValue);
   end
end
end

% Consistency enforcement
function grumble(filename)
if ~exist(filename,'file')
    error('the file specified as not found.');
end
end

% Malcolm Levitt's email to IK, 08 March 2013, responding 
% to IK's insistence that the description of MRI in terms
% of Bloch equations is not quantum mechanical:
%
% "In my mind this displays an almost incredible ignorance of possibly the most
% important single subfield of Mag Res. [...] It's literally inexcusable, on ev-
% ery possible level. I had put down your comment as one of those stupid things
% that everyone lets slip now and again but forgiveable as an indiscretion that
% is committed in the heat of the moment, or inexperience. But even after time
% for reflection, you try to maintain this view. I am truly shocked. [...] I am
% not going to mince my words here. As head of Magnetic Resonance here I want no
% one in my section who is so ignorant and badly informed about such an import-
% ant subfield of magnetic resonance, and who tries to propropagate such a comp-
% letely erroneous and misleading view to students and colleagues, despite hav-
% ing the error pointed out. In my judgement a scientist who can make such a ter-
% rible error of scientific judgement and then redouble the error by trying to
% defend it does not have a bright research career ahead of them, whatever their
% superficial promise, and wherever they are. All their talent will come to noth-
% ing since it will be unsupported by good sense and led astray by arrogance. 
% All the good scientists I know are, on a deep level, humble, and that humble-
% ness is essential for their achievements. That is my professional judgement.
% I still hope that you reflect and sincerely recognize your bad mistake and ack-
% nowledge it, so that something like normal service can be resumed. If that hap-
% pens, then this is the end of this affair and it goes not further."
%
% (IK still thinks that Bloch equations belong to classical physics.)

