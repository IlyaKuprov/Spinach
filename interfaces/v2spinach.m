% Imports NMR fid data in Varian format. Syntax:
%
%               vdata=v2spinach(inpath)
%
% Parameters:
%
%  inpath - character strong with the path to the 
%           Varian data folder
%
% Outputs:
%
%  vdata.fid      - complex free induction decay
%
% Modified for Spinach from the function by:
%
%   Dr. Mathias Nilsson
%   School of Chemistry, University of Manchester,
%   Oxford Road, Manchester M13 9PL, UK
%   Telephone: +44 (0) 161 306 4465
%   Fax: +44 (0) 161 275 4598
%   mathias.nilsson@manchester.ac.uk
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=v2spinach.m>

function vdata=v2spinach(inpath)

% Get the fid file
if nargin==0
    inpath=uigetdir(pwd,'Varian data directory');
end

% Check consistency
grumble(inpath);

% Open the file
fileid_fid=fopen([inpath filesep 'fid'],'r','b');

% Store the directory name
vdata.dirname=inpath;

% Read the header
vdata.nblocks=fread(fileid_fid,1,'int32');
vdata.ntraces=fread(fileid_fid,1,'int32');
vdata.np=fread(fileid_fid,1,'int32');
vdata.ebytes=fread(fileid_fid,1,'int32'); 
vdata.tbytes=fread(fileid_fid,1,'int32'); 
vdata.bbytes=fread(fileid_fid,1,'int32'); 
vdata.version_id=fread(fileid_fid,1,'int16'); 
vdata.status=fread(fileid_fid,1,'int16'); 
vdata.nbheaders=fread(fileid_fid,1,'int32');

% Preallocate fid array
vdata.fid=zeros(vdata.np/2,vdata.nblocks,'like',1i);

% Read in fid blocks
for m=1:vdata.nblocks
    
    % Read in block header
    vdata.block(m).scale=fread(fileid_fid,1,'int16');
    vdata.block(m).status=fread(fileid_fid,1,'int16');
    vdata.block(m).bitstatus=bitget(uint16(vdata.block(m).status),1:16);
    vdata.block(m).index=fread(fileid_fid,1,'int16'); 
    vdata.block(m).mode=fread(fileid_fid,1,'int16'); 
    vdata.block(m).ctcount=fread(fileid_fid,1,'int32');
    vdata.block(m).lpval=fread(fileid_fid,1,'float32'); 
    vdata.block(m).rpval=fread(fileid_fid,1,'float32'); 
    vdata.block(m).lvl=fread(fileid_fid,1,'float32');
    vdata.block(m).tlt=fread(fileid_fid,1,'float32'); 
    
    % Read in block fid
    if vdata.block(m).bitstatus(4)==1
        fid_data=fread(fileid_fid,vdata.np,'float32');
    elseif vdata.block(m).bitstatus(3)==1
        fid_data=fread(fileid_fid,vdata.np,'int32');
    elseif vdata.block(m).bitstatus(3)==0
        fid_data=fread(fileid_fid,vdata.np,'int16');
    else
        error('illegal combination in file header status.')
    end
    
    % Add to the fid array
    vdata.fid(:,m)=fid_data(1:2:vdata.np,:)+1i*fid_data(2:2:vdata.np);
    
end

% Close the fid file
fclose(fileid_fid);

end

% Consistency enforcement
function grumble(inpath)
if isempty(inpath)
    error('the path string cannot be empty');
end
if ~ischar(inpath)
    error('inpath must be a character string.');
end
end

% I do know how one must live, but I don't 
% want to live like that.
%
% Dmitry Bykov

