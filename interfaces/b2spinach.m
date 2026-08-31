% Imports time-domain NMR data recorded by Bruker instruments:
% reads the binary fid or ser file together with the acquisiti-
% on and processing parameter files from the numbered experi-
% ment directory. Syntax:
%
%                      bdata=b2spinach(inpath)
%
% Parameters:
%
%     inpath  -  character string with the path to the numbered
%                Bruker experiment directory containing the ac-
%                qus file and the fid or ser file
%
% Outputs:
%
%     bdata.fid         - matrix of complex free induction de-
%                         cays, one per column, in the order
%                         they are stored in the file; for data
%                         sets with three and more dimensions
%                         the loop order over the indirect di-
%                         mensions is determined by the AQSEQ
%                         parameter of the acqus structure
%
%     bdata.acqus       - structure with every parameter found
%                         in the acqus file: numeric parameters
%                         as scalars or column vectors, string
%                         parameters as character strings or
%                         cell arrays thereof
%
%     bdata.acqu2s      - same for the acqu2s file of data sets
%                         with two or more dimensions
%
%     bdata.acqu3s      - same for the acqu3s file of data sets
%                         with three or more dimensions
%
%     bdata.acqu4s      - same for the acqu4s file of four-di-
%                         mensional data sets
%
%     bdata.procs       - same for the procs file of the first
%                         processed data directory when present
%
%     bdata.dirname     - experiment directory name
%
%     bdata.ndims_data  - number of dimensions in the data set
%
%     bdata.arraydim    - number of fids the status parameter
%                         files declare as acquired
%
%     bdata.fids_in_file - number of fid slots held by the bina-
%                         ry file, and the column count of the
%                         fid matrix; exceeds arraydim when the
%                         acquisition was preallocated or inter-
%                         rupted, and falls short of it for non-
%                         uniformly sampled data sets, where the
%                         fids follow the acquisition order of
%                         the sampling schedule in nus_list
%
%     bdata.npoints     - complex points per fid
%
%     bdata.pulprog     - pulse programme name
%
%     bdata.nucleus     - observe nucleus, e.g. '1H'
%
%     bdata.gamma       - magnetogyric ratio of the observe
%                         nucleus, rad/(s*T)
%
%     bdata.sfrq        - spectrometer frequency, MHz
%
%     bdata.at          - acquisition time, seconds
%
%     bdata.sw_ppm      - spectral width, ppm
%
%     bdata.spec_start  - lower edge of the spectrum, ppm
%
%     bdata.digshift    - group delay of the digital filter in
%                         complex points; the first round(dig-
%                         shift) points of each fid precede the
%                         start of the true signal
%
%     bdata.grad_amps   - gradient amplitudes, T/m, present
%                         when the difflist file exists
%
%     bdata.vd_list     - variable delay list, seconds, present
%                         when the vdlist file exists
%
%     bdata.vc_list     - variable counter list, present when
%                         the vclist file exists
%
%     bdata.nus_list    - non-uniform sampling schedule, one
%                         row per sampled fid and one column
%                         per indirect dimension, present
%                         when the nuslist file exists
%
% Adapted from the brukerimport() function of the GNAT package by:
%
%   Dr. Mathias Nilsson
%   School of Chemistry, University of Manchester,
%   Oxford Road, Manchester M13 9PL, UK
%   mathias.nilsson@manchester.ac.uk
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=b2spinach.m>

function bdata=b2spinach(inpath)

% Check consistency
grumble(inpath);

% Store the directory name
bdata.dirname=inpath;

% Refuse data with more than four dimensions
if isfile([inpath filesep 'acqu5s'])
    error('data sets with more than four dimensions are not supported.');
end

% Detect data dimension from the status parameter files
if isfile([inpath filesep 'acqu4s'])
    bdata.ndims_data=4;
elseif isfile([inpath filesep 'acqu3s'])
    bdata.ndims_data=3;
elseif isfile([inpath filesep 'acqu2s'])
    bdata.ndims_data=2;
else
    bdata.ndims_data=1;
end

% Read acquisition status parameters for every dimension
bdata.acqus=read_jcamp([inpath filesep 'acqus']);
if bdata.ndims_data>=2
    bdata.acqu2s=read_jcamp([inpath filesep 'acqu2s']);
end
if bdata.ndims_data>=3
    bdata.acqu3s=read_jcamp([inpath filesep 'acqu3s']);
end
if bdata.ndims_data>=4
    bdata.acqu4s=read_jcamp([inpath filesep 'acqu4s']);
end

% Read processing parameters when present
if isfile([inpath filesep 'pdata' filesep '1' filesep 'procs'])
    bdata.procs=read_jcamp([inpath filesep 'pdata' filesep '1' filesep 'procs']);
end

% Count the fids in the data set
bdata.arraydim=1;
if bdata.ndims_data>=2, bdata.arraydim=bdata.arraydim*bdata.acqu2s.TD; end
if bdata.ndims_data>=3, bdata.arraydim=bdata.arraydim*bdata.acqu3s.TD; end
if bdata.ndims_data>=4, bdata.arraydim=bdata.arraydim*bdata.acqu4s.TD; end

% Complex points per fid
bdata.npoints=bdata.acqus.TD/2;

% Pulse programme name and observe nucleus
bdata.pulprog=bdata.acqus.PULPROG;
bdata.nucleus=bdata.acqus.NUC1;

% Magnetogyric ratio of the observe nucleus
bdata.gamma=spin(bdata.nucleus);

% Spectrometer frequency, MHz
bdata.sfrq=bdata.acqus.SFO1;

% Spectral width, ppm
bdata.sw_ppm=bdata.acqus.SW;

% Acquisition time, seconds
bdata.at=bdata.npoints/(bdata.sw_ppm*bdata.sfrq);

% Spectrum lower edge from processing or acquisition parameters
if isfield(bdata,'procs')
    bdata.spec_start=bdata.procs.OFFSET-bdata.sw_ppm;
else
    bdata.spec_start=1e6*(bdata.acqus.SFO1/bdata.acqus.BF1-1)-...
                     0.5*bdata.sw_ppm*bdata.acqus.SFO1/bdata.acqus.BF1;
end

% Digital filter group delay in complex points
if isfield(bdata.acqus,'DIGMOD')&&(bdata.acqus.DIGMOD==0)

    % Analogue filtering has no group delay
    bdata.digshift=0;

elseif isfield(bdata.acqus,'GRPDLY')&&(bdata.acqus.GRPDLY~=-1)

    % Modern hardware reports the group delay directly
    bdata.digshift=bdata.acqus.GRPDLY;

elseif isfield(bdata.acqus,'DECIM')&&(bdata.acqus.DECIM~=1)

    % Decimation factors of older DSP firmware
    decims=[2 3 4 6 8 12 16 24 32 48 64 96 128 192 ...
            256 384 512 768 1024 1536 2048];

    % Group delay table, one row per DECIM, columns DSPFVS 10 to 12
    delays=[44.7500000000  46.000  46.311; 33.5000000000  36.500  36.530;
            66.6250000000  48.000  47.870; 59.0833333333  50.167  50.229;
            68.5625000000  53.250  53.289; 60.3750000000  69.500  69.551;
            69.5312500000  72.250  71.600; 61.0208333333  70.167  70.184;
            70.0156250000  72.750  72.138; 61.3437500000  70.500  70.528;
            70.2578125000  73.000  72.348; 61.5052083333  70.667  70.700;
            70.3789062500  72.500  72.524; 61.5859375000  71.333     NaN;
            70.4394531250  72.250     NaN; 61.6263020833  71.667     NaN;
            70.4697265625  72.125     NaN; 61.6464843750  71.833     NaN;
            70.4848632813  72.063     NaN; 61.6565755208  71.917     NaN;
            70.4924316406  72.031     NaN];

    % Make sure the firmware version is present
    if ~isfield(bdata.acqus,'DSPFVS')
        error('DECIM without DSPFVS in acqus: cannot get the group delay.');
    end

    % Locate the decimation factor
    decim_row=find(decims==bdata.acqus.DECIM,1);
    if isempty(decim_row)
        error('unknown DECIM value in acqus.');
    end

    % Branch on the firmware version
    if (bdata.acqus.DSPFVS>=10)&&(bdata.acqus.DSPFVS<=12)

        % Tabulated group delays of the DMX generation
        bdata.digshift=delays(decim_row,bdata.acqus.DSPFVS-9);
        if isnan(bdata.digshift)
            error('unknown DECIM and DSPFVS combination in acqus.');
        end

    elseif bdata.acqus.DSPFVS==13

        % Closed form group delay of the 13 firmware family
        bdata.digshift=3-1/(2*bdata.acqus.DECIM);

    else

        % Complain and bomb out
        error('unknown DSPFVS value in acqus.');

    end

else

    % Data without digital filtering
    bdata.digshift=0;

end

% Binary data byte order
switch bdata.acqus.BYTORDA
    case 0
        byte_order='l';
    case 1
        byte_order='b';
    otherwise
        error('unknown BYTORDA value in acqus.');
end

% Binary data type and file block size in points
switch bdata.acqus.DTYPA
    case 0
        data_type='int32'; block_pts=256;
    case 1
        data_type='float32'; block_pts=256;
    case 2
        data_type='float64'; block_pts=128;
    otherwise
        error('unknown DTYPA value in acqus.');
end

% Stored points per fid, padded to whole 1024-byte blocks
padded_pts=block_pts*ceil(2*bdata.npoints/block_pts);

% Locate the binary data file
if bdata.ndims_data==1
    bin_path=[inpath filesep 'fid'];
else
    bin_path=[inpath filesep 'ser'];
end
if ~isfile(bin_path)
    error(['file not found: ' bin_path]);
end

% Read the binary data
bin_file=fopen(bin_path,'r',byte_order);
data_pts=fread(bin_file,inf,data_type);
fclose(bin_file);

% Get the per-fid stride, padded or unpadded
if mod(numel(data_pts),padded_pts)==0
    fid_stride=padded_pts;
elseif mod(numel(data_pts),2*bdata.npoints)==0
    fid_stride=2*bdata.npoints;
else
    error('the size of the data file does not match the acqus parameters.');
end

% Count the fids held by the file, at least the declared number
bdata.fids_in_file=numel(data_pts)/fid_stride;
if (bdata.fids_in_file<bdata.arraydim)&&(~isfile([inpath filesep 'nuslist']))
    error('the data file holds fewer fids than the status parameters declare.');
end

% Reshape and drop the block padding
data_pts=reshape(data_pts,[fid_stride bdata.fids_in_file]);
data_pts=data_pts(1:(2*bdata.npoints),:);

% Assemble complex fids using the Bruker sign convention
bdata.fid=data_pts(1:2:end,:)-1i*data_pts(2:2:end,:);

% Apply the binary scaling exponent
bdata.fid=(2^bdata.acqus.NC)*bdata.fid;

% Gradient amplitudes in T/m when a difflist file is present
if isfile([inpath filesep 'difflist'])
    bdata.grad_amps=0.01*read_list([inpath filesep 'difflist']);
elseif isfile([inpath filesep 'difflist.txt'])
    bdata.grad_amps=0.01*read_list([inpath filesep 'difflist.txt']);
end

% Variable delay list in seconds when present
if isfile([inpath filesep 'vdlist'])
    bdata.vd_list=read_list([inpath filesep 'vdlist']);
end

% Variable counter list when present
if isfile([inpath filesep 'vclist'])
    bdata.vc_list=read_list([inpath filesep 'vclist']);
end

% Non-uniform sampling schedule when present
if isfile([inpath filesep 'nuslist'])
    bdata.nus_list=read_list([inpath filesep 'nuslist']);
end

end

% Reads a Bruker JCAMP-DX parameter file into a structure
function params=read_jcamp(file_path)

% Read the file as a cell array of lines
all_lines=regexp(fileread(file_path),'\r?\n','split');

% Loop over the lines
params=struct(); line_num=1;
while line_num<=numel(all_lines)

    % Read the next line
    text_line=all_lines{line_num}; line_num=line_num+1;

    % Only private parameter lines are of interest
    if ~strncmp(text_line,'##$',3), continue; end

    % Cut off inline comments
    dollar_pos=strfind(text_line,'$$');
    if ~isempty(dollar_pos)
        text_line=text_line(1:(dollar_pos(1)-1));
    end

    % Split the name from the value
    eq_pos=find(text_line=='=',1);
    par_name=text_line(4:(eq_pos-1));
    par_value=strtrim(text_line((eq_pos+1):end));

    % Assemble multi-line array values
    if ~isempty(regexp(par_value,'^\(\d+\.\.\d+\)$','once'))
        par_value='';
        while line_num<=numel(all_lines)
            next_line=all_lines{line_num};
            if strncmp(next_line,'##',2)||strncmp(next_line,'$$',2)
                break;
            end
            par_value=[par_value ' ' next_line]; line_num=line_num+1; %#ok<AGROW>
        end
        par_value=strtrim(par_value);
    end

    % Branch on the value type
    if any(par_value=='<')

        % Strip the angle brackets from strings
        str_toks=regexp(par_value,'<([^<>]*)>','tokens');
        str_toks=[str_toks{:}];

        % Store single strings as characters, string lists as cells
        if isscalar(str_toks)
            params.(par_name)=str_toks{1};
        else
            params.(par_name)=str_toks;
        end

    else

        % Convert to numbers with text fallback
        [num_value,~,~,next_idx]=sscanf(par_value,'%f');
        if isempty(num_value)||(next_idx<=numel(par_value))
            params.(par_name)=par_value;
        else
            params.(par_name)=num_value;
        end

    end

end

end

% Reads a Bruker list file with optional time unit suffixes
function vals=read_list(file_path)

% Read the file as a cell array of lines
all_lines=regexp(fileread(file_path),'\r?\n','split');

% Drop empty lines
all_lines=strtrim(all_lines);
all_lines=all_lines(~cellfun(@isempty,all_lines));

% Refuse empty list files
if isempty(all_lines)
    error(['list file is empty: ' file_path]);
end

% Split the lines into whitespace-separated entries
all_lines=cellfun(@strsplit,all_lines,'UniformOutput',false);

% Require a consistent entry count per line
n_cols=numel(all_lines{1});
if any(cellfun(@numel,all_lines)~=n_cols)
    error(['inconsistent number of columns in ' file_path]);
end

% Preallocate the values
vals=zeros(numel(all_lines),n_cols);

% Loop over the entries
for n=1:numel(all_lines)
    for k=1:n_cols

        % Convert time unit suffixes into multipliers
        entry=all_lines{n}{k};
        switch entry(end)
            case 'n'
                multiplier=1e-9; entry=entry(1:(end-1));
            case 'u'
                multiplier=1e-6; entry=entry(1:(end-1));
            case 'm'
                multiplier=1e-3; entry=entry(1:(end-1));
            case 's'
                multiplier=1;    entry=entry(1:(end-1));
            otherwise
                multiplier=1;
        end

        % Convert the value into a number
        vals(n,k)=multiplier*str2double(entry);

        % Make sure the conversion succeeded
        if isnan(vals(n,k))
            error(['cannot interpret list file entry: ' all_lines{n}{k}]);
        end

    end
end

end

% Consistency enforcement
function grumble(inpath)
if ~ischar(inpath)
    error('inpath must be a character string.');
end
if ~isfolder(inpath)
    error(['directory ' inpath ' does not exist.']);
end
if ~isfile([inpath filesep 'acqus'])
    error(['no acqus file in ' inpath]);
end
end

% The purpose of computing is insight, not numbers.
%
% Richard Hamming

