% Imports time-domain NMR data recorded by Varian and Agilent inst-
% ruments: reads the binary fid file and the procpar parameter file
% from the experiment directory. Syntax:
%
%                      vdata=v2spinach(inpath)
%
% Parameters:
%
%     inpath  -  character string with the path to the experiment
%                directory containing fid and procpar files
%
% Outputs:
%
%     vdata.fid         - matrix of complex free induction decays,
%                         one per column, in the order they appear
%                         in the file: traces within a data block
%                         run first, data blocks run second
%
%     vdata.procpar     - structure with every parameter found in
%                         the procpar file: numeric parameters as
%                         column vectors, string parameters as
%                         character strings or cell arrays thereof
%
%     vdata.dirname     - experiment directory name
%
%     vdata.nblocks     - number of data blocks in the fid file
%
%     vdata.ntraces     - number of traces in each data block
%
%     vdata.np          - stored points per trace, real and imagi-
%                         nary parts counted separately
%
%     vdata.ebytes      - bytes per stored data point
%
%     vdata.tbytes      - bytes per trace
%
%     vdata.bbytes      - bytes per data block
%
%     vdata.version_id  - VnmrJ version identifier
%
%     vdata.status      - fid file status bitfield
%
%     vdata.nbheaders   - number of block headers per data block
%
%     vdata.block       - block header array with fields: scale,
%                         status, bitstatus, index, mode, ctcount,
%                         lpval, rpval, lvl, and tlt
%
%     vdata.npoints     - complex points per trace
%
%     vdata.arraydim    - number of elements in the parameter array
%
%     vdata.sfrq        - spectrometer frequency, MHz
%
%     vdata.at          - acquisition time, seconds
%
%     vdata.sw_ppm      - spectral width, ppm
%
%     vdata.spec_start  - lower edge of the spectrum, ppm
%
%     vdata.ni          - increment count of the indirect dimensi-
%                         on, present when procpar specifies it
%
%     vdata.grad_amps   - diffusion gradient amplitudes, T/m, pre-
%                         sent when procpar holds a gradient array
%                         and a gradient calibration factor
%
%     vdata.big_delta   - diffusion delay, seconds, present when
%                         procpar specifies it
%
%     vdata.small_delta - diffusion encoding gradient duration,
%                         seconds, present when procpar specifies it
%
%     vdata.gamma       - magnetogyric ratio used by the diffusion
%                         sequence, rad/(s*T), present when procpar
%                         specifies it
%
%     vdata.dosy_const  - Stejskal-Tanner time factor: gamma^2 mul-
%                         tiplied by the effective delta^2*(DELTA-
%                         delta/3) term of the pulse sequence, pre-
%                         sent when procpar specifies it; the sig-
%                         nal attenuation in a pulsed field gradi-
%                         ent experiment is exp(-dosy_const*g^2*D)
%                         where g is the gradient amplitude in T/m
%                         and D is the diffusion coefficient
%                         in m^2/s
%
% Adapted from the varianimport() function of the GNAT package by:
%
%   Dr. Mathias Nilsson
%   School of Chemistry, University of Manchester,
%   Oxford Road, Manchester M13 9PL, UK
%   mathias.nilsson@manchester.ac.uk
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=v2spinach.m>

function vdata=v2spinach(inpath)

% Check consistency
grumble(inpath);

% Store the directory name
vdata.dirname=inpath;

% Open the fid file, big-endian
fid_file=fopen([inpath filesep 'fid'],'r','b');

% Read the fid file header
vdata.nblocks=fread(fid_file,1,'int32');
vdata.ntraces=fread(fid_file,1,'int32');
vdata.np=fread(fid_file,1,'int32');
vdata.ebytes=fread(fid_file,1,'int32');
vdata.tbytes=fread(fid_file,1,'int32');
vdata.bbytes=fread(fid_file,1,'int32');
vdata.version_id=fread(fid_file,1,'int16');
vdata.status=fread(fid_file,1,'int16');
vdata.nbheaders=fread(fid_file,1,'int32');

% Check the header for internal consistency
if (vdata.tbytes~=vdata.np*vdata.ebytes)||...
   (vdata.bbytes~=28*vdata.nbheaders+vdata.ntraces*vdata.tbytes)
    error('inconsistent trace and block sizes in the fid file header.');
end

% Preallocate the fid array
vdata.fid=zeros(vdata.np/2,vdata.nblocks*vdata.ntraces,'like',1i);

% Preallocate the block header array
vdata.block(vdata.nblocks)=struct('scale',[],'status',[],'bitstatus',[],...
                                  'index',[],'mode',[],'ctcount',[],...
                                  'lpval',[],'rpval',[],'lvl',[],'tlt',[]);

% Loop over the data blocks
for n=1:vdata.nblocks

    % Read the block header
    vdata.block(n).scale=fread(fid_file,1,'int16');
    vdata.block(n).status=fread(fid_file,1,'int16');
    vdata.block(n).bitstatus=bitget(uint16(vdata.block(n).status),1:16);
    vdata.block(n).index=fread(fid_file,1,'int16');
    vdata.block(n).mode=fread(fid_file,1,'int16');
    vdata.block(n).ctcount=fread(fid_file,1,'int32');
    vdata.block(n).lpval=fread(fid_file,1,'float32');
    vdata.block(n).rpval=fread(fid_file,1,'float32');
    vdata.block(n).lvl=fread(fid_file,1,'float32');
    vdata.block(n).tlt=fread(fid_file,1,'float32');

    % Skip hypercomplex block headers
    fseek(fid_file,28*(vdata.nbheaders-1),'cof');

    % Get the data type from the block status bits
    if vdata.block(n).bitstatus(4)==1
        data_type='float32';
    elseif vdata.block(n).bitstatus(3)==1
        data_type='int32';
    else
        data_type='int16';
    end

    % Read the data points of the current block
    points=fread(fid_file,vdata.np*vdata.ntraces,data_type);

    % Make sure the block was complete
    if numel(points)<vdata.np*vdata.ntraces
        error('unexpected end of the fid file.');
    end

    % Assemble complex fids, one trace per column
    points=reshape(points,[vdata.np vdata.ntraces]);
    vdata.fid(:,(n-1)*vdata.ntraces+(1:vdata.ntraces))=...
                             points(1:2:end,:)+1i*points(2:2:end,:);

end

% Close the fid file
fclose(fid_file);

% Open the procpar file
par_file=fopen([inpath filesep 'procpar'],'rt');

% Loop over the parameters
while true

    % Read the parameter description line
    name_line=fgetl(par_file);

    % Exit the loop at the end of the file
    if ~ischar(name_line), break; end

    % Skip stray empty lines
    if isempty(strtrim(name_line)), continue; end

    % Get the parameter name and the basic type
    descr=textscan(name_line,'%s %f %f',1);
    par_name=descr{1}{1}; basic_type=descr{3};

    % Read the value line
    val_line=fgetl(par_file);

    % Branch on the basic type
    if basic_type==1

        % Read the value count and the values
        numbers=sscanf(val_line,'%f');
        val_count=numbers(1); par_vals=numbers(2:end);

        % Read further lines if the values continue
        while numel(par_vals)<val_count
            val_line=fgetl(par_file);
            par_vals=[par_vals; sscanf(val_line,'%f')]; %#ok<AGROW>
        end

        % Store the numeric parameter
        vdata.procpar.(par_name)=par_vals;

    elseif basic_type==2

        % Read the value count
        val_count=sscanf(val_line,'%d',1);

        % Loop over the string values
        str_vals=cell(val_count,1);
        for k=1:val_count

            % Move on to the next line
            if k>1, val_line=fgetl(par_file); end

            % Append lines until the closing quote arrives
            while nnz(val_line=='"')<2
                next_line=fgetl(par_file);
                if ~ischar(next_line)
                    error('unexpected end of the procpar file.');
                end
                val_line=[val_line newline next_line]; %#ok<AGROW>
            end

            % Extract the text between the outermost quotes
            quote_pos=find(val_line=='"');
            str_vals{k}=val_line((quote_pos(1)+1):(quote_pos(end)-1));

        end

        % Store single strings as characters, string arrays as cells
        if val_count==1
            vdata.procpar.(par_name)=str_vals{1};
        else
            vdata.procpar.(par_name)=str_vals;
        end

    else

        % Complain and bomb out
        error(['unknown basic type in procpar parameter ' par_name]);

    end

    % Skip the enumeration line
    fgetl(par_file);

end

% Close the procpar file
fclose(par_file);

% Complex points per trace
vdata.npoints=vdata.np/2;

% Array dimension
vdata.arraydim=vdata.procpar.arraydim;

% Spectrometer frequency, MHz
vdata.sfrq=vdata.procpar.sfrq;

% Acquisition time, seconds
vdata.at=vdata.procpar.at;

% Spectral width and lower edge, ppm
vdata.sw_ppm=vdata.procpar.sw/vdata.procpar.sfrq;
vdata.spec_start=(vdata.procpar.rfp-vdata.procpar.rfl)/vdata.procpar.sfrq;

% Indirect dimension increment count when present
if isfield(vdata.procpar,'ni')
    vdata.ni=vdata.procpar.ni;
end

% Gradient amplitudes in T/m when calibration data is present
if isfield(vdata.procpar,'gzlvl1')&&isfield(vdata.procpar,'gcal_')
    vdata.grad_amps=0.01*vdata.procpar.gcal_*vdata.procpar.gzlvl1;
elseif isfield(vdata.procpar,'gzlvl1')&&isfield(vdata.procpar,'DAC_to_G')
    vdata.grad_amps=0.01*vdata.procpar.DAC_to_G*vdata.procpar.gzlvl1;
end

% Diffusion delay and gradient duration when present
if isfield(vdata.procpar,'del'), vdata.big_delta=vdata.procpar.del; end
if isfield(vdata.procpar,'gt1'), vdata.small_delta=vdata.procpar.gt1; end

% Magnetogyric ratio and Stejskal-Tanner factor when present
if isfield(vdata.procpar,'dosygamma')
    vdata.gamma=vdata.procpar.dosygamma;
    if isfield(vdata.procpar,'dosytimecubed')
        vdata.dosy_const=(vdata.procpar.dosygamma^2)*vdata.procpar.dosytimecubed;
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
if ~isfile([inpath filesep 'fid'])
    error(['no fid file in ' inpath]);
end
if ~isfile([inpath filesep 'procpar'])
    error(['no procpar file in ' inpath]);
end
end

% I do know how one must live, but I don't
% want to live like that.
%
% Dmitry Bykov

