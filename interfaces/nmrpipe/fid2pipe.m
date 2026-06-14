% Exports phase-sensitive 2D Spinach free induction decays into native
% NMRPipe time-domain files. Syntax:
%
%          fid2pipe(spin_system,file_root,fid,parameters,nmrpipe_root)
%
% Parameters:
%
%    spin_system  - Spinach spin system structure
%
%    file_root    - output file name root without shell metacharacters;
%                   the function writes <file_root>_pos.fid and
%                   <file_root>_neg.fid
%
%    fid          - structure with fid.pos and fid.neg complex matrices;
%                   rows are direct-dimension points, and columns are
%                   indirect-dimension points
%
%    parameters   - pulse sequence parameters structure with fields
%                   sweep, npoints, offset, and spins
%
%    nmrpipe_root - NMRPipe installation root containing com/txt2pipe.tcl
%                   and nmrbin.linux212_64
%
% Outputs:
%
%    this function writes two NMRPipe files
%
% Note: the exported files preserve the two Spinach echo/anti-echo
%       components; recombination is done in the accompanying NMRPipe
%       processing script after the direct-dimension Fourier transform
%
% ilya.kuprov@weizmann.ac.il

function fid2pipe(spin_system,file_root,fid,parameters,nmrpipe_root)

% Check consistency
grumble(spin_system,file_root,fid,parameters,nmrpipe_root);

% Locate the NMRPipe program directories
nmrbin=fullfile(nmrpipe_root,'nmrbin.linux212_64');
nmrcom=fullfile(nmrpipe_root,'com');
nmrtcl=fullfile(nmrpipe_root,'nmrtcl');

% Set the NMRPipe runtime environment
setenv('PATH',[nmrbin pathsep nmrcom pathsep getenv('PATH')]);
setenv('LD_LIBRARY_PATH',[fullfile(nmrbin,'lib') pathsep getenv('LD_LIBRARY_PATH')]);
setenv('TCLPATH',nmrcom);
setenv('TCL_LIBRARY',fullfile(nmrtcl,'tcl8.4'));
setenv('TK_LIBRARY',fullfile(nmrtcl,'tk8.4'));
setenv('BLT_LIBRARY',fullfile(nmrtcl,'blt2.4'));
setenv('NMRBASE',nmrpipe_root);
setenv('NMRBIN',nmrbin);
setenv('NMRTXT',fullfile(nmrpipe_root,'nmrtxt'));

% Set the OpenWindows path if it is present
if isfolder(fullfile(nmrbin,'openwin'))
    setenv('OPENWINHOME',fullfile(nmrbin,'openwin'));
end


% Convert spectrometer frequencies into MHz
obs_f2=abs(spin(parameters.spins{2})*spin_system.inter.magnet/(2*pi))*1e-6;
obs_f1=abs(spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi))*1e-6;

% Convert carrier offsets into signed ppm
car_f2=hz2ppm(parameters.offset(2),spin_system.inter.magnet,parameters.spins{2});
car_f1=hz2ppm(parameters.offset(1),spin_system.inter.magnet,parameters.spins{1});

% Get acquisition dimensions
npts_f2=size(fid.pos,1);
npts_f1=size(fid.pos,2);

% Export echo and anti-echo components separately
branches={'pos','neg'};
for n=1:numel(branches)

    % Select the current component
    data=fid.(branches{n});

    % Write the text input for txt2pipe
    txt_file=[file_root '_' branches{n} '.txt'];
    file_id=fopen(txt_file,'w');
    if file_id<0
        error('could not open the temporary text file.');
    end
    for k=1:npts_f1
        for m=1:npts_f2
            fprintf(file_id,'%d %d %24.16E %24.16E\n',m,k,real(data(m,k)),imag(data(m,k)));
        end
    end
    fclose(file_id);

    % Convert the text file into NMRPipe format
    pipe_file=[file_root '_' branches{n} '.fid'];
    command=[fullfile(nmrcom,'txt2pipe.tcl') ' -in ' txt_file ...
             ' -xy -complex -time -xN ' num2str(2*npts_f2) ...
             ' -xT ' num2str(npts_f2) ' -xMODE Complex' ...
             ' -xSW ' num2str(parameters.sweep(2),16) ...
             ' -xOBS ' num2str(obs_f2,16) ...
             ' -xCAR ' num2str(car_f2,16) ...
             ' -xLAB ' parameters.spins{2} ...
             ' -yN ' num2str(npts_f1) ...
             ' -yT ' num2str(npts_f1) ' -yMODE Complex' ...
             ' -ySW ' num2str(parameters.sweep(1),16) ...
             ' -yOBS ' num2str(obs_f1,16) ...
             ' -yCAR ' num2str(car_f1,16) ...
             ' -yLAB ' parameters.spins{1} ...
             ' -ndim 2 -aq2D States > ' pipe_file];
    [status,cmdout]=system(command);

    % Remove the temporary text file
    delete(txt_file);

    % Report NMRPipe conversion failures
    if status~=0
        error(['NMRPipe conversion failed: ' cmdout]);
    end

end

end

% Consistency enforcement
function grumble(spin_system,file_root,fid,parameters,nmrpipe_root)
if ~isstruct(spin_system)
    error('spin_system must be a Spinach spin system structure.');
end
if (~isfield(spin_system,'inter'))||(~isfield(spin_system.inter,'magnet'))||...
   (~isnumeric(spin_system.inter.magnet))||(~isreal(spin_system.inter.magnet))||...
   (numel(spin_system.inter.magnet)~=1)||(~isfinite(spin_system.inter.magnet))||...
   (spin_system.inter.magnet==0)
    error('spin_system.inter.magnet must be a non-zero real scalar.');
end
if (~isfield(spin_system,'comp'))||(~isfield(spin_system.comp,'isotopes'))||...
   (~iscell(spin_system.comp.isotopes))
    error('spin_system.comp.isotopes must be a cell array.');
end
if (~ischar(file_root))||isempty(file_root)||any(isspace(file_root))||...
   any((~isstrprop(file_root,'alphanum'))&(~ismember(file_root,'/_.-')))
    error('file_root must be a non-empty character string without shell metacharacters.');
end
[out_dir,~,~]=fileparts(file_root);
if (~isempty(out_dir))&&(~isfolder(out_dir))
    error('file_root directory must exist.');
end
if (~isstruct(fid))||(~isfield(fid,'pos'))||(~isfield(fid,'neg'))||...
   (~isnumeric(fid.pos))||(~isnumeric(fid.neg))||(~ismatrix(fid.pos))||...
   (~ismatrix(fid.neg))||(~all(size(fid.pos)==size(fid.neg)))||...
   any(~isfinite(fid.pos(:)))||any(~isfinite(fid.neg(:)))
    error('fid must contain finite numeric fid.pos and fid.neg matrices of the same size.');
end
if (~isstruct(parameters))||(~isfield(parameters,'sweep'))||...
   (~isfield(parameters,'npoints'))||(~isfield(parameters,'offset'))||...
   (~isfield(parameters,'spins'))
    error('parameters must contain sweep, npoints, offset, and spins fields.');
end
if (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
   (numel(parameters.sweep)~=2)||any(~isfinite(parameters.sweep))||...
   any(parameters.sweep<=0)
    error('parameters.sweep must contain two positive real numbers.');
end
if (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
   (numel(parameters.npoints)~=2)||any(~isfinite(parameters.npoints))||...
   any(parameters.npoints<1)||any(mod(parameters.npoints,1)~=0)
    error('parameters.npoints must contain two positive integers.');
end
if (parameters.npoints(1)~=size(fid.pos,2))||...
   (parameters.npoints(2)~=size(fid.pos,1))
    error('parameters.npoints must match the fid.pos and fid.neg dimensions.');
end
if (~isnumeric(parameters.offset))||(~isreal(parameters.offset))||...
   (numel(parameters.offset)~=2)||any(~isfinite(parameters.offset))
    error('parameters.offset must contain two finite real numbers.');
end
if (~iscell(parameters.spins))||(numel(parameters.spins)~=2)||...
   (~ischar(parameters.spins{1}))||(~ischar(parameters.spins{2}))||...
   any(~ismember(parameters.spins,spin_system.comp.isotopes))
    error('parameters.spins must contain two isotope names present in spin_system.');
end
if (~ischar(nmrpipe_root))||isempty(nmrpipe_root)||any(isspace(nmrpipe_root))||...
   any((~isstrprop(nmrpipe_root,'alphanum'))&(~ismember(nmrpipe_root,'/_.-')))||...
   (~isfolder(nmrpipe_root))||...
   (~isfile(fullfile(nmrpipe_root,'com','txt2pipe.tcl')))||...
   (~isfile(fullfile(nmrpipe_root,'nmrbin.linux212_64','nmrPipe')))
    error('nmrpipe_root must be a valid NMRPipe installation root without shell metacharacters.');
end
end


