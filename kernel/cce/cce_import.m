% Import versioned CCE results, spin tables and Cartesian tensors. Syntax:
%
%                 imported=cce_import(filename,format,units)
%
% format: mat/json (cce_export envelopes returning model/result),
% csv (export sidecar or table
% id,isotope,x,y,z), xyz, tensor (CSV id1,id2,value or id1,id2,xx,...,zz),
% or cube. Scalars denote isotropic tensors.
% External coordinates require units.length='angstrom','nm','m','bohr'.
% Tensors require units.tensor='Hz','kHz','MHz','rad/s','T','mT','G'; field
% units also require units.gyro_rad_s_T. units.frame is mandatory ('lab'
% or a descriptive input frame); optional rotation maps into laboratory
% axes, r_lab=R*r_input and T_lab=R*T_input*R'. Canonical output is Angstrom
% and Hz. Cube output retains density samples and explicit voxel geometry;
% negative atom counts/orbital cubes are rejected rather than guessed.
%
% Inputs:
%   filename, format, units - input path; mat/json/csv/xyz/tensor/cube;
%   external length unit angstrom/nm/m/bohr, tensor unit
%   Hz/kHz/MHz/rad/s/T/mT/G. Field units also require gamma_rad_per_s_T;
%   optional frame is a proper rotation matrix.
%
% Outputs:
%   imported - recovered model/result or canonical spin/tensor/cube data with
%   explicit metadata.
%
% See also cce_export, cce_bath.
%
% ilya.kuprov@weizmann.ac.il

function imported=cce_import(filename,format,units)

% Validate external unit declarations before parsing physical quantities
grumble(filename,format,units);

% Recover the complete versioned result where the file carries its schema
if ismember(format,{'json','mat'})||...
   (strcmp(format,'csv')&&isfile([filename '.json']))
    if strcmp(format,'mat')
        loaded=load(filename,'envelope'); envelope=loaded.envelope;
    else
        if strcmp(format,'csv'), filename=[filename '.json']; end
        envelope=unpack(jsondecode(fileread(filename)));
    end
    if ~isstruct(envelope)||~all(isfield(envelope,{'schema','version','payload'}))||...
       ~strcmp(envelope.schema,'spinach-cce')||envelope.version~=1||...
       ~isfield(envelope,'model')||~isfield(envelope,'payload')
        error('cce_import:version','Unsupported CCE file schema.');
    end
    imported=struct('model',envelope.model,'result',envelope.payload); return;
end

% Require a frame for all non-self-describing physical input
if ~isfield(units,'frame')||~ischar(units.frame)||isempty(units.frame)
    error('External physical input requires a declared frame.');
end
rotation=eye(3);
if isfield(units,'rotation'), rotation=units.rotation; end
imported=struct('schema_version',1,'source',filename,'input_units',units,...
                'frame','lab','rotation',rotation);

% Parse strict tables rather than inferring isotope names or frequency units
switch format
    case 'xyz'
        fid=fopen(filename,'r'); cleanup=onCleanup(@()fclose(fid));
        n=str2double(fgetl(fid)); imported.comment=fgetl(fid);
        if ~isscalar(n)||~isfinite(n)||n<0||n~=fix(n), error('Invalid XYZ site count.'); end
        data=textscan(fid,'%s %f %f %f','CollectOutput',true);
        if numel(data{1})~=n||size(data{2},1)~=n, error('XYZ count mismatch.'); end
        imported.labels=data{1}.'; imported.ids=1:n;
        imported.coordinates_A=data{2}*length_factor(units)*rotation.';
    case 'csv'
        data=readtable(filename,'TextType','string');
        if ~all(ismember({'id','isotope','x','y','z'},data.Properties.VariableNames))
            error('Spin CSV requires id,isotope,x,y,z.');
        end
        imported.ids=data.id.'; imported.isotopes=cellstr(data.isotope).';
        imported.coordinates_A=[data.x data.y data.z]*length_factor(units)*rotation.';
        if numel(unique(imported.ids))~=height(data), error('Duplicate physical IDs.'); end
    case 'tensor'
        data=readmatrix(filename); data=data(~all(isnan(data),2),:);
        if ~ismember(size(data,2),[3 11])||any(~isfinite(data),'all')
            error('Tensor CSV needs id1,id2,scalar or nine row-major components.');
        end
        imported.pairs_Hz=struct('ids',{},'tensor',{}); scale=tensor_factor(units);
        for k=1:size(data,1)
            if size(data,2)==3, tensor=data(k,3)*eye(3);
            else, tensor=reshape(data(k,3:11),3,3).'; end
            imported.pairs_Hz(k)=struct('ids',data(k,1:2),...
                                       'tensor',scale*rotation*tensor*rotation.');
        end
    case 'cube'
        fid=fopen(filename,'r'); cleanup=onCleanup(@()fclose(fid));
        imported.comment={fgetl(fid),fgetl(fid)};
        line=sscanf(fgetl(fid),'%f').'; atoms=line(1);
        if atoms<0||atoms~=fix(atoms), error('Orbital cubes are not supported.'); end
        imported.origin_A=line(2:4)*length_factor(units)*rotation.';
        counts=zeros(1,3); vectors=zeros(3);
        for k=1:3
            line=sscanf(fgetl(fid),'%f').'; counts(k)=line(1); vectors(k,:)=line(2:4);
        end
        if any(counts<=0|counts~=fix(counts)), error('Cube counts must be positive; declare units explicitly.'); end
        atom_table=zeros(atoms,5);
        for k=1:atoms, atom_table(k,:)=sscanf(fgetl(fid),'%f').'; end
        values=fscanf(fid,'%f');
        if numel(values)~=prod(counts), error('Cube voxel count mismatch.'); end
        imported.atomic_numbers=atom_table(:,1); imported.charges=atom_table(:,2);
        imported.coordinates_A=atom_table(:,3:5)*length_factor(units)*rotation.';
        imported.voxel_vectors_A=vectors*length_factor(units)*rotation.';
        imported.grid_size=counts; imported.density=permute(reshape(values,fliplr(counts)),[3 2 1]);
        if ~isfield(units,'density'), error('Cube density units must be declared.'); end
        imported.density_units=units.density;
end

% Reject corrupt coordinates rather than silently dropping sites
if isfield(imported,'coordinates_A')&&any(~isfinite(imported.coordinates_A),'all')
    error('Nonfinite imported coordinates.');
end

end

function factor=length_factor(units)

% Convert explicitly declared lengths into the native Spinach boundary
if ~isfield(units,'length'), error('Length units are required.'); end
switch units.length
    case 'angstrom', factor=1;
    case 'nm', factor=10;
    case 'm', factor=1e10;
    case 'bohr', factor=0.529177210544;
    otherwise, error('Unsupported length units.');
end

end

function factor=tensor_factor(units)

% Convert frequency or explicitly gyromagnetic field tensors into Hz
if ~isfield(units,'tensor'), error('Tensor units are required.'); end
switch units.tensor
    case 'Hz', factor=1;
    case 'kHz', factor=1e3;
    case 'MHz', factor=1e6;
    case 'rad/s', factor=1/(2*pi);
    case {'T','mT','G'}
        if ~isfield(units,'gyro_rad_s_T')||~isscalar(units.gyro_rad_s_T)||...
           ~isfinite(units.gyro_rad_s_T), error('Field tensors need explicit rad/s/T gyro.'); end
        scale=1;
        if strcmp(units.tensor,'mT'), scale=1e-3; end
        if strcmp(units.tensor,'G'), scale=1e-4; end
        factor=scale*units.gyro_rad_s_T/(2*pi);
    otherwise, error('Unsupported tensor units.');
end

end

function value=unpack(encoded)

% Reconstitute each numeric component and the original MATLAB array shape
switch encoded.kind
    case 'integer'
        value=reshape(typecast(uint8(encoded.bytes),encoded.class),encoded.shape.');
    case 'array'
        if prod(encoded.shape)==0
            value=zeros(encoded.shape.',encoded.class); return;
        end
        re=str2double(encoded.real); im=str2double(encoded.imag);
        if all(im==0), data=re; else, data=complex(re,im); end
        value=reshape(cast(data,encoded.class),encoded.shape.');
    case 'char'
        value=encoded.value;
    case 'string'
        value=reshape(string(encoded.items),encoded.shape.');
    case 'cell'
        items=encoded.items;
        if isstruct(items), items=num2cell(items); end
        value=reshape(cellfun(@unpack,items,'UniformOutput',false),encoded.shape.');
    case 'struct'
        names=encoded.fields; items=encoded.items;
        if isstruct(items), items=num2cell(items); end
        value=repmat(cell2struct(cell(size(names)),names),encoded.shape.');
        for j=1:numel(names)
            for k=1:numel(value), value(k).(names{j})=unpack(items{k+(j-1)*numel(value)}); end
        end
    otherwise
        error('Unknown serialization type.');
end

end

function grumble(filename,format,units)

if ~ischar(filename)||~isfile(filename)||~ischar(format)||...
   ~ismember(format,{'mat','json','csv','xyz','tensor','cube'})||...
   ~isstruct(units)||~isscalar(units)
    error('Invalid CCE import arguments.');
end
if isfield(units,'rotation')&&(~isequal(size(units.rotation),[3 3])||...
   any(~isfinite(units.rotation),'all')||norm(units.rotation.'*units.rotation-eye(3),'fro')>100*eps||...
   abs(det(units.rotation)-1)>100*eps)
    error('Invalid frame rotation.');
end

end

