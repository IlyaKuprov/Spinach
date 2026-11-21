% Reads a literal spinsys block into native Spinach inputs. Syntax:
%
%                  [sys,inter]=s2spinach(filename)
%
% Parameters:
%
%     filename  - path to a text file containing one literal
%                 spinsys block, as a character row or scalar string
%
% Outputs:
%
%     sys       - Spinach system specification with sys.isotopes
%                 in nuclei record order; supply sys.magnet separately
%
%     inter     - Spinach interaction specification, with chemical
%                 shift tensors in ppm and coupling tensors in Hz
%
% Supports channels, nuclei, shift, dipole, jcoupling, and quadrupole,
% with one record per line and # comments. Channels identify RF
% addressing and are not additional spins. Shift amplitudes require
% a p (ppm) suffix, except for zero. Quadrupolar orders 1 to 3 are
% read as source annotations; the caller chooses the native Hamiltonian
% approximation. No pulse programs or expressions are run.
%
% Chemical shift anisotropy is zz minus isotropic shift. Dipolar and
% anisotropic J inputs are half the corresponding traceless zz value.
% Euler angles are passive ZYZ in degrees, and eta uses yy minus xx.

function [sys,inter]=s2spinach(filename)

% Check consistency
grumble(filename);

% Tokenise the literal block and retain source line numbers
lines=regexp(fileread(filename),'\r\n|\n|\r','split');
records=cell(size(lines));
opened=false;
closed=false;
for n=1:numel(lines)

    % Discard blank lines and comments
    line=strtrim(regexprep(lines{n},'#.*$',''));
    if isempty(line), continue; end

    % Identify the block and its records
    if ~opened&&~isempty(regexp(line,'^spinsys\s*\{$','once'))
        opened=true;
    elseif opened&&~closed&&strcmp(line,'}')
        closed=true;
    elseif opened&&~closed
        records{n}=regexp(line,'\s+','split');
    else
        error('s2spinach:syntax','%s:%d: expected a single literal spinsys block.',filename,n);
    end

end
if ~closed
    error('s2spinach:syntax','%s:%d: missing or unclosed spinsys block.',filename,numel(lines));
end

% Resolve spin numbers from the nuclei record before reading interactions
record_lines=find(~cellfun(@isempty,records));
nuclei=record_lines(cellfun(@(words)strcmp(words{1},'nuclei'),records(record_lines)));
if isempty(nuclei)
    error('s2spinach:syntax','%s:%d: missing nuclei record.',filename,numel(lines));
elseif numel(nuclei)>1
    error('s2spinach:syntax','%s:%d: repeated nuclei record.',filename,nuclei(2));
end
sys.isotopes=records{nuclei}(2:end);
nspins=numel(sys.isotopes);
if nspins==0
    error('s2spinach:syntax','%s:%d: empty nuclei record.',filename,nuclei);
end
mult=zeros(1,nspins);
for n=1:nspins

    % Look up physical nuclear spin multiplicities
    label=sys.isotopes{n};
    if isempty(regexp(label,'^\d+[A-Z][a-z]?(_m\d*)?$','once'))
        error('s2spinach:isotope','%s:%d: invalid nuclear isotope %s.',filename,nuclei,label);
    end
    try
        [~,mult(n)]=spin(label);
    catch err
        if ~ismember(err.identifier,{'spin:unknown_isotope','spin:data_unavailable'})
            rethrow(err);
        end
        error('s2spinach:isotope','%s:%d: %s',filename,nuclei,err.message);
    end
    if mult(n)<2
        error('s2spinach:isotope','%s:%d: isotope %s has no nuclear spin.',filename,nuclei,label);
    end

end

% Allocate native input cells and define the literal record sizes
inter.zeeman.matrix=cell(1,nspins);
inter.coupling.scalar=cell(nspins);
inter.coupling.matrix=cell(nspins);
dipole_seen=false(nspins);
channels_seen=false;
terms={'shift','dipole','jcoupling','quadrupole'};
counts={7,[6 7],8,7};
indices=[1 2 2 1];
literal='^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$';

% Resolve each interaction and accumulate it in the native cells
for n=record_lines

    % Handle non-interaction records
    words=records{n};
    tag=words{1};
    args=words(2:end);
    if strcmp(tag,'nuclei'), continue; end
    if strcmp(tag,'channels')
        if channels_seen||isempty(args)||~all(ismember(args,sys.isotopes))||...
           numel(unique(args))~=numel(args)
            error('s2spinach:channels','%s:%d: expected distinct channels from the nuclei list.',filename,n);
        end
        channels_seen=true;
        continue;
    end

    % Validate the record type and number of arguments
    term=find(strcmp(tag,terms));
    if isempty(term)
        error('s2spinach:unsupported','%s:%d: unsupported record %s.',filename,n,tag);
    elseif ~ismember(numel(args),counts{term})
        error('s2spinach:syntax','%s:%d: wrong number of arguments for %s.',filename,n,tag);
    end

    % Convert numeric literals and require field-independent shift values
    if term==1
        has_ppm=endsWith(args(2:3),'p');
        args(2:3)=regexprep(args(2:3),'p$','');
    end
    values=str2double(args);
    if any(cellfun(@isempty,regexp(args,literal,'once')))||any(~isfinite(values))
        error('s2spinach:syntax','%s:%d: expected finite numeric literals in %s.',filename,n,tag);
    end
    if term==1&&any(~has_ppm&(values(2:3)~=0))
        error('s2spinach:units','%s:%d: nonzero shift amplitudes require the p (ppm) suffix.',filename,n);
    end

    % Resolve one spin or a canonical, off-diagonal spin pair
    spins=sort(values(1:indices(term)));
    if any(spins<1|spins>nspins|spins~=fix(spins))||...
       (numel(spins)==2&&spins(1)==spins(2))
        error('s2spinach:index','%s:%d: invalid spin index or self-pair in %s.',filename,n,tag);
    end
    first=spins(1);
    last=spins(end);

    % Convert passive ZYZ degrees to active ZYZ radians
    eulers=-deg2rad(values(end:-1:end-2));

    % Convert to native units and store each pair only once
    switch tag
        case 'shift'
            if ~isempty(inter.zeeman.matrix{first})
                error('s2spinach:duplicate','%s:%d: repeated shift record for spin %d.',filename,n,first);
            end
            inter.zeeman.matrix{first}=values(2)*eye(3)+...
                                       pas_tensor(values(3),values(4),eulers);
        case 'dipole'
            if dipole_seen(first,last)
                error('s2spinach:duplicate','%s:%d: repeated dipole record.',filename,n);
            end
            dipole_seen(first,last)=true;
            if numel(values)==6
                eta=0;
            else
                eta=values(4);
            end
            pair_tensor=pas_tensor(2*values(3),eta,eulers);
        case 'jcoupling'
            if ~isempty(inter.coupling.scalar{first,last})
                error('s2spinach:duplicate','%s:%d: repeated jcoupling record.',filename,n);
            end
            inter.coupling.scalar{first,last}=values(3);
            pair_tensor=pas_tensor(2*values(4),values(5),eulers);
        case 'quadrupole'
            if ~isempty(inter.coupling.matrix{first,first})
                error('s2spinach:duplicate','%s:%d: repeated quadrupole record for spin %d.',filename,n,first);
            end
            if mult(first)<3||~ismember(values(2),1:3)
                error('s2spinach:quadrupole','%s:%d: quadrupole requires spin >= 1 and order 1, 2, or 3.',filename,n);
            end
            inter.coupling.matrix{first,first}=eeqq2nqi(values(3),-values(4),(mult(first)-1)/2,eulers);
    end

    % Add the traceless dipolar and J tensors without duplicating scalar J
    if numel(spins)==2
        if isempty(inter.coupling.matrix{first,last})
            inter.coupling.matrix{first,last}=pair_tensor;
        else
            inter.coupling.matrix{first,last}=inter.coupling.matrix{first,last}+pair_tensor;
        end
    end

end

end

% Rotate a traceless principal-axis tensor into the molecular frame
function T=pas_tensor(delta,eta,eulers)

% Set the principal values using the signed zz and yy-minus-xx convention
R=euler2dcm(eulers);
T=R*diag([-delta*(1+eta)/2 -delta*(1-eta)/2 delta])*R';

end

% Consistency enforcement
function grumble(filename)
if ~((ischar(filename)&&isrow(filename))||...
     (isstring(filename)&&isscalar(filename)&&~ismissing(filename)))||strlength(filename)==0
    error('filename must be a non-empty character row or scalar string.');
end
if ~isfile(filename)
    error('input file does not exist: %s',filename);
end
end


