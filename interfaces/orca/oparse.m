% A parser for ORCA text output logs, versions 2.6 to 6.1. Reads the
% geometry and every magnetic parameter that ORCA prints in the main
% output file. Syntax:
%
%                        props=oparse(file_name)
%
% Parameters:
%
%    file_name - a character string with the file path
%
% Outputs:
%
%   props.filename          - log file name
%   props.orca_version      - ORCA version string
%   props.symbols           - atomic symbols, 1 x natoms cell
%   props.atomic_numbers    - atomic numbers, 1 x natoms
%   props.std_geom          - atomic coordinates, natoms x 3, Angstrom
%   props.natoms            - number of atoms
%   props.charge            - total charge
%   props.multiplicity      - spin multiplicity
%   props.energy            - final single point energy, Hartree
%   props.dip_moment        - electric dipole moment, a.u.
%   props.mulliken_chg      - Mulliken atomic charges, natoms x 1
%   props.mulliken_spin     - Mulliken spin populations, natoms x 1
%   props.g_tensor.raw      - g-matrix as printed by ORCA
%   props.g_tensor.matrix   - symmetrised g-matrix
%   props.g_tensor.eigvals  - eigenvalues of the symmetrised g-matrix
%   props.g_tensor.eigvecs  - eigenvectors of the symmetrised g-matrix
%   props.zfs.matrix        - zero-field splitting tensor, cm^-1
%   props.zfs.eigvals       - ZFS tensor eigenvalues, cm^-1
%   props.zfs.eigvecs       - ZFS tensor eigenvectors
%   props.hfc.full.matrix   - hyperfine tensors, Gauss, natoms cell
%   props.hfc.full.eigvals  - hyperfine eigenvalues, Gauss, natoms cell
%   props.hfc.full.eigvecs  - hyperfine eigenvectors, natoms cell
%   props.hfc.iso           - isotropic hyperfine couplings, Gauss
%   props.efg               - EFG tensors, a.u.^-3, natoms cell
%   props.nqi               - quadrupolar tensors, Hz, natoms cell
%   props.isotopes          - isotopes used by ORCA, natoms cell
%   props.cst               - shielding tensors, ppm, natoms cell
%   props.j_couplings       - isotropic J-couplings, Hz, natoms x natoms
%   props.chi_temps         - susceptibility temperatures, K
%   props.chi_tensors       - molar magnetic susceptibility tensors,
%                             cm^3*K/mol, one cell per temperature
%
% Only the fields that ORCA has actually printed are returned; the
% caller should test for their presence with isfield.
%
% Notes: ORCA prints magnetic parameters only for the nuclei that were
%        requested in the input, and labels each of them with the zero
%        based index of the atom in the Cartesian coordinate table. All
%        per-atom outputs above are therefore indexed by the position of
%        the atom in props.std_geom, and are left empty for atoms whose
%        parameters were not printed.
%
%        When a log contains multiple geometries or multiple property
%        sections, for example a geometry optimisation or a relaxed
%        surface scan, the last one printed is returned.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=oparse.m>

function props=oparse(file_name)

% Check consistency
grumble(file_name);

% Read the file, preserving the leading whitespace
file_id=fopen(file_name,'r');
orca_log=textscan(file_id,'%s','delimiter','\n','whitespace','');
fclose(file_id); raw_lines=orca_log{1};

% Content matching is done on whitespace-trimmed lines
log_lines=strtrim(raw_lines); nlines=numel(log_lines);
props.filename=file_name;

% ORCA rules its section banners from the first column, and indents
% the rulers that separate the per-nucleus blocks inside a section
ruler_mask=~cellfun(@isempty,regexp(raw_lines,'^[-=]{5,}\s*$','once'));
title_mask=false(nlines,1);
title_mask(2:(nlines-1))=ruler_mask(1:(nlines-2))&ruler_mask(3:nlines)&...
                         (~cellfun(@isempty,log_lines(2:(nlines-1))));
banner_idx=find(title_mask); banner_titles=log_lines(banner_idx);

% Property sections that list nuclei run past the banners of the
% sub-programs that ORCA calls while filling them in, and end only
% where the next property section or the final report begins
prop_idx=banner_idx(startsWith(banner_titles,...
         {'CARTESIAN COORDINATES','ELECTRONIC G-MATRIX','ZERO-FIELD-SPLITTING',...
          'ELECTRIC AND MAGNETIC HYPERFINE','CHEMICAL SHIFTS','CHEMICAL SHIELDINGS',...
          'MULLIKEN ATOMIC CHARGES','LOEWDIN ATOMIC CHARGES','ORCA EULER ANGLE PROGRAM',...
          'ORCA PROPERTY CALCULATIONS','TIMINGS','ORCA TERMINATED NORMALLY'}));

% Read the ORCA version and refuse anything that is not an ORCA log
version_idx=find_lines(log_lines,'^Program Version\s+\d',1,nlines);
if isempty(version_idx)
    error('this file does not contain an ORCA version banner.');
end
props.orca_version=regexp(log_lines{version_idx(1)},'\d+(\.\d+)*','match','once');
version_major=str2double(regexp(props.orca_version,'^\d+','match','once'));

% Select the parser branch: ORCA 6 renamed the hyperfine and EFG
% blocks from "Raw" to "Total", and the shielding section from
% CHEMICAL SHIFTS to CHEMICAL SHIELDINGS. The spelling used by the
% other branch is retained as a fallback, so that patched builds
% and pre-release versions are still parsed
if version_major>=6
    hfc_hdr='^Total HFC matrix'; hfc_alt='^Raw HFC matrix';
    efg_hdr='^Total EFG matrix'; efg_alt='^Raw EFG matrix';
    shield_sect='CHEMICAL SHIELDINGS';
else
    hfc_hdr='^Raw HFC matrix'; hfc_alt='^Total HFC matrix';
    efg_hdr='^Raw EFG matrix'; efg_alt='^Total EFG matrix';
    shield_sect='CHEMICAL SHIFTS';
end

% Read the last Cartesian coordinate table in the log
sect=find(startsWith(banner_titles,'CARTESIAN COORDINATES (ANGSTROEM)'),1,'last');
if isempty(sect)
    error(['no Cartesian coordinate table in the ORCA log; atom-resolved '...
           'parsing is impossible, rerun ORCA without the miniprint option.']);
end
first=banner_idx(sect); last=sect_end(banner_idx,nlines,first);
symbols={}; std_geom=zeros(0,3); n=first+2;
while (n<=last)&&isempty(log_lines{n}), n=n+1; end
while n<=last
    atom_line=regexp(log_lines{n},'^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)$','tokens','once');
    if isempty(atom_line)||isnan(str2double(atom_line{2})), break; end
    symbols{end+1}=atom_line{1}; %#ok<AGROW>
    std_geom(end+1,:)=[str2double(atom_line{2}) str2double(atom_line{3})...
                       str2double(atom_line{4})]; %#ok<AGROW>
    n=n+1;
end
props.symbols=symbols; props.std_geom=std_geom;
props.natoms=numel(symbols); natoms=props.natoms;
disp('ORCA import: found atomic coordinates.');

% ORCA decorates ghost, dummy, and capped ECP centres with punctuation
% that the per-nucleus property headers do not repeat
bare_symbols=regexprep(symbols,'[^A-Za-z]','');

% Assign atomic numbers, leaving ghost and dummy centres at zero
periodic_table={'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar',...
                'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',...
                'Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',...
                'I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',...
                'Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',...
                'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',...
                'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg'};
props.atomic_numbers=zeros(1,natoms);
for n=1:natoms
    element=find(strcmpi(periodic_table,symbols{n}),1);
    if ~isempty(element), props.atomic_numbers(n)=element; end
end

% Read the charge and the spin multiplicity from the calculation
% settings block; later parts of the log reprint both quantities in
% formats that do not always refer to the system as a whole
charge_idx=find_lines(log_lines,'^Total Charge\s+Charge\s+\.+\s+-?\d+$',1,nlines);
if ~isempty(charge_idx)
    props.charge=str2double(regexp(log_lines{charge_idx(end)},'-?\d+$','match','once'));
end
mult_idx=find_lines(log_lines,'^Multiplicity\s+Mult\s+\.+\s+-?\d+$',1,nlines);
if ~isempty(mult_idx)
    props.multiplicity=str2double(regexp(log_lines{mult_idx(end)},'-?\d+$','match','once'));
end

% Read the final single point energy
energy_idx=find_lines(log_lines,'^FINAL SINGLE POINT ENERGY',1,nlines);
if ~isempty(energy_idx)
    props.energy=str2double(regexp(log_lines{energy_idx(end)},...
                                   '[-+]?\d+\.\d+([eEdD][-+]?\d+)?','match','once'));
    disp('ORCA import: found the total energy.');
end

% Read the electric dipole moment
dipole_idx=find_lines(log_lines,'^Total Dipole Moment',1,nlines);
if ~isempty(dipole_idx)
    numbers=str2double(regexp(log_lines{dipole_idx(end)},...
                              '[-+]?\d+\.\d+([eEdD][-+]?\d+)?','match'));
    if numel(numbers)>=3
        props.dip_moment=numbers((end-2):end);
        disp('ORCA import: found electric dipole moment.');
    end
end

% Read the last Mulliken population analysis
sect=find(startsWith(banner_titles,'MULLIKEN ATOMIC CHARGES'),1,'last');
if ~isempty(sect)
    first=banner_idx(sect); last=sect_end(banner_idx,nlines,first);
    props.mulliken_chg=NaN(natoms,1); spin_pop=NaN(natoms,1);
    for n=(first+2):last
        atom_line=regexp(log_lines{n},'^(\d+)\s+\S+\s*:\s*(\S+)\s*(\S*)','tokens','once');
        if isempty(atom_line), continue; end
        atom=str2double(atom_line{1})+1;
        if atom>natoms, continue; end
        props.mulliken_chg(atom)=str2double(atom_line{2});
        spin_pop(atom)=str2double(atom_line{3});
    end
    disp('ORCA import: found Mulliken population analysis.');

    % Spin populations are only printed for open-shell systems
    if ~all(isnan(spin_pop)), props.mulliken_spin=spin_pop; end
end

% Read the last electronic g-matrix, avoiding the individual
% contribution sections that ORCA prints after the total
sect=find(ismember(banner_titles,{'ELECTRONIC G-MATRIX',...
                                  'ELECTRONIC G-MATRIX FROM EFFECTIVE HAMILTONIAN'}),1,'last');
if ~isempty(sect)
    first=banner_idx(sect); last=sect_end(banner_idx,nlines,first);
    matrix_idx=find_lines(log_lines,'^(The\s+)?g-matrix:',first,last);
    if ~isempty(matrix_idx)
        [g_matrix,found]=read_tensor(log_lines,matrix_idx(1)+1,last);
        if found
            props.g_tensor.raw=g_matrix;
            props.g_tensor.matrix=(g_matrix+g_matrix')/2;
            [eigvecs,eigvals]=eig(props.g_tensor.matrix);
            props.g_tensor.eigvals=diag(eigvals);
            props.g_tensor.eigvecs=eigvecs;
            disp('ORCA import: found the g-tensor.');
        end
    end
end

% Read the last zero-field splitting tensor
sect=find(startsWith(banner_titles,'ZERO-FIELD-SPLITTING TENSOR'),1,'last');
if ~isempty(sect)
    first=banner_idx(sect); last=sect_end(banner_idx,nlines,first);
    matrix_idx=find_lines(log_lines,'^raw-matrix',first,last);
    if ~isempty(matrix_idx)
        [zfs_matrix,found]=read_tensor(log_lines,matrix_idx(1)+1,last);
        if found
            props.zfs.matrix=zfs_matrix;
            [eigvecs,eigvals]=eig((zfs_matrix+zfs_matrix')/2);
            props.zfs.eigvals=diag(eigvals);
            props.zfs.eigvecs=eigvecs;
            disp('ORCA import: found the zero-field splitting tensor.');
        end
    end
end

% Read the last hyperfine and electric field gradient section
sect=find(startsWith(banner_titles,'ELECTRIC AND MAGNETIC HYPERFINE STRUCTURE'),1,'last');
if ~isempty(sect)
    first=banner_idx(sect); last=sect_end(prop_idx,nlines,first);
    nucleus_idx=find_lines(log_lines,'^Nucleus\s+\d+',first,last);
    hfc_found=false; efg_found=false;
    hfc_matrix=cell(1,natoms); hfc_eigvals=cell(1,natoms);
    hfc_eigvecs=cell(1,natoms); hfc_iso=zeros(natoms,1);
    efg_matrix=cell(1,natoms); nqi_matrix=cell(1,natoms);
    isotopes=cell(1,natoms);
    for n=1:numel(nucleus_idx)

        % Each nucleus owns the lines up to the next nucleus header
        block_top=nucleus_idx(n);
        if n<numel(nucleus_idx), block_bot=nucleus_idx(n+1)-1; else, block_bot=last; end

        % ORCA labels the nucleus by its zero-based position in the
        % coordinate table, and only prints the requested ones
        header=regexp(log_lines{block_top},'^Nucleus\s+(\d+)\s*([A-Za-z]{1,2})','tokens','once');
        if isempty(header), continue; end
        atom=str2double(header{1})+1;
        if (atom<1)||(atom>natoms)
            error(['ORCA import: nucleus index ' header{1} ' is outside the coordinate table.']);
        end
        if ~strcmpi(header{2},bare_symbols{atom})
            error(['ORCA import: element mismatch for atom ' num2str(atom) '.']);
        end

        % The isotope, the spin, and the quadrupole moment are printed
        % in the two header lines of the nucleus block
        istp=regexp(log_lines{block_top},'(ISTP|Isotope)=\s*(\d+)','tokens','once');
        if ~isempty(istp), isotopes{atom}=[istp{2} symbols{atom}]; end
        spin_qm=NaN; quad_mom=NaN;
        if block_top<block_bot
            quad_line=regexp(log_lines{block_top+1},...
                             'I=\s*(\S+)\s+Q=\s*(\S+)\s+barn','tokens','once');
            if ~isempty(quad_line)
                spin_qm=str2double(quad_line{1}); quad_mom=str2double(quad_line{2});
            end
        end

        % Read the hyperfine coupling tensor and convert it to Gauss
        header_idx=find_lines(log_lines,hfc_hdr,block_top,block_bot);
        if isempty(header_idx)
            header_idx=find_lines(log_lines,hfc_alt,block_top,block_bot);
        end
        if ~isempty(header_idx)
            [tensor,found]=read_tensor(log_lines,header_idx(1)+1,block_bot);
            if found
                tensor=mhz2gauss(tensor); hfc_matrix{atom}=tensor;
                [eigvecs,eigvals]=eig(tensor);
                hfc_eigvals{atom}=diag(eigvals); hfc_eigvecs{atom}=eigvecs;
                hfc_iso(atom)=trace(tensor)/3; hfc_found=true;
            end
        end

        % Read the electric field gradient tensor
        header_idx=find_lines(log_lines,efg_hdr,block_top,block_bot);
        if isempty(header_idx)
            header_idx=find_lines(log_lines,efg_alt,block_top,block_bot);
        end
        if ~isempty(header_idx)
            [tensor,found]=read_tensor(log_lines,header_idx(1)+1,block_bot);
            if found
                efg_matrix{atom}=tensor; efg_found=true;

                % Quadrupolar tensors follow the gparse convention, in
                % which the spin-dependent denominator is not applied
                if (spin_qm>=1)&&(~isnan(quad_mom))
                    nqi=castep2nqi(tensor,quad_mom,spin_qm)*2*spin_qm*(2*spin_qm-1);
                    nqi_matrix{atom}=nqi-eye(3)*trace(nqi)/3;
                end
            end
        end

    end
    if hfc_found
        props.hfc.full.matrix=hfc_matrix; props.hfc.full.eigvals=hfc_eigvals;
        props.hfc.full.eigvecs=hfc_eigvecs; props.hfc.iso=hfc_iso;
        disp('ORCA import: found hyperfine couplings.');
    end
    if efg_found
        props.efg=efg_matrix;
        disp('ORCA import: found electric field gradients.');
        if any(~cellfun(@isempty,nqi_matrix))
            props.nqi=nqi_matrix;
            disp('ORCA import: found nuclear quadrupole tensors.');
        end
    end
    if any(~cellfun(@isempty,isotopes)), props.isotopes=isotopes; end
end

% Read the last nuclear shielding section
sect=find(startsWith(banner_titles,shield_sect),1,'last');
if ~isempty(sect)
    first=banner_idx(sect); last=sect_end(prop_idx,nlines,first);
    nucleus_idx=find_lines(log_lines,'^Nucleus\s+\d+',first,last);
    cst_matrix=cell(1,natoms); cst_found=false;
    for n=1:numel(nucleus_idx)

        % Each nucleus owns the lines up to the next nucleus header
        block_top=nucleus_idx(n);
        if n<numel(nucleus_idx), block_bot=nucleus_idx(n+1)-1; else, block_bot=last; end

        % ORCA labels the nucleus by its zero-based position in the
        % coordinate table, and only prints the requested ones
        header=regexp(log_lines{block_top},'^Nucleus\s+(\d+)\s*([A-Za-z]{1,2})','tokens','once');
        if isempty(header), continue; end
        atom=str2double(header{1})+1;
        if (atom<1)||(atom>natoms)
            error(['ORCA import: nucleus index ' header{1} ' is outside the coordinate table.']);
        end
        if ~strcmpi(header{2},bare_symbols{atom})
            error(['ORCA import: element mismatch for atom ' num2str(atom) '.']);
        end

        % Read the total shielding tensor
        header_idx=find_lines(log_lines,'^Total shielding tensor',block_top,block_bot);
        if ~isempty(header_idx)
            [tensor,found]=read_tensor(log_lines,header_idx(1)+1,block_bot);
            if found, cst_matrix{atom}=tensor; cst_found=true; end
        end

    end
    if cst_found
        props.cst=cst_matrix;
        disp('ORCA import: found nuclear shielding tensors.');
    end
end

% Read the last spin-spin coupling section, the banner of which
% ORCA centres on the page rather than ruling from the first column
sect_idx=find_lines(log_lines,'^NMR SPIN-SPIN COUPLING CONSTANTS',1,nlines);
if ~isempty(sect_idx)
    first=sect_idx(end); last=sect_end(prop_idx,nlines,first);
    pair_idx=find_lines(log_lines,'^NUCLEUS A =',first,last);
    j_couplings=zeros(natoms,natoms); j_found=false;
    for n=1:numel(pair_idx)

        % Each nuclear pair owns the lines up to the next pair header
        block_top=pair_idx(n);
        if n<numel(pair_idx), block_bot=pair_idx(n+1)-1; else, block_bot=last; end

        % Both nuclei are labelled by their zero-based coordinate index
        header=regexp(log_lines{block_top},...
                      '^NUCLEUS A =\s*\S+\s+(\d+)\s+NUCLEUS B =\s*\S+\s+(\d+)','tokens','once');
        if isempty(header), continue; end
        atom_a=str2double(header{1})+1; atom_b=str2double(header{2})+1;
        if (atom_a>natoms)||(atom_b>natoms)
            error('ORCA import: coupled nucleus index is outside the coordinate table.');
        end

        % The isotropic coupling is one third of the tensor trace
        header_idx=find_lines(log_lines,'^Total spin-spin coupling tensor',block_top,block_bot);
        if ~isempty(header_idx)
            [tensor,found]=read_tensor(log_lines,header_idx(1)+1,block_bot);
            if found
                j_couplings(atom_a,atom_b)=trace(tensor)/3;
                j_couplings(atom_b,atom_a)=trace(tensor)/3; j_found=true;
            end
        end

    end
    if j_found
        props.j_couplings=j_couplings;
        disp('ORCA import: found J-couplings.');
    end
end

% Read the last molar magnetic susceptibility section, the banner
% of which ORCA also centres on the page
sect_idx=find_lines(log_lines,'^TEMPERATURE DEPENDENT MOLAR MAGNETIC SUSCEPTIBILITY',1,nlines);
if ~isempty(sect_idx)
    first=sect_idx(end); last=sect_end(prop_idx,nlines,first);
    temp_idx=find_lines(log_lines,'^TEMPERATURE/K:',first,last);
    chi_temps=zeros(numel(temp_idx),1); chi_tensors=cell(1,numel(temp_idx));
    for n=1:numel(temp_idx)
        chi_temps(n)=str2double(regexp(log_lines{temp_idx(n)},...
                                       '[-+]?\d+\.\d+([eEdD][-+]?\d+)?','match','once'));
        [chi_tensors{n},found]=read_tensor(log_lines,temp_idx(n)+1,last);
        if ~found
            error('ORCA import: malformed magnetic susceptibility tensor.');
        end
    end
    if ~isempty(temp_idx)
        props.chi_temps=chi_temps; props.chi_tensors=chi_tensors;
        disp('ORCA import: found magnetic susceptibility tensors.');
    end
end

end

% Returns the indices of the lines within the specified range that
% match the supplied regular expression
function idx=find_lines(log_lines,pattern,first,last)
if first>last, idx=[]; return; end
hits=regexp(log_lines(first:last),pattern,'once');
idx=first-1+find(~cellfun(@isempty,hits));
end

% Returns the last line of the section that starts at the given line,
% that is the line just above the next banner in the supplied index
function last=sect_end(banner_idx,nlines,first)
next_banner=banner_idx(banner_idx>first);
if isempty(next_banner), last=nlines; else, last=next_banner(1)-1; end
end

% Returns the first 3x3 matrix printed within the specified range of
% lines, located as three consecutive lines that contain nothing but
% three numbers each; ORCA inserts rulers, blank lines, and comments
% between a block header and its matrix in a version-dependent way
function [tensor,found]=read_tensor(log_lines,first,last)
tensor=zeros(3,3); found=false; nrows=0;
for n=first:min(last,first+15)
    [numbers,count,~,next_char]=sscanf(log_lines{n},'%f');
    if (count==3)&&(next_char>numel(log_lines{n}))
        nrows=nrows+1; tensor(nrows,:)=numbers';
        if nrows==3, found=true; return; end
    elseif nrows>0
        nrows=0;
    end
end
tensor=zeros(3,3);
end

% Consistency enforcement
function grumble(file_name)
if ~ischar(file_name)
    error('file_name must be a character string.');
end
end

% I can't think that it would be terrible of me to say - and it
% is occasionally true - that I need physics more than friends.
%
% J. Robert Oppenheimer

