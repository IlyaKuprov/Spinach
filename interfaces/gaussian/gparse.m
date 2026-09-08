% A parser for Gaussian (03, 09, 16) calculation logs. Ex-
% tracts all potentially useful information. Syntax:
%
%             props=gparse(filename,options)
%
% Parameters:
%
%    file_name - a character strong with a file name
%
%    options   - symmetrisation of the interaction
%                tensors. By default all tensors are
%                symmetrised. The symmetrisation may
%                be turned off by adding the following
%                strings to the options cell array:
%
%                      'g_nosymm', 'cst_nosymm',
%                            'hfc_nosymm'
%
% Outputs:
%
%   props.inp_geom         - input geometry (Angstrom)
%   props.std_geom         - standard geometry (Angstrom)
%   props.natoms           - number of atoms
%   props.method           - energy method
%   props.energy           - SCF energy (Hartree)
%   props.hfc.iso          - isotropic hyperfines (Gauss)
%   props.hfc.full.eigvals - HFC eigenvalues (Gauss)
%   props.hfc.full.eigvecs - HFC eigenvectors
%   props.hfc.full.matrix  - HFC tensors (Gauss)
%   props.g_tensor.eigvecs - g-tensor eigenvectors
%   props.g_tensor.eigvals - g-tensor eigenvalues
%   props.g_tensor.matrix  - g-tensor
%   props.cst              - absolute shielding tensors
%   props.k_couplings      - isotropic K-couplings (Hz)
%   props.j_couplings      - isotropic J-couplings (Hz)
%   props.srt              - spin-rotation tensors (Hz)
%   props.nqi              - nuclear quadrupolar tensors (Hz)
%   props.chi              - susceptibility tensor
%   props.gibbs            - Gibbs free energy (Hartree)
%   props.symbols          - atomic symbols; 'Bq' for ghost
%                            atoms, 'X' for dummy atoms, 'TV'
%                            for translation vectors
%   props.isotopes         - mass numbers of the isotopes for
%                            which hyperfine couplings are
%                            reported
%   props.mass_numbers     - mass numbers of the isotopes used
%                            by Gaussian for atomic masses
%   props.atom_masses      - atomic masses used by Gaussian
%                            (atomic mass units)
%   props.nuc_spins        - spin quantum numbers of the mass
%                            isotopes
%   props.nuc_qmom         - quadrupole moments of the mass
%                            isotopes (fm^2)
%   props.nuc_mmom         - magnetic moments of the mass
%                            isotopes (nuclear magnetons)
%   props.atomic_numbers   - atomic numbers
%   props.charge           - overall charge
%   props.el_dip_std       - electric dipole moment, Debye
%   props.multiplicity     - overall multiplicity
%   props.filename         - log file name
%   props.error            - true if the calculation
%                            contains an error of any type
%
% Notes: the following keywords must be added to the route
%        section of the Gaussian input file to produce a
%        useful log:
%
%          #p nmr=(giao,spinspin,susceptibility)
%             output=pickett pop=minimal IOp(6/82=1)
%
%        Gaussian divides its isotropic Fermi contact couplings
%        by 2S=multiplicity-1, but prints the anisotropic spin
%        dipole couplings without that normalisation; the two
%        blocks therefore disagree by 2S for anything above a
%        doublet. This is corrected here, and the hyperfine
%        tensors returned are the ones that enter the spin
%        Hamiltonian as S*A*I, in agreement with oparse.m
%
%        Gaussian prints the spin-rotation and quadrupole ten-
%        sors of the output=pickett block in the principal axis
%        frame of the inertia tensor, not in the standard orien-
%        tation used for everything else. They are rotated here
%        into the standard orientation using a rotation fitted
%        between the principal axis coordinates that Gaussian
%        prints before them and the current orientation. The
%        rotation matrix Gaussian prints is not used because it
%        is an identity matrix in some logs where the principal
%        axis coordinates are visibly permuted.
%
%        When the calculation is run with NoSymm, Gaussian does
%        not reorient the molecule and prints no standard orien-
%        tation; the input orientation is then returned as the
%        standard one because all tensors refer to it.
%
%        In multi-job (Link1) logs, the last occurrence of each
%        quantity is returned. A warning is printed when the
%        number of atoms changes between jobs.
%
%        Parsed g-tensors, shielding tensors, and hyperfine
%        couplings are checked against the principal g-shifts,
%        isotropic shieldings, and unit conversions that Gaus-
%        sian prints alongside them; a disagreement is an error.
%
% gareth.charnock@oerc.ox.ac.uk
% jennifer.handsel@stx.ox.ac.uk
% janm@umbc.edu
% luke.ward@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=gparse.m>

function props=gparse(filename,options)

% Check consistency
if ~exist('options','var'), options={}; end
grumble(filename,options);

% Read the file
file_id=fopen(filename,'r');
g03_output=textscan(file_id,'%s','delimiter','\n');
fclose(file_id); g03_output=g03_output{1};

% Deblank all lines
for n=1:numel(g03_output), g03_output(n)=deblank(g03_output(n)); end

% Locate the route echo lines that start each job and warn about jobs without #p option
route_lines=false(size(g03_output));
for n=2:numel(g03_output)
    current_line=char(g03_output(n)); previous_line=char(g03_output(n-1));
    route_lines(n)=(~isempty(current_line))&&(current_line(1)=='#')&&...
                   (~isempty(previous_line))&&all(previous_line=='-');
    if route_lines(n)&&(~strncmpi(current_line,'#p',2))
        warning('Gaussian import: detailed printing (#p) is off, some blocks may be missing from this log.');
    end
end

% Set the default error flag
props.error=0;

% Set the default complete flag
props.complete=0;

% Parse the file
atoms=[]; atomic_numbers=[]; charge_read=false; new_job=false;
for n=1:length(g03_output)

   % Note the start of each job
   current_line=char(g03_output(n));
   if route_lines(n), new_job=true; charge_read=false; end

   % Read the charge and multiplicity of the whole system, printed first in each job
   if strncmp(current_line,'Charge =',8)&&(~charge_read)
       tokens=regexp(current_line,'Charge =\s*(-?\d+)\s+Multiplicity =\s*(-?\d+)','tokens','once');
       props.charge=str2double(tokens{1});
       props.multiplicity=str2double(tokens{2}); charge_read=true;
       disp('Gaussian import: found charge and multiplicity.');
   end

   % Read the input or standard orientation and atomic numbers
   if strcmp(current_line,'Input orientation:')||strcmp(current_line,'Z-Matrix orientation:')||...
      strcmp(current_line,'Standard orientation:')
      if new_job, atoms=[]; atomic_numbers=[]; new_job=false; end
      k=n+5; m=0; S=sscanf(char(g03_output(k)),'%f');
      while numel(S)>=4
         m=m+1; atoms(m,:)=S((end-2):end)'; atomic_numbers(m)=S(2); %#ok<AGROW>
         k=k+1; S=sscanf(char(g03_output(k)),'%f');
      end
      if isfield(props,'natoms')&&(props.natoms~=m)
          warning('Gaussian import: atom count changed between jobs, tensors may refer to different molecules.');
      end
      natoms=m; props.natoms=natoms;
      if ~strcmp(current_line,'Standard orientation:')
          props.inp_geom=atoms; disp('Gaussian import: found input orientation.');
      else
          props.std_geom=atoms; disp('Gaussian import: found standard orientation.');
      end
   end

   % Read the SCF energy
   if length(current_line)>10 && strcmp(current_line(1:9),'SCF Done:')
       scan_data=textscan(current_line,'SCF Done: %s = %f','Delimiter',' ','MultipleDelimsAsOne',1);
       props.method=char(scan_data{1}); props.method=props.method(3:(end-1)); props.energy=scan_data{2};
       disp('Gaussian import: found SCF energy.');
   end

   % Read the total electron spin
   if length(current_line)>4 && strcmp(current_line(1:4),'<Sx>')
       scan_data=textscan(current_line,'<Sx>= %f <Sy>= %f <Sz>= %f <S**2>= %f S= %f',...
                                       'Delimiter',' ','MultipleDelimsAsOne',1);
       props.spin=[scan_data{1} scan_data{2} scan_data{3}]; props.s_sq=scan_data{4};
       disp('Gaussian import: found total electron spin.');
   end

   % Read the isotope table
   if strcmp(current_line,'Isotopes and Nuclear Properties:')
       [props.mass_numbers,props.atom_masses,props.nuc_spins,props.nuc_qmom,props.nuc_mmom]=deal(zeros(0,1));
       k=n+1; while ~strncmp(char(g03_output(k)),'Atom',4), k=k+1; end
       while true
           table_line=char(g03_output(k)); k=k+1;
           if isempty(table_line), continue; end
           table_values=sscanf(regexprep(table_line,'^[A-Za-z]+=?',''),'%f');
           if (~strncmp(table_line,'Atom',4))&&(numel(table_values)~=numel(atom_index))
               table_values=nan(size(atom_index));
           end
           switch regexp(table_line,'^[A-Za-z]+','match','once')
               case 'Atom',   atom_index=table_values;
               case 'IAtWgt', props.mass_numbers(atom_index,1)=table_values;
               case 'AtmWgt', props.atom_masses(atom_index,1)=table_values;
               case 'NucSpn', props.nuc_spins(atom_index,1)=table_values/2;
               case 'NQMom',  props.nuc_qmom(atom_index,1)=table_values;
               case 'NMagM',  props.nuc_mmom(atom_index,1)=table_values;
               case {'AtZEff','AtZNuc'}
               otherwise, break;
           end
       end
       disp('Gaussian import: found isotope table.');
   end

   % Read isotropic hyperfine couplings and check the unit columns
   if strcmp(current_line,'Isotropic Fermi Contact Couplings')
      props.hfc.iso=zeros(natoms,1); props.isotopes=zeros(natoms,1); k=n+2;
      while true
         tokens=regexp(char(g03_output(k)),'^(\d+)\s+[A-Za-z]+\((\d+)\)?\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)$','tokens');
         if isempty(tokens), break; end
         tokens=str2double(tokens{1}); k=k+1;
         if abs(tokens(4)-2.802495*tokens(5))>(1e-3*abs(tokens(4))+1e-3)
             error('Gaussian import: hyperfine coupling unit columns disagree.');
         end
         props.isotopes(tokens(1),1)=tokens(2); props.hfc.iso(tokens(1),1)=tokens(5);
      end
      disp('Gaussian import: found isotropic hyperfine couplings.');
      disp('Gaussian import: found isotope numbers.');
   end

   % Read, renormalise, and symmetrize anisotropic hyperfine couplings
   if strcmp(current_line,'Anisotropic Spin Dipole Couplings in Principal Axis System')&&(~isfield(props,'multiplicity'))
      error('Gaussian import: hyperfine normalisation needs the multiplicity, which this log does not print.');
   elseif strcmp(current_line,'Anisotropic Spin Dipole Couplings in Principal Axis System')&&(props.multiplicity==1)
      if isfield(props,'hfc'), props=rmfield(props,'hfc'); end
      warning('Gaussian import: hyperfine couplings of a closed-shell system ignored.');
   elseif strcmp(current_line,'Anisotropic Spin Dipole Couplings in Principal Axis System')
      props.hfc.full.eigvals=cell(natoms,1);
      props.hfc.full.eigvecs=cell(natoms,1);
      props.hfc.full.matrix=cell(natoms,1);
      k=n+1; while ~strncmp(char(g03_output(k)),'Baa',3), k=k+1; end
      while strncmp(char(g03_output(k)),'Baa',3)
         baa=sscanf(regexprep(char(g03_output(k)),'^.*B[abc][abc]',''),'%f')';
         bbb=sscanf(regexprep(char(g03_output(k+1)),'^.*B[abc][abc]',''),'%f')';
         bcc=sscanf(regexprep(char(g03_output(k+2)),'^.*B[abc][abc]',''),'%f')';
         atom_num=sscanf(char(g03_output(k+1)),'%d');
         if any(abs([baa(2) bbb(2) bcc(2)]-2.802495*[baa(3) bbb(3) bcc(3)])>(1e-2*abs([baa(2) bbb(2) bcc(2)])+1e-2))
             error('Gaussian import: hyperfine coupling unit columns disagree.');
         end

         % Renormalise the dipolar part and add the isotropic part in the
         % laboratory basis, where it stays exactly isotropic even though
         % Gaussian's four-decimal eigenvectors are not exactly orthonormal
         dip_part=[baa(3) bbb(3) bcc(3)]/(props.multiplicity-1);
         props.hfc.full.eigvals{atom_num}=dip_part+props.hfc.iso(atom_num);
         props.hfc.full.eigvecs{atom_num}=[baa(5:7)' bbb(5:7)' bcc(5:7)'];
         props.hfc.full.matrix{atom_num}=props.hfc.full.eigvecs{atom_num}*diag(dip_part)*...
                                         props.hfc.full.eigvecs{atom_num}'+props.hfc.iso(atom_num)*eye(3);
         if ~ismember('hfc_nosymm',options)
            props.hfc.full.matrix{atom_num}=(props.hfc.full.matrix{atom_num}+props.hfc.full.matrix{atom_num}')/2;
         end
         k=k+3; while isempty(char(g03_output(k))), k=k+1; end
      end
      disp('Gaussian import: found anisotropic hyperfine couplings.');
   end

   % Read the g-tensor and check it against the printed g-shifts
   if strcmp(current_line,'g tensor [g = g_e + g_RMC + g_DC + g_OZ/SOC]:')||strcmp(current_line,'g tensor (ppm):')
      disp('Gaussian import: found g-tensor.'); g=zeros(3);
      for k=1:3
          g(k,:)=sscanf(regexprep(strrep(char(g03_output(n+k)),'D','E'),'[A-Za-z]+=',' '),'%f')';
      end
      if strcmp(char(g03_output(n+4)),'g shifts relative to the free electron (ppm):')
          g_shifts=sscanf(regexprep(char(g03_output(n+5)),'[a-z]+=',' '),'%f')';
          if max(abs(sort(eig((g+g')/2))'-sort(2.0023193043+1e-6*g_shifts)))>1e-6
              error('Gaussian import: parsed g-tensor disagrees with the printed g-shifts.');
          end
      end
      if ~ismember('g_nosymm',options)
         g=(g+g')/2; disp('Gaussian import: g-tensor symmetrized.');
      else
         disp('Gaussian import: g-tensor not symmetrized.');
      end
      [V,D]=eig(g);
      props.g_tensor.eigvecs=V;
      props.g_tensor.eigvals=diag(D)';
      props.g_tensor.matrix=g;
   end

   % Read chemical shielding tensors and check them against the printed isotropic values
   if strcmp(current_line,'SCF GIAO Magnetic shielding tensor (ppm):')||...
      strcmp(current_line,'MP2 GIAO Magnetic shielding tensor (ppm):')||...
      strcmp(current_line,'Magnetic shielding (ppm):')
       disp('Gaussian import: found chemical shielding tensors.');
       cst=cell(natoms,1); k=n+1;
       while ~isempty(regexp(char(g03_output(k)),'^\d+\s+\S+\s+Isotropic =','once'))
           atom_num=sscanf(char(g03_output(k)),'%d'); cst{atom_num}=zeros(3);
           for m=1:3
               cst{atom_num}(m,:)=sscanf(regexprep(char(g03_output(k+m)),{'[A-Za-z]+=','\*+'},{' ',' NaN '}),'%f')';
           end
           iso_value=sscanf(regexprep(char(g03_output(k)),{'^.*Isotropic =','\*+'},{'',' NaN '}),'%f');
           if abs(trace(cst{atom_num})/3-iso_value(1))>1e-3
               error('Gaussian import: parsed shielding tensor disagrees with the printed isotropic value.');
           end
           k=k+5+4*strcmp(char(g03_output(k+5)),'Eigenvectors:');
       end
       if ~ismember('cst_nosymm',options)
           disp('Gaussian import: shielding tensors symmetrized.');
           for k=1:natoms, cst{k}=(cst{k}+cst{k}')/2; end
       else
           disp('Gaussian import: shielding tensors not symmetrized.');
       end
       props.cst=cst;
   end

   % Read reduced K-couplings or isotope-specific J-couplings using the printed atom indices
   if strcmp(current_line,'Total nuclear spin-spin coupling K (Hz):')||...
      strcmp(current_line,'Total nuclear spin-spin coupling J (Hz):')
       couplings=zeros(natoms); k=n+1;
       while true
           block_line=char(g03_output(k)); k=k+1;
           block_values=sscanf(strrep(block_line,'D','E'),'%f')';
           if isempty(block_values), break; end
           if contains(block_line,'D')
               couplings(block_values(1),columns(1:(numel(block_values)-1)))=block_values(2:end);
           else
               columns=block_values;
           end
       end
       couplings=couplings+couplings';
       if contains(current_line,'K (Hz)')
           props.k_couplings=couplings; disp('Gaussian import: found isotropic K-couplings.');
       else
           props.j_couplings=couplings; disp('Gaussian import: found isotropic J-couplings.');
       end
   end

   % Fit the rotation from the inertial frame into the current orientation using the printed coordinates
   if strcmp(current_line,'Principal axis orientation:')
       k=n+5; abc_geom=zeros(0,3); S=sscanf(char(g03_output(k)),'%f');
       while numel(S)>=4
           abc_geom(end+1,:)=S((end-2):end)'; k=k+1; S=sscanf(char(g03_output(k)),'%f'); %#ok<AGROW>
       end
       cur_geom=atoms; if size(abc_geom,1)~=size(atoms,1), cur_geom=atoms(atomic_numbers>0,:); end
       cur_geom=cur_geom-mean(cur_geom,1); abc_geom=abc_geom-mean(abc_geom,1);
       [U,~,V]=svd(abc_geom'*cur_geom); abc_rot=V*diag([1 1 det(V*U')])*U';
       if norm(cur_geom-abc_geom*abc_rot','fro')>(1e-3*norm(cur_geom,'fro')+1e-3)
           error('Gaussian import: inertial frame coordinates do not match the molecule.');
       end
       disp('Gaussian import: found the inertial frame orientation.');
   end

   % Read spin-rotation or quadrupole tensors and rotate them into the standard orientation
   if strcmp(current_line,'nuclear spin - molecular rotation tensor [C] (MHz):')||...
      strcmp(current_line,'Nuclear quadrupole coupling constants [Chi] (MHz):')
       tensors=cell(natoms,1); k=n+1;
       while ~isempty(regexp(char(g03_output(k)),'^\d+\s+[A-Za-z]','once'))
           atom_num=sscanf(char(g03_output(k)),'%d'); T=zeros(3);
           for m=1:3
               T(m,:)=sscanf(regexprep(char(g03_output(k+m)),'[a-c]{2}=',' '),'%f')';
           end
           tensors{atom_num}=1e6*abc_rot*T*abc_rot'; k=k+4;
       end
       if contains(current_line,'[C]')
           props.srt=tensors; disp('Gaussian import: found nuclear spin-rotation tensors.');
       else
           for k=1:natoms
               if ~isempty(tensors{k}), tensors{k}=tensors{k}-eye(3)*trace(tensors{k})/3; end
           end
           props.nqi=tensors; disp('Gaussian import: found nuclear quadrupole tensors.');
       end
   end

   % Read the magnetic susceptibility tensor (GIAO or CSGT)
   if strcmp(current_line,'Magnetic susceptibility tensor (cgs-ppm):')||...
      strcmp(current_line,'Magnetic susceptibility (cgs-ppm):')
       k=n+strcmp(current_line,'Magnetic susceptibility (cgs-ppm):'); chi=zeros(3);
       for m=1:3
           chi(m,:)=sscanf(regexprep(char(g03_output(k+m)),'[A-Za-z]+=',' '),'%f')';
       end
       props.chi=cgsppm2ang(chi); props.chi=(props.chi+props.chi')/2;
       disp('Gaussian import: found and symmetrised magnetic susceptibility tensor.');
   end

   % Read Gibbs free energy
   if strncmp(current_line,'Sum of electronic and thermal Free Energies',43)
       props.gibbs=sscanf(current_line((strfind(current_line,'=')+1):end),'%f');
       disp('Gaussian import: found Gibbs free energy.');
   end

   % Read electric dipole moment
   if strcmp(current_line,'Dipole moment (field-independent basis, Debye):')
       dip_values=regexp(char(g03_output(n+1)),'([-+]?[0-9]*\.?[0-9]+)','match');
       dip_values=str2double(dip_values); props.el_dip_std=dip_values(1:3);
       disp('Gaussian import: found electric dipole moment.');
   end

   % Check for error flags
   if strncmp(current_line,'Error termination',17)||strncmp(current_line,'Erroneous write.',16)
       props.error=1; warning('Gaussian import: error message detected.');
   end

   % Check for incomplete calculations
   if strncmp(current_line,'Normal termination',18)
       props.complete=1;
   end
   if strncmp(current_line,'Entering Gaussian System',24)
       props.complete=0;
   end

end

% Warn about incomplete logs
if (~props.error)&&(~props.complete)
   warning('Gaussian import: incomplete log detected.');
end

% Use the input orientation when Gaussian did not reorient the molecule
if (~isfield(props,'std_geom'))&&isfield(props,'inp_geom')
    props.std_geom=props.inp_geom;
    disp('Gaussian import: no standard orientation found, using the input orientation.');
end

% Assign atomic symbols, including translation vectors, dummies, and ghosts
periodic_table={'TV','X','Bq','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar',...
                'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',...
                'Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',...
                'I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',...
                'Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',...
                'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',...
                'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg'};
if ~isempty(atomic_numbers)
    props.symbols=periodic_table(atomic_numbers+3);
    props.atomic_numbers=atomic_numbers;
end

% Assign source information
props.filename=filename;

end

% Consistency enforcement
function grumble(filename,options)
if (~ischar(filename))||isempty(filename)
    error('filename must be a non-empty character string.');
end
if ~exist(filename,'file')
    error('the file specified was not found.');
end
if (~iscell(options))||any(~cellfun(@ischar,options))
    error('options must be a cell array of character strings.');
end
if ~all(ismember(options,{'hfc_nosymm','g_nosymm','cst_nosymm'}))
    error('invalid option specification.');
end
end

% Moral indignation is jealousy with a halo.
%
% H.G. Wells

