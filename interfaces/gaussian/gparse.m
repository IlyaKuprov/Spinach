% A parser for Gaussian (03, 09, 16) calculation logs. Ex-
% tracts all potentially useful information. The following
% options must be added to the route section of the input
% file to produce a useful log:
%
%       #p nmr=(giao,spinspin,susceptibility) 
%          output=pickett pop=minimal IOp(6/82=1)
%
% The options parameter controls the symmetrization of the
% interaction tensors after they have been read in. By de-
% fault all tensors are symmetrized. The symmetrization may
% be turned off by adding the following strings to the op-
% tions cell array:
%
%         'g_nosymm', 'cst_nosymm', 'hfc_nosymm'
%
% The following output fields are returned if the corres-
% ponding information is found in the file:
%             
%   props.inp_geom         - input geometry (Angstrom)
%   props.std_geom         - standard geometry (Angstrom)
%   props.natoms           - number of atoms
%   props.method           - energy method
%   props.energy           - SCF energy (Hartree)
%   props.hfc.iso          - isotropic hyperfines
%   props.hfc.full.eigvals - HFC eigenvalues (Gauss)
%   props.hfc.full.eigvecs - HFC eigenvectors
%   props.hfc.full.matrix  - HFC tensors (Gauss)
%   props.g_tensor.eigvecs - g-tensor eigenvectors
%   props.g_tensor.eigvals - g-tensor eigenvalues
%   props.g_tensor.matrix  - g-tensor
%   props.cst              - absolute shielding tensors
%   props.k_couplings      - isotropic K-couplings (Hz)
%   props.j_couplings      - isotropic J-couplings (Hz)
%   props.srt              - spin-rotation tensor
%   props.nqi              - nuclear quadrupolar tensors
%   props.chi              - susceptibility tensor
%   props.gibbs            - Gibbs free energy
%   props.symbols          - atomic symbols
%   props.isotopes         - nuclear isotopes used by Gaussian
%   props.atomic_numbers   - atomic numbers
%   props.charge           - overall charge
%   props.multiplicity     - overall multiplicity
%   props.filename         - log file name
%   props.error            - true if the calculation
%                            contains an error of any type
%
% gareth.charnock@oerc.ox.ac.uk
% jennifer.handsel@stx.ox.ac.uk
% janm@umbc.edu
% i.kuprov@soton.ac.uk
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

% Refuse to process files without #p option
for n=1:numel(g03_output)
    current_line=char(g03_output(n));
    if (numel(current_line)>=2)&&strcmp(current_line(1:2),'# ')
        error('This Gaussian log has the detailed printing option switched off. Re-run with #p in the route line.');
    end
end

% Set the default error flag
props.error=0;

% Set the default complete flag
props.complete=0;

% Parse the file
for n=1:length(g03_output)
    
   % Read charge and multiplicity
   current_line=char(g03_output(n));
   if length(current_line)>8 && strcmp(current_line(1:8),'Charge =')
       scan_data=textscan(current_line,'%*s %*s %f %*s %*s %f',...
                                       'Delimiter',' ','MultipleDelimsAsOne',1);
       props.charge=scan_data{1}; 
       props.multiplicity=scan_data{2};
       disp('Gaussian import: found charge and multiplicity.');
   end
   
   % Read the input orientation and atomic numbers
   if strcmp(g03_output(n),'Input orientation:')
      k=n+5; m=1;
      while ~strcmp(deblank(g03_output(k)),'---------------------------------------------------------------------')
         S=eval(['[' char(g03_output(k)) ']']); atoms(m,:)=S(4:6); atomic_numbers(m)=S(2); k=k+1; m=m+1; %#ok<AGROW>
      end
      natoms=m-1; props.inp_geom=atoms; props.natoms=natoms;
      disp('Gaussian import: found input orientation.');
   end

   % Read the standard orientation and atomic numbers
   if strcmp(g03_output(n),'Standard orientation:')
      k=n+5; m=1;
      while ~strcmp(deblank(g03_output(k)),'---------------------------------------------------------------------')
         S=eval(['[' char(g03_output(k)) ']']); atoms(m,:)=S(4:6); atomic_numbers(m)=S(2); k=k+1; m=m+1;
      end
      natoms=m-1; props.std_geom=atoms; props.natoms=natoms;
      disp('Gaussian import: found standard orientation.');
   end

   % Read the SCF energy
   current_line=char(g03_output(n));
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

   % Read isotropic hyperfine couplings
   if strcmp(g03_output(n),'Isotropic Fermi Contact Couplings')
      props.hfc.iso=zeros(natoms,1);
      props.isotopes=zeros(natoms,1);
      for k=1:natoms
         S=char(g03_output(n+1+k));
         first_pass=textscan(S,'%*f %s %*f %*f %f %*f','Delimiter',...
                               ' ','MultipleDelimsAsOne',1);
         props.hfc.iso(k)=first_pass{2};
         if ~strcmp(first_pass{1}{1}(end),')')
             first_pass{1}{1}(end+1)=')'; % Gaussian forgets a bracket for Ab(123)
         end
         second_pass=textscan(first_pass{1}{1},'%*[^(]%*[(]%f%*[)]',...
                              'Delimiter',' ','MultipleDelimsAsOne',1);
         props.isotopes(k)=second_pass{1};
      end
      disp('Gaussian import: found isotropic hyperfine couplings.');
      disp('Gaussian import: found found isotope numbers.');
   end

   % Read and symmetrize anisotropic hyperfine couplings
   if strcmp(g03_output(n),'Anisotropic Spin Dipole Couplings in Principal Axis System')
      props.hfc.full.eigvals=cell(natoms,1);
      props.hfc.full.eigvecs=cell(natoms,1);
      props.hfc.full.matrix=cell(natoms,1);
      for k=1:natoms
         baa=g03_output(n+4*k+1); baa=char(baa); baa=eval(['[' baa(4:end) ']']);
         bbb=g03_output(n+4*k+2); bbb=char(bbb); bbb=eval(['[' bbb(15:end) ']']);
         bcc=g03_output(n+4*k+3); bcc=char(bcc); bcc=eval(['[' bcc(4:end) ']']);
         props.hfc.full.eigvals{k}=[baa(3) bbb(3) bcc(3)]+props.hfc.iso(k);
         props.hfc.full.eigvecs{k}=[baa(5:7)' bbb(5:7)' bcc(5:7)'];
         props.hfc.full.matrix{k}=props.hfc.full.eigvecs{k}*diag(props.hfc.full.eigvals{k})*props.hfc.full.eigvecs{k}';
         if (~exist('options','var'))||(~ismember('hfc_nosymm',options))
            props.hfc.full.matrix{k}=(props.hfc.full.matrix{k}+props.hfc.full.matrix{k}')/2;
         end
      end
      disp('Gaussian import: found anisotropic hyperfine couplings.');
   end

   % Read the g-tensor (Gaussian03 for Unix)
   if strcmp(g03_output(n),'g tensor [g = g_e + g_RMC + g_DC + g_OZ/SOC]:')
      disp('Gaussian import: found g-tensor.');
      line1=char(deblank(g03_output(n+1))); 
      line2=char(deblank(g03_output(n+2))); 
      line3=char(deblank(g03_output(n+3)));
      g=eval(['[' line1(5:18) '   ' line1(26:39) '   ' line1(46:60) ';   '...
                  line2(5:18) '   ' line2(26:39) '   ' line2(46:60) ';   '...
                  line3(5:18) '   ' line3(26:39) '   ' line3(46:60) ']']);
      if (~exist('options','var'))||(~ismember('g_nosymm',options))
         g=(g+g')/2; disp('Gaussian import: g-tensor symmetrized.');
      else
         disp('Gaussian import: g-tensor not symmetrized.');
      end
      [V,D]=eig(g);
      props.g_tensor.eigvecs=V;
      props.g_tensor.eigvals=diag(D)';
      props.g_tensor.matrix=g;
   end
   
   % Read the g-tensor (Gaussian03 for Windows)
   if strcmp(g03_output(n),'g tensor (ppm):')
      disp('Gaussian import: found g-tensor.');
      line1=char(deblank(g03_output(n+1))); 
      line2=char(deblank(g03_output(n+2))); 
      line3=char(deblank(g03_output(n+3)));
      g=eval(['[' line1(4:15) '   ' line1(23:35) '   ' line1(42:54) ';   '...
                  line2(4:15) '   ' line2(23:35) '   ' line2(42:54) ';   '...
                  line3(4:15) '   ' line3(23:35) '   ' line3(42:54) ']']);
      if (~exist('options','var'))||(~ismember('g_nosymm',options))
         g=(g+g')/2; disp('g-tensor symmetrized.');
      else
         disp('Gaussian import: g-tensor not symmetrized.');
      end
      [V,D]=eig(g);
      props.g_tensor.eigvecs=V;
      props.g_tensor.eigvals=diag(D)';
      props.g_tensor.matrix=g;
   end
   
   % Read chemical shielding tensors
   if strcmp(g03_output(n),'SCF GIAO Magnetic shielding tensor (ppm):')||...
      strcmp(g03_output(n),'MP2 GIAO Magnetic shielding tensor (ppm):')||...
      strcmp(g03_output(n),'Magnetic shielding (ppm):')
       disp('Gaussian import: found chemical shielding tensors.');
       cst=cell(natoms,1); 
       for k=1:natoms
           line1=char(g03_output(n+5*k-3)); 
           line2=char(g03_output(n+5*k-2)); 
           line3=char(g03_output(n+5*k-1)); cst{k}=zeros(3,3);
           cst{k}(1,:)=cell2mat(textscan(line1,'XX= %f YX= %f ZX= %f','Delimiter',...
                                               ' ','MultipleDelimsAsOne',1));
           cst{k}(2,:)=cell2mat(textscan(line2,'XY= %f YY= %f ZY= %f','Delimiter',...
                                               ' ','MultipleDelimsAsOne',1));
           cst{k}(3,:)=cell2mat(textscan(line3,'XZ= %f YZ= %f ZZ= %f','Delimiter',...
                                               ' ','MultipleDelimsAsOne',1)); 
       end
       if (~exist('options','var'))||(~ismember('cst_nosymm',options))
           disp('Gaussian import: shielding tensors symmetrized.');
           for k=1:natoms, cst{k}=(cst{k}+cst{k}')/2; end
       else
           disp('Gaussian import: shielding tensors not symmetrized.');
       end
       props.cst=cst;
   end
   
   % Read isotope-independent reduced scalar couplings
   if strcmp(g03_output(n),'Total nuclear spin-spin coupling K (Hz):')
       k_couplings=zeros(natoms,ceil(natoms/5)*5);
       for k=1:ceil(natoms/5)
           lines_to_read=natoms-(k-1)*5;
           current_block=zeros(natoms,5);
           for m=(n+2):(n+1+lines_to_read)
               current_line=eval(['[' char(g03_output(m)) ']']);
               current_block(current_line(1),1:(length(current_line)-1))=current_line(2:end);
           end
           k_couplings(:,(5*(k-1)+1):(5*k))=current_block;
           n=n+lines_to_read+1; %#ok<FXSET>
       end
       k_couplings=k_couplings(1:natoms,1:natoms); k_couplings=k_couplings+k_couplings';
       props.k_couplings=k_couplings;
       disp('Gaussian import: found isotropic K-couplings.');
   end
   
   % Read isotope-specific scalar couplings
   if strcmp(g03_output(n),'Total nuclear spin-spin coupling J (Hz):')
       j_couplings=zeros(natoms,ceil(natoms/5)*5);
       for k=1:ceil(natoms/5)
           lines_to_read=natoms-(k-1)*5;
           current_block=zeros(natoms,5);
           for m=(n+2):(n+1+lines_to_read)
               current_line=eval(['[' char(g03_output(m)) ']']);
               current_block(current_line(1),1:(length(current_line)-1))=current_line(2:end);
           end
           j_couplings(:,(5*(k-1)+1):(5*k))=current_block;
           n=n+lines_to_read+1; %#ok<FXSET>
       end
       j_couplings=j_couplings(1:natoms,1:natoms); j_couplings=j_couplings+j_couplings';
       props.j_couplings=j_couplings;
       disp('Gaussian import: found isotropic J-couplings.');
   end
   
   % Read spin-rotation couplings
   if strcmp(g03_output(n),'nuclear spin - molecular rotation tensor [C] (MHz):')
       props.srt=cell(natoms,1);
       for k=1:natoms
           line3=char(g03_output(n+4*k));   line2=char(g03_output(n+4*k-1));
           line1=char(g03_output(n+4*k-2)); line0=char(g03_output(n+4*k-3));
           atom_num=regexp(line0,'^[0-9]*','match'); atom_num=eval(atom_num{1});
           props.srt{atom_num}=[eval(['[' line1(4:14) '     ' line1(21:31) '     ' line1(38:48) ']']);
                                eval(['[' line2(4:14) '     ' line2(21:31) '     ' line2(38:48) ']']);
                                eval(['[' line3(4:14) '     ' line3(21:31) '     ' line3(38:48) ']'])]*1e6;
           if strcmp(g03_output(n+4*k+1),'Dipole moment (Debye):'), break; end
           if strcmp(g03_output(n+4*k+1),'Nuclear quadrupole coupling constants [Chi] (MHz):'), break; end
       end
       disp('Gaussian import: found nuclear spin-rotation tensors.');
   end
   
   % Read quadrupole couplings and kill their trace
   if strcmp(g03_output(n),'Nuclear quadrupole coupling constants [Chi] (MHz):')
       props.nqi=cell(natoms,1);
       for k=1:natoms
           line3=char(g03_output(n+4*k));   line2=char(g03_output(n+4*k-1)); 
           line1=char(g03_output(n+4*k-2)); line0=char(g03_output(n+4*k-3));
           atom_num=regexp(line0,'^[0-9]*','match'); atom_num=eval(atom_num{1});
           props.nqi{atom_num}=[eval(['[' line1(4:14) '     ' line1(21:31) '     ' line1(38:48) ']']);
                                eval(['[' line2(4:14) '     ' line2(21:31) '     ' line2(38:48) ']']);
                                eval(['[' line3(4:14) '     ' line3(21:31) '     ' line3(38:48) ']'])]*1e6;
           props.nqi{atom_num}=props.nqi{atom_num}-eye(3)*trace(props.nqi{atom_num})/3;
           if strcmp(g03_output(n+4*k+1),'Dipole moment (Debye):'), break; end
       end
       disp('Gaussian import: found nuclear quadrupole tensors.');
   end
   
   % Read magnetic susceptibility tensor (GIAO)
   if strcmp(g03_output(n),'Magnetic susceptibility tensor (cgs-ppm):')
       line1=char(g03_output(n+1)); line2=char(g03_output(n+2)); line3=char(g03_output(n+3));
       props.chi=[eval(['[' line1(4:19) '     ' line1(25:40) '     ' line1(46:end) ']']);
                  eval(['[' line2(4:19) '     ' line2(25:40) '     ' line2(46:end) ']']);
                  eval(['[' line3(4:19) '     ' line3(25:40) '     ' line3(46:end) ']'])];
       props.chi=cgsppm2ang(props.chi); props.chi=(props.chi+props.chi')/2;
       disp('Gaussian import: found and symmetrised magnetic susceptibility tensor.');
   end
   
   % Read magnetic susceptibility tensor (CSGT)
   if strcmp(g03_output(n),'Magnetic susceptibility (cgs-ppm):')
       line1=char(g03_output(n+2)); line2=char(g03_output(n+3)); line3=char(g03_output(n+4));
       props.chi=[eval(['[' line1(4:14) '     ' line1(21:31) '     ' line1(38:48) ']']);
                  eval(['[' line2(4:14) '     ' line2(21:31) '     ' line2(38:48) ']']);
                  eval(['[' line3(4:14) '     ' line3(21:31) '     ' line3(38:48) ']'])];
       props.chi=cgsppm2ang(props.chi); props.chi=(props.chi+props.chi')/2;
       disp('Gaussian import: found and symmetrised magnetic susceptibility tensor.');
   end
   
   % Read Gibbs free energy
   if (numel(g03_output{n})>43)&&strcmp(g03_output{n}(1:43),'Sum of electronic and thermal Free Energies')
       props.gibbs=eval(g03_output{n}(45:end));
       disp('Gaussian import: found Gibbs free energy.');
   end
   
   % Read electric dipole moment
   if strcmp(g03_output(n),'Electric dipole moment (input orientation):')
       S=char(g03_output(n+3));
       S=textscan(S,'%s %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1);
       props.electric_dip=S{3};
       disp('Gaussian import: found electric dipole moment.');
   end
   
   % Check for error flags
   if ((numel(g03_output{n})>17)&&strcmp(g03_output{n}(1:17),'Error termination'))||...
      ((numel(g03_output{n})>16)&&strcmp(g03_output{n}(1:16),'Erroneous write.'))
       props.error=1; warning('Gaussian import: error message detected.');
   end
   
   % Check for incomplete calculations
   if (numel(g03_output{n})>18)&&strcmp(g03_output{n}(1:18),'Normal termination')
       props.complete=1;
   end
   if (numel(g03_output{n})>24)&&strcmp(g03_output{n}(1:24),'Entering Gaussian System')
       props.complete=0;
   end
   
end
   
% Assign atomic symbols
periodic_table={'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar',...
                'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',...
                'Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',...
                'I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',...
                'Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',...
                'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',...
                'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg'};

if exist('atomic_numbers','var')
    props.symbols=periodic_table(atomic_numbers);
    props.atomic_numbers=atomic_numbers;
end

% Assign source information
props.filename=filename;

end

% Consistency enforcement
function grumble(filename,options)
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

