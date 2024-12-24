% A parser for ORCA logs. Syntax:
%
%                    props=oparse(filename)
%
% The following output fields are returned if the corres-
% ponding information is found in the file:
%             
%   props.coords           - atomic coordinates (Angstrom)
%   props.symbols          - atomic symbols
%   props.natoms           - number of atoms
%   props.filename         - log file name
%   props.hfc.full.eigvals - eigenvalues of hyperfine tensors, Gauss
%   props.hfc.full.eigvecs - eigenvectors of hyperfine tensors
%   props.hfc.full.matrix  - hyperfine tensor matrices, Gauss
%   props.efg              - EFG tensors, a.u.^-3
%   props.g_tensor.matrix  - g-tensor matrix, Bohr magneton units
%
% Note: make sure that magnetic parameters are calculated and printed
%       in the log for ALL atoms in the molecule, otherwise the parser
%       would get confused and throw an error.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=oparse.m>

function props=oparse(file_name)

% Check consistency
grumble(file_name);

% Read the file
file_id=fopen(file_name,'r');
orca_log=textscan(file_id,'%s','delimiter','\n');
fclose(file_id); orca_log=orca_log{1};
props.filename=file_name;

% Deblank all lines
for n=1:numel(orca_log), orca_log(n)=deblank(orca_log(n)); end

% Parse the file
for n=1:length(orca_log)
   
   % Read the input orientation and atomic symbols
   if strcmp(orca_log{n},'CARTESIAN COORDINATES (ANGSTROEM)')
      k=n+2; m=1;
      while ~isempty(orca_log{k})
         S=textscan(orca_log{k},'%s %f %f %f'); 
         props.symbols{m}=S{1}{1};
         props.std_geom(m,:)=[S{2:4}];
         k=k+1; m=m+1;
      end
      props.natoms=m-1;
      disp('ORCA import: found atomic coordinates.');
   end
   
   % Read the g-tensor
   if strcmp(orca_log{n},'The g-matrix:')
      props.g_tensor.matrix=zeros(3,3);
      S=textscan(orca_log{n+1},'%f %f %f'); props.g_tensor.matrix(1,:)=cell2mat(S);
      S=textscan(orca_log{n+2},'%f %f %f'); props.g_tensor.matrix(2,:)=cell2mat(S);
      S=textscan(orca_log{n+3},'%f %f %f'); props.g_tensor.matrix(3,:)=cell2mat(S);
      props.g_tensor.matrix=(props.g_tensor.matrix+...
                             props.g_tensor.matrix')/2;
      disp('ORCA import: found the g-tensor.');
   end
   
   % Read HFCs and EFGs
   if strfind(orca_log{n},'Nucleus ')==1
       atom=textscan(orca_log{n},'%*s %f%*s');
       idx=atom{:}(1)+1;
       k=n+1;
       while isempty(strfind(orca_log{k},'Nucleus '))&&k<length(orca_log)
           k=k+1;
           if strcmp(orca_log{k},'Raw HFC matrix (all values in MHz):')==1
               S=textscan(orca_log{k+2},'%f %f %f'); props.hfc.full.matrix{idx}(1,:)=[S{1} S{2} S{3}]; 
               S=textscan(orca_log{k+3},'%f %f %f'); props.hfc.full.matrix{idx}(2,:)=[S{1} S{2} S{3}];
               S=textscan(orca_log{k+4},'%f %f %f'); props.hfc.full.matrix{idx}(3,:)=[S{1} S{2} S{3}];
               props.hfc.full.matrix{idx}=mhz2gauss(props.hfc.full.matrix{idx});
               [A,B]=eig(props.hfc.full.matrix{idx});
               props.hfc.full.eigvals{idx}=diag(B);
               props.hfc.full.eigvecs{idx}=A;               
           end
           if strcmp(orca_log{k},'Raw EFG matrix (all values in a.u.**-3):')==1
               S=textscan(orca_log{k+1},'%f %f %f'); props.efg{idx}(1,:)=[S{1} S{2} S{3}];
               S=textscan(orca_log{k+2},'%f %f %f'); props.efg{idx}(2,:)=[S{1} S{2} S{3}];
               S=textscan(orca_log{k+3},'%f %f %f'); props.efg{idx}(3,:)=[S{1} S{2} S{3}];
           end
       end
   end
   
end

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

