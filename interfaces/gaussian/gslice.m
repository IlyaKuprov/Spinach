% Slices a Gaussian geometry scan log into property calculation 
% inputs at the energy minimum geometries. The function asks for
% the log and for a text file containing the property calcula-
% tion header. Some paths are hard-coded; edit as appropriate.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=gslice.m>

function gslice()

% Assign atomic symbols
periodic_table={'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar',...
                'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',...
                'Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',...
                'I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',...
                'Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',...
                'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',...
                'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg'};

% Read the log file
[g03_output,g03_output_path]=uigetfile('*.*','Relaxed geometry scan log');
g03_output=textread([g03_output_path g03_output],'%s','bufsize',65536,'delimiter','\n'); %#ok<DTXTRD>

% Locate and the stationary point reports
minima=[];
for n=1:length(g03_output)
    if strcmp(deblank(g03_output(n)),'-- Stationary point found.')
        minima=[minima n]; %#ok<AGROW>
    end
end
disp(['slice: ' num2str(length(minima)) ' stationary point flags found.']);
               
% Locate standard orientation entry points
std_geom_start=[];
for n=minima
    for k=n:-1:1
        if strcmp(deblank(g03_output(k)),'Standard orientation:')
            std_geom_start=[std_geom_start k+5]; %#ok<AGROW>
            break
        end
    end
end
disp(['slice: ' num2str(length(std_geom_start)) ' standard orientation headers located.']);

% Locate standard orientation end points
std_geom_end=[];
for n=std_geom_start
    for k=n:length(g03_output)
        if strcmp(deblank(g03_output(k)),'---------------------------------------------------------------------')
            std_geom_end=[std_geom_end k-1]; %#ok<AGROW>
            break
        end
    end
end
disp(['slice: ' num2str(length(std_geom_start)) ' standard orientation footers located.']);

% Read the standard orientations
coord_blocks=cell(length(std_geom_start),1);
for n=1:length(std_geom_start)
    coord_blocks{n}=zeros(std_geom_end(n)-std_geom_start(n),6);
    for k=std_geom_start(n):std_geom_end(n)
        coord_blocks{n}(k-std_geom_start(n)+1,:)=cell2mat(textscan(g03_output{k},'%n %n %n %n %n %n'));
    end
end

% Get the header
[g03_header,g03_header_path]=uigetfile('*.*','Gaussian header');
header=textread([g03_header_path g03_header],'%s','delimiter','\n'); %#ok<DTXTRD>

% Write the inputs
runscript=fopen('compute.bat','a');
fprintf(runscript,'%s\n','@ECHO OFF');
fprintf(runscript,'%s\n','set GAUSS_EXEDIR=C:\G16W');
fprintf(runscript,'%s\n','set GAUSS_SCRDIR=C:\Temp');
for n=1:length(coord_blocks)
    g03_input=fopen(['g16_input_' num2str(n) '.gjf'],'a');
    dump(g03_input,header);
    for k=1:size(coord_blocks{n},1)
        current_line=[periodic_table{coord_blocks{n}(k,2)} '       ' num2str(coord_blocks{n}(k,4:6),'%1.8f  ')];
        fprintf(g03_input,'%s\n',current_line);
    end
    fprintf(g03_input,'%s\n','');
    fclose(g03_input);
    fprintf(runscript,'%s\n',['C:\G16W\g16.exe g16_input_' num2str(n) '.gjf']);
end
fclose(runscript);

end

% File writer subroutine
function dump(file_id,cell_array)
for p=1:length(cell_array)
    fprintf(file_id,'%s\n',cell_array{p});
end
end

% Don't show a half-done job to a fool. 
%
% A Russian saying

% #NHEAD #NGRUM