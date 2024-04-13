% Reads GISSMO files and forms Spinach data structures. Syntax:
%
%         [sys,inter]=gissmo2spinach(file_name,subsystem)
%
% Parameters: 
%
%    file_name  -  character string with the name of 
%                  the GISSMO XML file
%
%    subsystem  -  which of the coupling matrices to 
%                  to import
%
% Outputs:
%
%    sys, inter -  Spinach data structures, ready for
%                  calling create.m
%
% Note: GISSMO only provides chemical shifts, J-couplings, the
%       non-selective line width, and the magnet field. You may
%       want to add further parameters by editing sys and inter
%       data structures manually.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=gissmo2spinach.m>

function [sys,inter]=gissmo2spinach(filename,subsystem)

% Check consistency
grumble(filename);

% Parse the XML
disp('parsing GISSMO xml...');
xml=parsexml(filename);

% Set essential flags
magnet=false(); lw=false();
shifts=false(); jc=false();
current_subsystem=1;

% Look inside
for n=1:numel(xml.children)
    
    % Look for magnetic induction
    if strcmpi(xml.children(n).name,'field_strength')
        
        % Get magnetic induction
        sys.magnet=2*pi*1e6*str2double(xml.children(n).children.data)/spin('1H');
        
        % Inform the user
        disp('found magnet induction'); magnet=true();
        
    end
    
    % Look for a coupling matrix
    if strcmpi(xml.children(n).name,'coupling_matrix')&&...
              (current_subsystem==subsystem)
        
        % Grab a shorthand
        cm=xml.children(n);
        
        % Look inside
        for k=1:numel(cm.children)
            
            % Look for line width
            if strcmpi(cm.children(k).name,'lw')
                
                % Get the line width
                inter.relaxation={'damp'};
                inter.rlx_keep='diagonal';
                inter.equilibrium='zero';
                inter.damp_rate=fwhm2rlx(str2double(cm.children(k).children.data));
                
                % Inform the user
                disp('found line width'); lw=true();
                
            end
            
            % Look for the spin list
            if strcmpi(cm.children(k).name,'spin_names')
                
                % Grab a shorthand
                sn=cm.children(k);
                
                % Look inside
                for m=1:numel(sn.children)
                    
                    % Look for spin names
                    if strcmpi(sn.children(m).name,'spin')
                        
                        % Assign labels
                        idx=strcmpi('index',{sn.children(m).attributes.name});
                        lab=strcmpi('name',{sn.children(m).attributes.name});
                        idx=str2double(sn.children(m).attributes(idx).value);
                        lab=sn.children(m).attributes(lab).value;
                        sys.labels{idx}=['Atom ' lab];
                        
                    end
                    
                end
                
                % Inform the user
                disp('found atom labels');
                
            end
            
            % Look for the shift list
            if strcmpi(cm.children(k).name,'chemical_shifts_ppm')
                
                % Grab a shorthand
                cs=cm.children(k);
                
                % Look inside
                for m=1:numel(cs.children)
                    
                    % Look for spin names
                    if strcmpi(cs.children(m).name,'cs')
                        
                        % Assign chemical shifts
                        idx=strcmpi('index',{cs.children(m).attributes.name});
                        ppm=strcmpi('ppm',{cs.children(m).attributes.name});
                        idx=str2double(cs.children(m).attributes(idx).value);
                        ppm=str2double(cs.children(m).attributes(ppm).value);
                        inter.zeeman.scalar{idx}=ppm;

                    end
                    
                end
                
                % Inform the user
                disp('found chemical shifts'); shifts=true();
                
            end
            
            % Look for the coupling list
            if strcmpi(cm.children(k).name,'couplings_hz')
                
                % Grab a shorthand
                jc=cm.children(k);
                
                % Generate an empty array
                inter.coupling.scalar=cell(numel(inter.zeeman.scalar),...
                                           numel(inter.zeeman.scalar));
                
                % Look inside
                for m=1:numel(jc.children)
                    
                    % Look for spin names
                    if strcmpi(jc.children(m).name,'coupling')
                        
                        % Assign chemical shifts
                        source=strcmpi('from_index',{jc.children(m).attributes.name});
                        destin=strcmpi('to_index',{jc.children(m).attributes.name});
                        value= strcmpi('value',{jc.children(m).attributes.name});
                        source=str2double(jc.children(m).attributes(source).value);
                        destin=str2double(jc.children(m).attributes(destin).value);
                        value= str2double(jc.children(m).attributes(value ).value);
                        inter.coupling.scalar{source,destin}=value;

                    end
                    
                end
                
                % Inform the user
                disp('found scalar couplings'); jc=true();
                
            end
                    
        end
          
    elseif strcmpi(xml.children(n).name,'coupling_matrix')
        
        % Increment subsystem counter
        current_subsystem=current_subsystem+1;
        
    end
    
end

% Make sure the data is in place
if ~(magnet&&lw&&shifts&&jc)
    error('essential infomation is missing from the XML file');
end

% Set isotopes to protons
sys.isotopes=cell(size(inter.zeeman.scalar));
sys.isotopes(:)={'1H'};

end

% Consistency enforcement
function grumble(filename)
if ~exist(filename,'file')
    error('the file specified was not found.');
end
end

% Perhaps this habit goes back to the primitive belief 
% that the word and the thing, the name and the object,
% are identical. 
%
% Enoch Powell

