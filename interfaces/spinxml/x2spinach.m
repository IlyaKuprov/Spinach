% Reads SpinXML files and forms Spinach data structures. Syntax:
%
%         [sys,inter]=x2spinach(file_name,shielding_refs)
%
% Parameters: 
%
%   shielding_refs - when absolute shielding is provided, this
%                    array must give absolute shielding of the
%                    corresponding nuclei in the reference sub-
%                    stance, e.g. {{'1H',31.5},{'13C',189.7}}.
%                    This is necessary because Spinach requires
%                    chemical shifts rather than shieldings.
%                    If chemical shifts are provided, use an
%                    empty cell array.
%
% WARNING: this function assumes that the SpinXML file has passed
%          the validation against the schema, which may be obtain-
%          ed from http://spindynamics.org/SpinXML.php
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=x2spinach.m>

function [sys,inter]=x2spinach(filename,shielding_refs)

% Check consistency
grumble(filename,shielding_refs);

% Get the data structures started
spin_ids=[]; sys.isotopes={}; 
sys.labels={}; inter.coordinates={};

% Parse the XML
xml=parsexml(filename);

% First pass - discover spins
for n=1:numel(xml.children)
    
    % Look for spin elements
    if strcmpi(xml.children(n).name,'spin')
        
        % Labels and coordinates are optional
        current_spin_label=''; xyz=[];
        
        % Look for spin attributes
        for k=1:numel(xml.children(n).attributes)
            
            % Process attributes
            if strcmpi(xml.children(n).attributes(k).name,'id')
                
                % Process spin id number
                current_spin_id=str2double(xml.children(n).attributes(k).value);
                
                % Check for zero-base indexing
                if current_spin_id==0
                    error('zero-base indexing detected - renumber your spins from 1.');
                end
                
            elseif strcmpi(xml.children(n).attributes(k).name,'isotope')
                
                % Process isotope specifications
                current_spin_isotope=xml.children(n).attributes(k).value;
                
            elseif strcmpi(xml.children(n).attributes(k).name,'label')
                
                % Process text labels
                current_spin_label=xml.children(n).attributes(k).value;
                
            end
            
        end
        
        % Look for coordinate information
        for k=1:numel(xml.children(n).children)
            
            % Parse coordinate information
            if strcmpi(xml.children(n).children(k).name,'coordinates')
                
                % Process Cartesian vector components
                for m=1:numel(xml.children(n).children(k).attributes)
                    if strcmpi(xml.children(n).children(k).attributes(m).name,'x')
                        xyz(1)=str2double(xml.children(n).children(k).attributes(m).value);
                    elseif strcmpi(xml.children(n).children(k).attributes(m).name,'y')
                        xyz(2)=str2double(xml.children(n).children(k).attributes(m).value);
                    elseif strcmpi(xml.children(n).children(k).attributes(m).name,'z')
                        xyz(3)=str2double(xml.children(n).children(k).attributes(m).value);
                    end
                end
                
            end
            
        end
        
        % Form Spinach structures
        spin_ids=[spin_ids current_spin_id]; %#ok<AGROW>
        sys.isotopes=[sys.isotopes {current_spin_isotope}];
        sys.labels=[sys.labels {current_spin_label}];
        inter.coordinates=[inter.coordinates; {xyz}];
        
        % Clear temporary variables
        clear('xyz','current_spin_id','current_spin_isotope','current_spin_label');
        
    end
    
end

% Sort arrays by spin index number
[spin_ids,index]=sort(spin_ids,'ascend');
sys.isotopes=sys.isotopes(index);
sys.labels=sys.labels(index);
inter.coordinates=inter.coordinates(index);

% Check for missing spin ids
if norm(diff(diff(spin_ids)),1)>1e-6
    error('spin id numbers are not sequential integers.');
end
if min(spin_ids)~=1
    error('spin id numbers must start from 1.');
end

% Preallocate interaction arrays
inter.zeeman.matrix=cell(1,numel(sys.isotopes));
inter.zeeman.reference=cell(1,numel(sys.isotopes));
inter.zeeman.label=cell(1,numel(sys.isotopes));
inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes));
inter.coupling.label=cell(numel(sys.isotopes),numel(sys.isotopes));
inter.spinrot.matrix=cell(1,numel(sys.isotopes));
inter.spinrot.label=cell(1,numel(sys.isotopes));

% Get the id list started
interaction_ids=[];

% Second pass - discover interactions
for n=1:numel(xml.children)
    
    % Locate an interaction
    if strcmpi(xml.children(n).name,'interaction')
        
        % Some attributes are optional
        inter_spin_b=[]; inter_label=''; 
        inter_reference=''; dcm=eye(3,3);
        
        % Process interaction attributes
        for k=1:numel(xml.children(n).attributes)
            
            % Match attributes
            if strcmpi(xml.children(n).attributes(k).name,'kind')
                
                % Process interaction type
                inter_kind=xml.children(n).attributes(k).value;
                
            elseif strcmpi(xml.children(n).attributes(k).name,'id')
                
                % Process interaction id
                inter_id=str2double(xml.children(n).attributes(k).value);
                
                % Check for zero-base indexing
                if inter_id==0
                    error('zero-base indexing detected - renumber your interactions from 1.');
                end
                
            elseif strcmpi(xml.children(n).attributes(k).name,'spin_a')
                
                % Process first spin
                inter_spin_a=str2double(xml.children(n).attributes(k).value);
                
            elseif strcmpi(xml.children(n).attributes(k).name,'spin_b')
                
                % Process second spin
                inter_spin_b=str2double(xml.children(n).attributes(k).value);
                
            elseif strcmpi(xml.children(n).attributes(k).name,'units')
                
                % Process units
                inter_units=xml.children(n).attributes(k).value;
                
             elseif strcmpi(xml.children(n).attributes(k).name,'label')
                
                % Process the label
                inter_label=xml.children(n).attributes(k).value;
                
             elseif strcmpi(xml.children(n).attributes(k).name,'reference')
                
                % Process the reference
                inter_reference=xml.children(n).attributes(k).value;
                
            end
            
        end
            
        % Process interaction specification
        for k=1:numel(xml.children(n).children)
            
            % Match specification type
            if strcmpi(xml.children(n).children(k).name,'tensor')
                
                % Loop over attributes
                for m=1:numel(xml.children(n).children(k).attributes)
                    
                    % Get the eigenframe matrix
                    if strcmpi(xml.children(n).children(k).attributes(m).name,'xx')
                        A(1,1)=str2double(xml.children(n).children(k).attributes(m).value);
                    elseif strcmpi(xml.children(n).children(k).attributes(m).name,'xy')
                        A(1,2)=str2double(xml.children(n).children(k).attributes(m).value);
                    elseif strcmpi(xml.children(n).children(k).attributes(m).name,'xz')
                        A(1,3)=str2double(xml.children(n).children(k).attributes(m).value);
                    elseif strcmpi(xml.children(n).children(k).attributes(m).name,'yx')
                        A(2,1)=str2double(xml.children(n).children(k).attributes(m).value);
                    elseif strcmpi(xml.children(n).children(k).attributes(m).name,'yy')
                        A(2,2)=str2double(xml.children(n).children(k).attributes(m).value);
                    elseif strcmpi(xml.children(n).children(k).attributes(m).name,'yz')
                        A(2,3)=str2double(xml.children(n).children(k).attributes(m).value);
                    elseif strcmpi(xml.children(n).children(k).attributes(m).name,'zx')
                        A(3,1)=str2double(xml.children(n).children(k).attributes(m).value);
                    elseif strcmpi(xml.children(n).children(k).attributes(m).name,'zy')
                        A(3,2)=str2double(xml.children(n).children(k).attributes(m).value);
                    elseif strcmpi(xml.children(n).children(k).attributes(m).name,'zz')
                        A(3,3)=str2double(xml.children(n).children(k).attributes(m).value);
                    end
                    
                end
                
            elseif strcmpi(xml.children(n).children(k).name,'scalar')
                
                % Loop over attributes
                for m=1:numel(xml.children(n).children(k).attributes)
                    
                    % Get the eigenframe matrix
                    if strcmpi(xml.children(n).children(k).attributes(m).name,'iso')
                        A=eye(3)*str2double(xml.children(n).children(k).attributes(m).value);
                    end
                    
                end
                
            elseif strcmpi(xml.children(n).children(k).name,'eigenvalues')
                
                % Loop over attributes
                for m=1:numel(xml.children(n).children(k).attributes)
                    
                    % Get the eigenframe matrix
                    if strcmpi(xml.children(n).children(k).attributes(m).name,'xx')
                        A(1,1)=str2double(xml.children(n).children(k).attributes(m).value);
                    elseif strcmpi(xml.children(n).children(k).attributes(m).name,'yy')
                        A(2,2)=str2double(xml.children(n).children(k).attributes(m).value);
                    elseif strcmpi(xml.children(n).children(k).attributes(m).name,'zz')
                        A(3,3)=str2double(xml.children(n).children(k).attributes(m).value);
                    end
                    
                end
                
            elseif strcmpi(xml.children(n).children(k).name,'aniso_asymm')
                
                % Loop over attributes
                for m=1:numel(xml.children(n).children(k).attributes)
                    
                    % Extract amplitude parameters
                    if strcmpi(xml.children(n).children(k).attributes(m).name,'iso')
                        iso=str2double(xml.children(n).children(k).attributes(m).value);
                    elseif strcmpi(xml.children(n).children(k).attributes(m).name,'aniso')
                        an=str2double(xml.children(n).children(k).attributes(m).value);
                    elseif strcmpi(xml.children(n).children(k).attributes(m).name,'asymm')
                        as=str2double(xml.children(n).children(k).attributes(m).value);
                    end
                    
                end
                
                % Get the eigenframe matrix
                A=anas2mat(iso,an,as,0,0,0);
                
                % Clear temporary variables
                clear('iso','an','as');
                
            elseif strcmpi(xml.children(n).children(k).name,'axiality_rhombicity')
                
                % Loop over attributes
                for m=1:numel(xml.children(n).children(k).attributes)
                    
                    % Extract amplitude parameters
                    if strcmpi(xml.children(n).children(k).attributes(m).name,'iso')
                        iso=str2double(xml.children(n).children(k).attributes(m).value);
                    elseif strcmpi(xml.children(n).children(k).attributes(m).name,'ax')
                        ax=str2double(xml.children(n).children(k).attributes(m).value);
                    elseif strcmpi(xml.children(n).children(k).attributes(m).name,'rh')
                        rh=str2double(xml.children(n).children(k).attributes(m).value);
                    end
                    
                end
                
                % Get the eigenframe matrix
                A=axrh2mat(iso,ax,rh,0,0,0);
                
                % Clear temporary variables
                clear('iso','ax','rh')
                
            elseif strcmpi(xml.children(n).children(k).name,'span_skew')
                
                % Loop over attributes
                for m=1:numel(xml.children(n).children(k).attributes)
                    
                    % Extract amplitude parameters
                    if strcmpi(xml.children(n).children(k).attributes(m).name,'iso')
                        iso=str2double(xml.children(n).children(k).attributes(m).value);
                    elseif strcmpi(xml.children(n).children(k).attributes(m).name,'span')
                        sp=str2double(xml.children(n).children(k).attributes(m).value);
                    elseif strcmpi(xml.children(n).children(k).attributes(m).name,'skew')
                        sk=str2double(xml.children(n).children(k).attributes(m).value);
                    end
                    
                end
                
                % Get the eigenframe matrix
                A=spsk2mat(iso,sp,sk,0,0,0);
                
                % Clear temporary variables
                clear('iso','sp','sk');
            
            % Presence of the orientation tag causes DCM to be built
            elseif strcmpi(xml.children(n).children(k).name,'orientation')
                
                % Process rotation specifications
                for m=1:numel(xml.children(n).children(k).children)
                    
                    % Match specification type
                    if strcmpi(xml.children(n).children(k).children(m).name,'euler_angles')
                        
                        % Loop over attributes
                        for h=1:numel(xml.children(n).children(k).children(m).attributes)
                    
                            % Extract Euler angles
                            if strcmpi(xml.children(n).children(k).children(m).attributes(h).name,'alpha')
                                alp=str2double(xml.children(n).children(k).children(m).attributes(h).value);
                            elseif strcmpi(xml.children(n).children(k).children(m).attributes(h).name,'beta')
                                bet=str2double(xml.children(n).children(k).children(m).attributes(h).value);
                            elseif strcmpi(xml.children(n).children(k).children(m).attributes(h).name,'gamma')
                                gam=str2double(xml.children(n).children(k).children(m).attributes(h).value);
                            end
                            
                        end
                    
                        % Get the DCM
                        dcm=euler2dcm(alp,bet,gam);
                        
                        % Clear temporary variables
                        clear('alp','bet','gam');
                        
                    elseif strcmpi(xml.children(n).children(k).children(m).name,'quaternion')
                        
                        % Loop over attributes
                        for h=1:numel(xml.children(n).children(k).children(m).attributes)
                    
                            % Extract the quaternion
                            if strcmpi(xml.children(n).children(k).children(m).attributes(h).name,'re')
                                q.u=str2double(xml.children(n).children(k).children(m).attributes(h).value);
                            elseif strcmpi(xml.children(n).children(k).children(m).attributes(h).name,'i')
                                q.i=str2double(xml.children(n).children(k).children(m).attributes(h).value);
                            elseif strcmpi(xml.children(n).children(k).children(m).attributes(h).name,'j')
                                q.j=str2double(xml.children(n).children(k).children(m).attributes(h).value);
                            elseif strcmpi(xml.children(n).children(k).children(m).attributes(h).name,'k')
                                q.k=str2double(xml.children(n).children(k).children(m).attributes(h).value);
                            end
                            
                        end
                    
                        % Get the DCM
                        dcm=quat2dcm([q.u q.i q.j q.k]);
                        
                        % Clear temporary variables
                        clear('q');
                        
                    elseif strcmpi(xml.children(n).children(k).children(m).name,'dcm')
                        
                        % Loop over attributes
                        for h=1:numel(xml.children(n).children(k).children(m).attributes)
                    
                            % Extract the DCM
                            if strcmpi(xml.children(n).children(k).children(m).attributes(h).name,'xx')
                                dcm(1,1)=str2double(xml.children(n).children(k).children(m).attributes(h).value);
                            elseif strcmpi(xml.children(n).children(k).children(m).attributes(h).name,'xy')
                                dcm(1,2)=str2double(xml.children(n).children(k).children(m).attributes(h).value);
                            elseif strcmpi(xml.children(n).children(k).children(m).attributes(h).name,'xz')
                                dcm(1,3)=str2double(xml.children(n).children(k).children(m).attributes(h).value);
                            elseif strcmpi(xml.children(n).children(k).children(m).attributes(h).name,'yx')
                                dcm(2,1)=str2double(xml.children(n).children(k).children(m).attributes(h).value);
                            elseif strcmpi(xml.children(n).children(k).children(m).attributes(h).name,'yy')
                                dcm(2,2)=str2double(xml.children(n).children(k).children(m).attributes(h).value);
                            elseif strcmpi(xml.children(n).children(k).children(m).attributes(h).name,'yz')
                                dcm(2,3)=str2double(xml.children(n).children(k).children(m).attributes(h).value);
                            elseif strcmpi(xml.children(n).children(k).children(m).attributes(h).name,'zx')
                                dcm(3,1)=str2double(xml.children(n).children(k).children(m).attributes(h).value);
                            elseif strcmpi(xml.children(n).children(k).children(m).attributes(h).name,'zy')
                                dcm(3,2)=str2double(xml.children(n).children(k).children(m).attributes(h).value);
                            elseif strcmpi(xml.children(n).children(k).children(m).attributes(h).name,'zz')
                                dcm(3,3)=str2double(xml.children(n).children(k).children(m).attributes(h).value);
                            end
                            
                        end
                        
                    elseif strcmpi(xml.children(n).children(k).children(m).name,'angle_axis')
                        
                        % Process angle-axis specification
                        for h=1:numel(xml.children(n).children(k).children(m).children)
                            
                            % Match children
                            if strcmpi(xml.children(n).children(k).children(m).children(h).name,'angle')
                                
                                % Loop over attributes
                                for j=1:numel(xml.children(n).children(k).children(m).children(h).attributes)
                                    
                                    % Get the phi angle
                                    if strcmpi(xml.children(n).children(k).children(m).children(h).attributes(j).name,'phi')
                                        phi=str2double(xml.children(n).children(k).children(m).children(h).attributes(j).value);
                                    end
                                    
                                end
                                
                            elseif strcmpi(xml.children(n).children(k).children(m).children(h).name,'axis')
                                
                                % Loop over attributes
                                for j=1:numel(xml.children(n).children(k).children(m).children(h).attributes)
                                    
                                    % Get the rotation axis
                                    if strcmpi(xml.children(n).children(k).children(m).children(h).attributes(j).name,'x')
                                        rx=str2double(xml.children(n).children(k).children(m).children(h).attributes(j).value);
                                    elseif strcmpi(xml.children(n).children(k).children(m).children(h).attributes(j).name,'y')
                                        ry=str2double(xml.children(n).children(k).children(m).children(h).attributes(j).value);
                                    elseif strcmpi(xml.children(n).children(k).children(m).children(h).attributes(j).name,'z')
                                        rz=str2double(xml.children(n).children(k).children(m).children(h).attributes(j).value);
                                    end
                                    
                                end
                                
                            end
                            
                        end
                                
                        % Get the DCM
                        dcm=anax2dcm([rx ry rz],phi);
                        
                        % Clear temporary variables
                        clear('rx','ry','rz','phi');
                        
                    end
                       
                end
                
            end

        end
        
        % Check the DCM
        if norm(dcm'*dcm-eye(3),'fro')>1e-2, error('the DCM specified is not an orthogonal matrix.'); end
        if abs(det(dcm)-1)>1e-2,             error('the DCM specified has a non-unit determinant.'); end
        
        % Rotate the tensor
        A=dcm*A*dcm';
        
        % Process interaction type
        switch inter_kind
            
            case 'hfc'
                
                % Make sure coordinates are not specified
                if (~isempty(inter.coordinates{inter_spin_a}))&&...
                   (~isempty(inter.coordinates{inter_spin_b}))
                    warning('both particles have coordinates, you might be counting the hyperfine coupling twice.');
                end
                
                % Make sure spins are different
                if inter_spin_a==inter_spin_b
                    error('hyperfine coupling specified between a spin and itself.');
                end
                
                % Start with zero if previously empty
                if isempty(inter.coupling.matrix{inter_spin_a,inter_spin_b})
                    inter.coupling.matrix{inter_spin_a,inter_spin_b}=zeros(3,3);
                end
                
                % Assign hyperfine coupling
                switch inter_units
                    case 'Hz'
                        inter.coupling.matrix{inter_spin_a,inter_spin_b}=inter.coupling.matrix{inter_spin_a,inter_spin_b}+A;
                    case 'MHz'
                        inter.coupling.matrix{inter_spin_a,inter_spin_b}=inter.coupling.matrix{inter_spin_a,inter_spin_b}+1e6*A;
                    case 'gauss'
                        inter.coupling.matrix{inter_spin_a,inter_spin_b}=inter.coupling.matrix{inter_spin_a,inter_spin_b}+1e6*gauss2mhz(A);
                    otherwise
                        error('unknown hyperfine coupling units.');
                end
                
                % Assign label
                inter.coupling.label{inter_spin_a,inter_spin_b}=inter_label;
            
            case 'shielding'
                
                % Only accept ppm
                if ~strcmp(inter_units,'ppm')
                    error('unknown chemical shielding units.');
                end
                
                % Catch double spec
                if ~isempty(inter.zeeman.matrix{inter_spin_a})
                    error('chemical shielding specified twice.');
                end
                
                % Catch non-empty B spin
                if ~isempty(inter_spin_b)
                    error('B spin specified for a Zeeman interaction.');
                end
                
                % Translate shielding into shift
                isotope=sys.isotopes{inter_spin_a};
                for k=1:numel(shielding_refs)
                    if strcmp(isotope,shielding_refs{k}{1})
                        A=eye(3)*shielding_refs{k}{2}-A;
                    end
                end
                
                % Assign chemical shift
                inter.zeeman.matrix{inter_spin_a}=A;
                
                % Assign label
                inter.zeeman.label{inter_spin_a}=inter_label;
                
                % Assign reference
                inter.zeeman.reference{inter_spin_a}=inter_reference;
            
            case 'shift'
                
                % Only accept ppm
                if ~strcmp(inter_units,'ppm')
                    error('unknown chemical shift units.');
                end
                
                % Catch double spec
                if ~isempty(inter.zeeman.matrix{inter_spin_a})
                    error('chemical shift specified twice.');
                end
                
                % Catch non-empty B spin
                if ~isempty(inter_spin_b)
                    error('B spin specified for a Zeeman interaction.');
                end
                
                % Assign chemical shift
                inter.zeeman.matrix{inter_spin_a}=A;
                
                % Assign label
                inter.zeeman.label{inter_spin_a}=inter_label;
                
                % Assign reference
                inter.zeeman.reference{inter_spin_a}=inter_reference;
            
            case 'dipolar'
                
                % Make sure coordinates are not specified
                if (~isempty(inter.coordinates{inter_spin_a}))&&...
                   (~isempty(inter.coordinates{inter_spin_b}))
                    warning('both particles have coordinates, you might be counting the dipolar coupling twice.');
                end
                
                % Make sure spins are different
                if inter_spin_a==inter_spin_b
                    error('dipolar coupling specified between a spin and itself.');
                end
                
                % Make sure the matrix is traceless
                if abs(trace(A))>1e-2, error('coupling declared to be dipolar is not traceless.'); end
                
                % Start with zero if previously empty
                if isempty(inter.coupling.matrix{inter_spin_a,inter_spin_b})
                    inter.coupling.matrix{inter_spin_a,inter_spin_b}=zeros(3,3);
                end
                
                % Assign the dipolar coupling
                switch inter_units
                    case 'Hz'
                        inter.coupling.matrix{inter_spin_a,inter_spin_b}=inter.coupling.matrix{inter_spin_a,inter_spin_b}+A;
                    case 'kHz'
                        inter.coupling.matrix{inter_spin_a,inter_spin_b}=inter.coupling.matrix{inter_spin_a,inter_spin_b}+1e3*A;
                    otherwise
                        error('unknown dipolar coupling units.');
                end
                
                % Assign label
                inter.coupling.label{inter_spin_a,inter_spin_b}=inter_label;
            
            case 'quadrupolar'
                
                % Make sure the spin is big enough
                [~,mult]=spin(sys.isotopes{inter_spin_a});
                if mult<3, error('quadratic coupling specified for a spin-1/2 particle.'); end
                
                % Catch double spec
                if ~isempty(inter.coupling.matrix{inter_spin_a,inter_spin_a})
                    error('quadrupolar interaction specified twice.');
                end
                
                % Catch non-empty B spin
                if ~isempty(inter_spin_b)
                    error('B spin specified for a quadrupolar coupling.');
                end
                
                % Make sure the matrix is traceless
                if abs(trace(A))>1e-2, error('coupling declared to be quadrupolar is not traceless.'); end
                
                % Assign the quadrupolar coupling
                switch inter_units
                    case 'kHz'
                        inter.coupling.matrix{inter_spin_a,inter_spin_a}=1e3*A;
                    case 'MHz'
                        inter.coupling.matrix{inter_spin_a,inter_spin_a}=1e6*A;
                    otherwise
                        error('unknown quadrupolar coupling units.');
                end
                
                % Assign label
                inter.coupling.label{inter_spin_a,inter_spin_a}=inter_label;
            
            case 'jcoupling'
                
                % Make sure spins are different
                if inter_spin_a==inter_spin_b
                    error('J-coupling specified between a spin and itself.');
                end
                
                % Start with zero if previously empty
                if isempty(inter.coupling.matrix{inter_spin_a,inter_spin_b})
                    inter.coupling.matrix{inter_spin_a,inter_spin_b}=zeros(3,3);
                end
                
                % Assign the J-coupling
                switch inter_units
                    case 'Hz'
                        inter.coupling.matrix{inter_spin_a,inter_spin_b}=A;
                    otherwise
                        error('unknown J-coupling units.');
                end
                
                % Assign label
                inter.coupling.label{inter_spin_a,inter_spin_b}=inter_label; 
            
            case 'gtensor'
                
                % Make sure the particle is an electron
                if ~regexp(sys.isotopes{inter_spin_a},'^E\d')
                    error('g-tensor specified for something that is not an electron.');
                end
                
                % Catch double spec
                if ~isempty(inter.zeeman.matrix{inter_spin_a})
                    error('g-tensor specified twice.');
                end
                
                % Catch non-empty B spin
                if ~isempty(inter_spin_b)
                    error('B spin specified for a Zeeman interaction.');
                end
                
                % Assign the g-tensor
                switch inter_units
                    case 'bohr'
                        inter.zeeman.matrix{inter_spin_a}=A;
                    otherwise
                        error('unknown g-tensor units.');
                end
                
                % Assign label
                inter.zeeman.label{inter_spin_a}=inter_label;
            
            case 'zfs'
                
                % Make sure the particle is an electron
                if ~regexp(sys.isotopes{inter_spin_a},'^E\d')
                    error('ZFS specified for something that is not an electron.');
                end
                
                % Make sure the spin is big enough
                [~,mult]=spin(sys.isotopes{inter_spin_a});
                if mult<3, error('quadratic coupling specified for a spin-1/2 particle.'); end
                
                % Catch double spec
                if ~isempty(inter.coupling.matrix{inter_spin_a,inter_spin_a})
                    error('ZFS specified twice.');
                end
                
                % Catch non-empty B spin
                if ~isempty(inter_spin_b)
                    error('B spin specified for a zero-field splitting.');
                end
                
                % Make sure the matrix is traceless
                if abs(trace(A))>1e-2, error('coupling declared to be ZFS is not traceless.'); end
                
                % Assign ZFS
                switch inter_units
                    case 'MHz'
                        inter.coupling.matrix{inter_spin_a,inter_spin_a}=1e6*A;
                    case 'GHz'
                        inter.coupling.matrix{inter_spin_a,inter_spin_a}=1e9*A;
                    otherwise
                        error('unknown ZFS units.');
                end
                
                % Assign label
                inter.coupling.label{inter_spin_a,inter_spin_a}=inter_label;
            
            case 'exchange'
                
                % Make sure both particles are electrons
                if (~regexp(sys.isotopes{inter_spin_a},'^E\d'))||...
                   (~regexp(sys.isotopes{inter_spin_b},'^E\d'))
                    error('exchange coupling specified for something that is not an electron.');
                end
                
                % Make sure spins are different
                if inter_spin_a==inter_spin_b
                    error('exchange coupling specified between a spin and itself.');
                end
                
                % Start with zero if previously empty
                if isempty(inter.coupling.matrix{inter_spin_a,inter_spin_b})
                    inter.coupling.matrix{inter_spin_a,inter_spin_b}=zeros(3,3);
                end
                
                % Assign the exchange coupling
                switch inter_units
                    case 'MHz'
                        inter.coupling.matrix{inter_spin_a,inter_spin_b}=1e6*A;
                    case 'GHz'
                        inter.coupling.matrix{inter_spin_a,inter_spin_b}=1e9*A;
                    otherwise
                        error('unknown exchange coupling units.');
                end
                
                % Assign label
                inter.coupling.label{inter_spin_a,inter_spin_b}=inter_label;
            
            case 'spinrotation'
                
                % Catch double spec
                if ~isempty(inter.spinrot.matrix{inter_spin_a})
                    error('spin-rotation tensor specified twice.');
                end
                
                % Catch non-empty B spin
                if ~isempty(inter_spin_b)
                    error('B spin specified for a spin-rotation coupling.');
                end
                
                % Assign the spin-rotation tensor
                switch inter_units
                    case 'Hz'
                        inter.spinrot.matrix{inter_spin_a}=A;
                    case 'kHz'
                        inter.spinrot.matrix{inter_spin_a}=1e3*A;
                    case 'GHz'
                        inter.spinrot.matrix{inter_spin_a}=1e6*A;
                    otherwise
                        error('unknown spin-rotation tensor units.');
                end
            
            otherwise
                
                % Complain and bomb out
                error('unrecognized spin interaction.');
        
        end
        
        % Add the id to the list
        interaction_ids=[interaction_ids inter_id]; %#ok<AGROW>
        
        % Clear temporary variables
        clear('inter_kind','inter_id','inter_units','inter_spin_a','inter_units',...
              'inter_spin_b','inter_label','inter_reference','A','dcm');

    end
    
    % Check for missing interaction ids
    if norm(diff(diff(sort(interaction_ids,'ascend'))),1)>1e-6
        warning('interaction id numbers are not sequential integers.');
    end

end

end

% Consistency enforcement
function grumble(filename,shielding_refs)
if ~exist(filename,'file')
    error('the specified file was not found.');
end
if ~iscell(shielding_refs)
    error('shielding_refs must be a cell array.');
end
for n=1:numel(shielding_refs)
    if ~iscell(shielding_refs{n})
        error('elements of shielding_refs must be cell arrays.');
    end
    if numel(shielding_refs{n})~=2
        error('elements of shielding_refs must have exactly two sub-elements each.');
    end
    if ~ischar(shielding_refs{n}{1})
        error('the first part of each element of shielding_refs must be a character string.');
    end
    if ~isnumeric(shielding_refs{n}{2})
        error('the second part of each element of parameters.rframes must be a number.');
    end
end
end

% People who think they know everything are a great 
% annoyance to those of us who do.
%
% Isaac Asimov

