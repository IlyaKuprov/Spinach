% Makes Spinach data structures from parsed outputs of electronic
% structure theory packages, such as Gaussian and ORCA. Syntax:
%
%     [sys,inter]=g2spinach(props,nuclei,references,options)
%
% where props is the output of gparse and nuclei is a cell array 
% of the following form:
%
%                    {{'H','1H'},{'N','15N'}...}
%
% giving the list of elements and isotopes that should be imported. 
%
% References is a vector of absolute shielding values for your refe-
% rence substances that are to be placed at zero ppm chemical shift:
% you need to run separate electronic structure theory calculations
% for those substances wuth the same method. Absolute isotropic shi-
% elding values for tetramethylsilane in vacuum are:
% 
%   GIAO                  13C       1H 
%   B3LYP/6-31G*        189.6621  32.1833 
%   B3LYP/6-311+G(2d,p) 182.4485  31.8201 
%   HF/6-31G*           199.9711  32.5957 
%   HF/6-311+G(2d,p)    192.5828  32.0710
%
%   CSGT                  13C       1H 
%   B3LYP/6-31G*        188.5603  29.1952 
%   B3LYP/6-311+G(2d,p) 182.1386  31.7788 
%   HF/6-31G*           196.8670  29.5517 
%   HF/6-311+G(2d,p)    192.5701  31.5989 
%
% If the isotope list contains an electron, e.g.
%
%                     {{'E','E'},{'H','1H'}...}
%
% then EPR mode is assumed - chemical shielding and scalar couplings
% are ignored, but g-tensor and hyperfine couplings are included. A
% spin-1/2 electron is assumed in this case.
%
% The following options are currently available:
%
%     options.min_j    -  scalar coupling threshold in Hz. J-coup-
%                         lings smaller than this value will be 
%                         ignored in the NMR mode.
%
%     options.min_hfc  -  hyperfine coupling threshold in Hz. Hy-
%                         perfine tensors with a Frobenius norm 
%                         smaller than this value will be ignored
%                         in the EPR mode.
%
%     options.purge    -  if set to 'on' in EPR mode, removes the
%                         spins with hyperfine coupling below 
%                         options.min_hfc from the spin system.
%
%     options.no_xyz   -  if set to 1, causes the function to ig-
%                         nore the coordinate information and
%                         only keep the interaction tensors
%
% The following parameters are returned if the corresponding informa-
% tion is found in the input:
%
%     sys.isotopes           Nspins x 1 cell array of strings
%
%     inter.coordinates      Nspins x 3 dense matrix, Angstrom. Not
%                            returned if there is an electron in the
%                            isotope list (in the EPR case it is not
%                            a good idea to use the molecular 
%                            coordinates for spins).
%
%     inter.zeeman.matrix    Nspins x 1 cell array of 3x3 matrices,
%                            ppm for nuclei, g-tensor for electrons.
%                            Zero interactions have zero matrices.
%
%     inter.coupling.matrix  Nspins x Nspins cell array of 3x3 mat-
%                            rices, all in Hz. Zero interactions have
%                            zero matrices.
%
%     inter.coupling.scalar  Nspins x Nspins cell array of scalar
%                            couplings, all in Hz. Zero couplings are
%                            returned as zeros.
%
%     inter.spinrot.matrix   spin-rotation coupling tensors for
%                            each nucleus
%
% janm@umbc.edu
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=g2spinach.m>

function [sys,inter]=g2spinach(props,particles,references,options)

% Set default options
if ~exist('options','var'), options=[]; end

% Check consistency
grumble(props,particles,references,options);

% Fundamental constants
nuclear_magneton=7.6225932291E6;

% Index the particles to include
sys.isotopes={}; index=[]; ref_index=[];
for n=1:length(props.symbols) %#ok<*AGROW>
    for k=1:length(particles) 
        if strcmp(props.symbols{n},particles{k}{1})
            sys.isotopes=[sys.isotopes particles{k}{2}];
            index=[index n]; ref_index=[ref_index k];
        end
    end
end
nspins=length(index);

% Decide whether to include coordinates
if nargin<4
    include_xyz=true();
else
    if ~isfield(options,'no_xyz')
        include_xyz=true();
    elseif ~options.no_xyz
        include_xyz=true();
    else
        include_xyz=false();
    end
end

% Process coordinates
if include_xyz
    inter.coordinates=props.std_geom(index,:);
    inter.coordinates=num2cell(inter.coordinates,2);
end

% Decide which magnetic parameters to return
switch ismember('E',[particles{:}])
    
    % EPR parameterization
    case 1
        
        % Add the electron as the last spin in the isotope list
        sys.isotopes=[sys.isotopes, 'E']; nspins=nspins+1;
        
        % Electron coordinates should be treated as unknown
        if include_xyz, inter.coordinates{end+1}=[]; end
        
        % All Zeeman tensors are zero except for the g-tensor of the electron
        inter.zeeman.matrix=cell(1,nspins);
        inter.zeeman.matrix{nspins}=props.g_tensor.matrix;
        
        % All couplings are zero except for the hyperfine couplings to the electron
        inter.coupling.matrix=mat2cell(zeros(3*nspins,3*nspins),3*ones(nspins,1),3*ones(nspins,1));
        for n=1:(nspins-1)
            inter.coupling.matrix{n,end}=1e6*gauss2mhz(props.hfc.full.matrix{index(n)}/2);
            inter.coupling.matrix{end,n}=1e6*gauss2mhz(props.hfc.full.matrix{index(n)}/2);
        end
        
        % Remove small hyperfine couplings
        if (nargin==4)&&isfield(options,'min_hfc')
            for n=1:nspins
                for k=1:nspins
                    if norm(inter.coupling.matrix{n,k},'fro')<options.min_hfc
                        inter.coupling.matrix{n,k}=[];
                    end
                end
            end
        end
        
        % Remove the nuclei that are not coupled to the electron
        if (nargin==4)&&isfield(options,'purge')&&strcmp(options.purge,'on')
            killing_pattern=cellfun(@isempty,inter.coupling.matrix(:,nspins)); killing_pattern(end)=0;
            inter.coupling.matrix(:,killing_pattern)=[];
            inter.coupling.matrix(killing_pattern,:)=[];
            inter.zeeman.matrix(killing_pattern)=[];
            sys.isotopes(killing_pattern)=[];
        end
        
    % NMR parameterization    
    case 0
        
        % Reference to bare nuclei if the user did not supply any reference values
        if (nargin<3)||(~any(references)), references=zeros(size(particles)); end
        
        % Absorb Zeeman tensors
        if isfield(props,'cst')
            inter.zeeman.matrix=cell(1,nspins);
            for n=1:nspins
                inter.zeeman.matrix{n}=-props.cst{index(n)}+eye(3)*references(ref_index(n));
            end
        end
        
        % Absorb quadrupolar couplings
        if isfield(props,'nqi')
            inter.coupling.matrix=cell(nspins,nspins);
            for n=1:nspins
                [~,mult]=spin(sys.isotopes{n}); I=(mult-1)/2;
                if mult>2
                    inter.coupling.matrix{n,n}=props.nqi{index(n)}/(2*I*(2*I-1));
                end
            end
        end
        
        % Absorb spin-rotation couplings
        if isfield(props,'src')
            inter.spinrot.matrix=cell(nspins);
            for n=1:nspins
                inter.spinrot.matrix{n}=props.src{index(n)};
            end
        end
        
        % Absorb and prune scalar couplings
        if isfield(props,'k_couplings')
            
            % Get conversion prefactor
            prefactor=1/(2*pi*nuclear_magneton)^2;
            
            % Loop over A spins
            for n=1:nspins
                
                % Magnetogyric ratio, A
                gamma_a=spin(sys.isotopes{n});
                
                % Loop over B spins
                for m=1:nspins
                    
                    % Magnetogyric ratio, B
                    gamma_b=spin(sys.isotopes{m});
                    
                    % Convert isotope-independent K-couplings into isotope-specific J-couplings
                    inter.coupling.scalar(n,m)=prefactor*gamma_a*gamma_b*props.k_couplings(index(n),index(m))/2;
                    
                end
                
            end
            
            % Drop insignificant J-couplings
            if (nargin==4)&&isfield(options,'min_j')
                inter.coupling.scalar=inter.coupling.scalar.*(abs(inter.coupling.scalar)>options.min_j);
            end
            
            % Convert into a cell array
            inter.coupling.scalar=num2cell(inter.coupling.scalar);
            
        elseif isfield(props,'j_couplings')

            % Warn the user about isotope scaling
            disp('g2spinach: WARNING - J-couplings are isotope-specific, make sure');
            disp('g2spinach: your electronic structure theory package had assumed');
            disp('g2spinach: correct isotopes; scale by gamma ratio if not.');
            
            % Absorb J-couplings and drop insignificant ones
            inter.coupling.scalar=props.j_couplings(index,index)/2;
            if (nargin==4)&&isfield(options,'min_j')
                inter.coupling.scalar=inter.coupling.scalar.*(abs(inter.coupling.scalar)>options.min_j);
            end

            % Convert into a cell array
            inter.coupling.scalar=num2cell(inter.coupling.scalar);

        end
       
end

end

% Consistency enforcement
function grumble(props,nuclei,references,options) %#ok<INUSD>
if ~isstruct(props)
    error('the first argument must be a structure returned by gparse().');
end
if ~iscell(nuclei)
    error('the second argument must have the form {{''H'',''1H''},{''C'',''13C''},...}');
end
if (~isnumeric(references))||(numel(nuclei)~=numel(references))
    error('references must be a numerical array with the same number of entries as nuclei.');
end
end

% ACHTUNG! ALLES LOOKENSPEEPERS
% Das Computermachine ist nicht fur gefingerpoken und mittengrabben.
% Ist easy schnappen der Springenwerk, blowenfusen und poppencorken
% mit Spitzensparken. Ist nicht fur gewerken bei das Dumpkopfen. Die
% rubbernecken Sichtseeren keepen Hands in die Pockets muss, relaxen
% und watchen die Blinkenlichten.

