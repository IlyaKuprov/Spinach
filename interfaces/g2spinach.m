% Makes Spinach data structures from parsed outputs of electronic
% structure theory packages, such as Gaussian and ORCA. Syntax:
%
%    [sys,inter]=g2spinach(props,elements,references,options)
%
% Parameters:
%
%    props      - the output of gparse() or oparse(); isotopes
%                 are taken from props.isotopes, which reports
%                 what the log specifies or assumes for each
%                 atom; interaction tensors are imported as
%                 the log reports them for those isotopes.
%                 When the data cannot be interpreted without
%                 a guess (no isotope or a non-magnetic one
%                 for a requested atom, a quadrupolar tensor
%                 on an isotope without a quadrupole moment,
%                 printed J-couplings that K-couplings and
%                 the isotopes do not reproduce), an error
%                 names the field of props to be edited.
%                 Use isoswap() to replace isotopes after the
%                 import.
%
%    elements   - a cell array of element symbols:
%
%                             {'H','N',...}
%
%                 giving the list of elements that should be
%                 imported. If the list contains an electron,
%                 e.g. {'E','H',...}, then EPR mode is assumed:
%                 chemical shielding and scalar couplings are
%                 ignored, but g-tensor and hyperfine couplings
%                 are included.
%
%    references - a vector of absolute shielding values, one
%                 per element, for the reference substances
%                 that are to be placed at zero ppm; you ne-
%                 ed to run separate electronic structure 
%                 theory calculations for those substances
%                 wuth the same method. Absolute isotropic 
%                 shielding values for tetramethylsilane in
%                 vacuum are:
% 
%                  GIAO                  13C       1H 
%                  B3LYP/6-31G*        189.6621  32.1833 
%                  B3LYP/6-311+G(2d,p) 182.4485  31.8201 
%                  HF/6-31G*           199.9711  32.5957 
%                  HF/6-311+G(2d,p)    192.5828  32.0710
%
%                  CSGT                  13C       1H 
%                  B3LYP/6-31G*        188.5603  29.1952 
%                  B3LYP/6-311+G(2d,p) 182.1386  31.7788 
%                  HF/6-31G*           196.8670  29.5517 
%                  HF/6-311+G(2d,p)    192.5701  31.5989
%
%                 This setting is ignored when electrons are
%                 present in the system.
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
% Outputs:
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

function [sys,inter]=g2spinach(props,elements,references,options)

% Set default options
if ~exist('options','var'), options=[]; end

% Check consistency
grumble(props,elements,references,options);

% Fundamental constants
nuclear_magneton=7.6225932291E6;

% Index the atoms to include and take their isotopes from the log
[selected,ref_index]=ismember(props.symbols,elements);
index=find(selected); ref_index=ref_index(selected);
sys.isotopes=reshape(props.isotopes(index),1,[]);
nspins=numel(index);

% Warn about zero-field splitting
if isfield(props,'zfs')
    disp('WARNING: zero-field splitting in the log has not been imported.');
end

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
switch ismember('E',elements)
    
    % EPR parameterization
    case 1
        
        % Add the electron as the last spin in the isotope list
        sys.isotopes=[sys.isotopes, 'E']; nspins=nspins+1;
        
        % Electron coordinates should be treated as unknown
        if include_xyz, inter.coordinates{end+1}=[]; end
        
        % All Zeeman tensors are zero except for electron g-tensor
        inter.zeeman.matrix=cell(1,nspins);
        inter.zeeman.matrix{nspins}=props.g_tensor.matrix;
        disp('WARNING: nuclear chemical shifts were set to zero.');
        disp('WARNING: offset specification has been ignored.');
        
        % All couplings are zero except for the hyperfine couplings to the electron
        inter.coupling.matrix=mat2cell(zeros(3*nspins,3*nspins),3*ones(nspins,1),3*ones(nspins,1));
        for n=1:(nspins-1)
            inter.coupling.matrix{n,end}=1e6*gauss2mhz(props.hfc.full.matrix{index(n)}/2);
            inter.coupling.matrix{end,n}=1e6*gauss2mhz(props.hfc.full.matrix{index(n)}/2);
        end
        
        % Warn about nuclear interactions that are ignored in EPR mode
        if isfield(props,'nqi')&&any(cellfun(@(x)any(x,'all'),props.nqi(index)))
            disp('WARNING: nuclear quadrupolar couplings in the log have been ignored.');
        end
        if isfield(props,'k_couplings')||isfield(props,'j_couplings')
            disp('WARNING: scalar couplings in the log have been ignored.');
        end
        if isfield(props,'srt')
            disp('WARNING: spin-rotation couplings in the log have been ignored.');
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
            if include_xyz, inter.coordinates(killing_pattern)=[]; end
        end
        
    % NMR parameterization    
    case 0
        
        % Reference to bare nuclei if the user did not supply any reference values
        if (nargin<3)||(~any(references)), references=zeros(size(elements)); end
        
        % Warn about electron interactions that are ignored in NMR mode
        if isfield(props,'hfc')||isfield(props,'g_tensor')
            disp('WARNING: hyperfine couplings and g-tensor in the log have been ignored.');
        end
        
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
        if isfield(props,'srt')
            inter.spinrot.matrix=cell(1,nspins);
            for n=1:nspins
                inter.spinrot.matrix{n}=props.srt{index(n)};
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
            
            % Refuse printed J-couplings that differ by more than the scatter of gamma tables (0.14%)
            if isfield(props,'j_couplings')
                printed=props.j_couplings(index,index)/2; rebuilt=inter.coupling.scalar;
                stated=(printed~=0)&(rebuilt~=0);
                if any(abs(printed(stated)-rebuilt(stated))>3e-3*abs(printed(stated)))
                    error(['J-couplings printed in the log do not follow from its K-couplings for the isotopes in ' ...
                           'props.isotopes: set props.isotopes to the isotopes of the calculation, or remove ' ...
                           'props.j_couplings to rebuild J-couplings for the isotopes in props.isotopes.']);
                end
            end
            
            % Drop insignificant J-couplings
            if (nargin==4)&&isfield(options,'min_j')
                inter.coupling.scalar=inter.coupling.scalar.*(abs(inter.coupling.scalar)>options.min_j);
            end
            
            % Convert into a cell array
            inter.coupling.scalar=num2cell(inter.coupling.scalar);
            
        elseif isfield(props,'j_couplings')

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
function grumble(props,elements,references,options) %#ok<INUSD>
if ~isstruct(props)
    error('the first argument must be a structure returned by gparse() or oparse().');
end
if (~iscell(elements))||any(~cellfun(@ischar,elements))
    error('elements must be a cell array of element symbols, e.g. {''H'',''C''}.');
end
if (~isnumeric(references))||(numel(elements)~=numel(references))
    error('references must be a numerical array with the same number of entries as elements.');
end
if (~isfield(props,'isotopes'))||(~iscell(props.isotopes))||(numel(props.isotopes)~=numel(props.symbols))
    error('props.isotopes must be a cell array of isotope strings, e.g. ''13C'', with one entry per atom of props.symbols.');
end
selected=find(ismember(props.symbols,elements));
unknown=selected(cellfun(@isempty,props.isotopes(selected)));
if ~isempty(unknown)
    error(['the log does not identify the isotope of atom ' int2str(unknown(1)) ' (' props.symbols{unknown(1)} '): set props.isotopes{' ...
           int2str(unknown(1)) '}, e.g. to ''13C'', for every such atom.']);
end
silent=selected(cellfun(@spin,props.isotopes(selected))==0);
if ~isempty(silent)
    error(['the isotope of atom ' int2str(silent(1)) ' (' props.isotopes{silent(1)} ') is non-magnetic: set props.isotopes{' ...
           int2str(silent(1)) '} to the magnetic isotope to be imported, for every such atom.']);
end
if (~ismember('E',elements))&&isfield(props,'nqi')
    for n=selected(:)'
        [~,mult]=spin(props.isotopes{n});
        if (mult<3)&&any(props.nqi{n},'all')
            error(['atom ' int2str(n) ' has a quadrupolar tensor in props.nqi, but ' props.isotopes{n} ' has no quadrupole moment: set ' ...
                   'props.isotopes{' int2str(n) '} to the isotope of that tensor, or empty props.nqi{' int2str(n) '}.']);
        end
    end
end
end

% ACHTUNG! ALLES NONTECHNISCHEN LOOKENSPEEPERS!
% Das Computermachine ist nicht fur gefingerpoken und mittengrabben.
% Ist easy schnappen der Springenwerk, blowenfusen und poppencorken
% mit Spitzensparken. Ist nicht fur gewerken bei das Dumpkopfen. Die
% rubbernecken Sichtseeren keepen Hands in die Pockets muss, relaxen
% und watchen die Blinkenlichten.
%
% A sign, first reported in 1955
% at an IBM computing facility


