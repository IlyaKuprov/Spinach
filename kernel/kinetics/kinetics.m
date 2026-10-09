% Direct-sum chemical kinetics from explicit reaction records. Syntax:
%
%                  K=kinetics(spin_system)
%                  K=kinetics(spin_system,'report')
%
% Parameters:
%
%    spin_system - Spinach object with chem.reactions and a compiled basis
%
%    mode        - optional 'report' prints the reaction maps and closures
%
% Outputs:
%
%    K           - sparse first-order superoperator, or @(t,eta) returning
%                  a sparse generator for mass action or time-dependent
%                  rates. eta is a space-times-spin column; the handle
%                  assembles chemistry independently in each voxel.
%                  In zeeman-hilb with reactions, K(t,rho) instead returns
%                  the block-diagonal matrix derivative; only first-order
%                  reactions are supported. This is not a step generator.
%
% In Liouville space, K enters as 1i*K: L=H+1i*R+1i*K. Nonlinear propagation
% uses step(spin_system,{@(t,eta)1i*K(t,eta),t,'RKMK4'},eta,dt).
% Reaction maps are compiled once per kinetics call. Additive closure
% carries the unit once and each reactant's internal orders; product
% closure also carries cross-reactant orders. No concentration division
% is used. Spinless reactants are dynamic pools, not fixed reservoirs.
% Named selectors implement Haberkorn or Jones-Hore loss and projected
% product arrival. User selector pairs are local left/right product super-
% operators, including in zeeman-hilb (Liouville block dimensions). Zeeman
% maps use the complete spherical-tensor compiler and explicit dense local
% basis transformations; this reference path is intended for small systems.
% Each time-rate callback is evaluated once per generator evaluation;
% the resulting rate is shared by all spatial voxels at that stage time.
%
% ilya.kuprov@weizmann.ac.il
% ledwards@cbs.mpg.de
% hannah.hogben@chem.ox.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=kinetics.m>

function K=kinetics(spin_system,mode)

% Preserve the ordinary call and the explicitly requested reporting form
if nargin==1, mode=''; end
grumble(spin_system,mode);

% Leave chemistry-free systems available in every formalism
if isempty(spin_system.chem.reactions)
    K=mprealloc(spin_system,0); return
end

% Use the same reaction compiler in complete spherical-tensor coordinates
if ismember(spin_system.bas.formalism,{'zeeman-liouv','zeeman-hilb'})
    bas.formalism='sphten-liouv';
    bas.approximation=repmat({'none'},1,spin_system.bas.nsubst);
    source=basis(spin_system,bas); P=sphten2zeeman(source);
    for n=1:source.bas.nsubst
        idx=(source.bas.offsets(n)+1):source.bas.offsets(n+1);
        dim=prod(source.comp.mults(source.chem.parts{n}));
        P(idx,:)=P(idx,:)/dim;
    end
    Q=P\speye(size(P));

    % Express caller-supplied Zeeman product superoperators in the compiler basis
    for n=1:numel(source.chem.reactions)
        reaction=source.chem.reactions{n};
        if isfield(reaction,'selector')&&~ischar(reaction.selector{1})
            subst=reaction.reactants(1);
            idx=(source.bas.offsets(subst)+1):source.bas.offsets(subst+1);
            if ~isequal(size(reaction.selector{1}),[numel(idx) numel(idx)])
                error('Spinach:kinetics:selectorSize','selector matrices must have the reactant Liouville block dimensions.');
            end
            for k=1:2
                reaction.selector{k}=Q(idx,idx)*reaction.selector{k}*P(idx,idx);
            end
            source.chem.reactions{n}=reaction;
        end
    end

    % Transform fixed maps once and state-dependent maps at each stage state
    K=kinetics(source,mode);
    if isnumeric(K)
        K=P*K*Q;
    else
        zeeman=spin_system; zeeman.bas.formalism='zeeman-liouv';
        zeeman.bas.nstates=source.bas.nstates; zeeman.bas.offsets=source.bas.offsets;
        K=@(t,eta)zeeman_gen(zeeman,K,P,Q,t,eta);
    end
    if strcmp(spin_system.bas.formalism,'zeeman-hilb')
        K=@(t,rho)hilb_action(spin_system,K,t,rho);
    end
    return
end

% Compile immutable maps and selector products outside the time-step loop
reactions=spin_system.chem.reactions; dim=spin_system.bas.offsets(end);
for n=1:numel(reactions)
    reaction=reactions{n};
    reaction.maps=react_gen(spin_system,reaction);
    reaction.drain=cell(size(reaction.reactants));
    for k=1:numel(reaction.reactants)
        reactant=reaction.reactants(k);
        block=(spin_system.bas.offsets(reactant)+1):spin_system.bas.offsets(reactant+1);
        reaction.drain{k}=sparse(block,block,ones(size(block)),dim,dim);
    end
    if isfield(reaction,'selector')

        % Build substance-local singlet or triplet projector products
        reactant=reaction.reactants(1);
        block=(spin_system.bas.offsets(reactant)+1):spin_system.bas.offsets(reactant+1);
        unit=speye(numel(block)); selector=reaction.selector;
        if ischar(selector{1})
            electrons=num2cell(selector{2}); projectors=cell(1,2);
            for k=1:2
                sides={'left','right'};
                P=operator(spin_system,{'Lz','Lz'},electrons,sides{k});
                P=P+0.5*operator(spin_system,{'L+','L-'},electrons,sides{k});
                P=P+0.5*operator(spin_system,{'L-','L+'},electrons,sides{k});
                projectors{k}=unit/4-P(block,block);
                if endsWith(selector{1},'triplet'), projectors{k}=unit-projectors{k}; end
            end
        else
            projectors=selector;
            if ~isequal(size(projectors{1}),size(unit))
                error('Spinach:kinetics:selectorSize','selector matrices must have the reactant block dimensions.');
            end
        end

        % Separate selective loss from the projected arrival source
        drain=(projectors{1}+projectors{2})/2;
        if ischar(selector{1})&&startsWith(selector{1},'jones-hore')
            drain=unit-(unit-projectors{1})*(unit-projectors{2});
        end
        [rows,cols,values]=find(drain); offset=spin_system.bas.offsets(reactant);
        reaction.drain{1}=sparse(rows+offset,cols+offset,values,dim,dim);

        % Compile projected arrival into the full direct sum once
        reaction.fill=cell(size(reaction.maps));
        for p=1:numel(reaction.maps)
            rows=reaction.maps{p};
            fill=sparse(rows(:,1),rows(:,2)-offset,ones(size(rows,1),1),dim,numel(block));
            [rows,cols,values]=find(fill*projectors{1}*projectors{2});
            reaction.fill{p}=sparse(rows,cols+offset,values,dim,dim);
        end
    end
    reactions{n}=reaction;
end

% Print the declared network and its atom transport
if strcmp(mode,'report')
    summary_chemistry(spin_system);
    for n=1:numel(reactions)
        reaction=reactions{n};
        source_spins=cellfun(@(x)x(:)',spin_system.chem.parts(reaction.reactants),'UniformOutput',false);
        source_spins=[source_spins{:}];
        report(spin_system,['reaction ' num2str(n) ': matched spins ' ...
               mat2str(reaction.matching(:,1)') ', traced spins ' ...
               mat2str(setdiff(source_spins,reaction.matching(:,1))) ...
               ', closure ' reaction.closure]);
    end
end

% Numeric zero-rate records do not make first-order chemistry nonlinear
if all(cellfun(@(r)isnumeric(r.rate)&&(isscalar(r.reactants)||r.rate==0),reactions))
    K=assemble(spin_system,reactions,0,unit_state(spin_system));
else
    K=@(t,eta)assemble(spin_system,reactions,t,eta);
end

end

% Transform the shared reaction generator independently in every voxel
function K=zeeman_gen(spin_system,generator,P,Q,t,eta)

% Validate the physical state and determine its spatial multiplicity
concs=chem_concs(spin_system,eta); nvoxels=size(concs,1);
source=Q*reshape(eta,size(Q,1),nvoxels);

% Apply the direct-sum basis transformation on both sides of the map
K=generator(t,source(:));
K=kron(speye(nvoxels),P)*K*kron(speye(nvoxels),Q);

end

% Assemble polynomial drains and product-row fills at each stage state
function K=assemble(spin_system,reactions,t,eta)

% Read all local concentrations without normalising any state
concs=chem_concs(spin_system,eta);
dim=spin_system.bas.offsets(end);
eta=reshape(eta,dim,[]);
blocks=cell(size(concs,1),1);

% Evaluate each spatially uniform schedule once at the shared stage time
rates=zeros(size(reactions));
for n=1:numel(reactions)
    rate=reactions{n}.rate;
    if isa(rate,'function_handle'), rate=rate(t); end
    if ~isnumeric(rate)||~isscalar(rate)||~isreal(rate)||~isfinite(rate)||rate<0
        error('Spinach:kinetics:rateValue','time-dependent rates must return a finite non-negative scalar.');
    end
    rates(n)=rate;
end

% Reuse the stage rates in every voxel
for v=1:size(concs,1)
    K=sparse(dim,dim);
    for n=1:numel(reactions)
        reaction=reactions{n}; reactants=reaction.reactants;
        rate=rates(n);

        % Each occurrence loses its own state times the other concentrations
        for k=1:numel(reactants)
            others=[1:k-1 k+1:numel(reactants)];
            coeff=rate*prod(concs(v,reactants(others)));
            K=K-coeff*reaction.drain{k};
        end
        units=spin_system.bas.offsets(reactants)+1;
        for p=1:numel(reaction.maps)
            rows=reaction.maps{p};
            if isfield(reaction,'selector')

                % Trace the projected reactant into each product block
                K=K+rate*reaction.fill{p};
            elseif strcmp(reaction.closure,'product')

                % Gather all but the last factor into a state-dependent matrix
                coeff=rate*ones(size(rows,1),1);
                for k=1:numel(reactants)-1
                    coeff=coeff.*eta(rows(:,k+1),v);
                end
                K=K+sparse(rows(:,1),rows(:,end),coeff,dim,dim);
            else

                % Retain internal orders and share the unit among occurrences
                active=rows(:,2:end)~=units(:)';
                for k=1:numel(reactants)
                    others=[1:k-1 k+1:numel(reactants)];
                    present=~any(active(:,others),2);
                    coeff=rate*prod(concs(v,reactants(others)))*ones(nnz(present),1);
                    coeff(~any(active(present,:),2))=coeff(~any(active(present,:),2))/numel(reactants);
                    K=K+sparse(rows(present,1),rows(present,k+1),coeff,dim,dim);
                end
            end
        end
    end
    blocks{v}=K;
end

% Spatial transport is separate from the block-diagonal local chemistry
K=blkdiag(blocks{:});

end

% Consistency enforcement
function grumble(spin_system,mode)
for field={'rates','flux_rate','flux_type','rp_theory','rp_rates','rp_electrons'}
    if isfield(spin_system.chem,field{1})
        error('Spinach:kinetics:retiredChemistry',...
              ['chem.' field{1} ' is retired; use chem.reactions records.']);
    end
end
if isfield(spin_system.bas,'basis')&&~iscell(spin_system.bas.basis)
    error('Spinach:basis:retiredGlobalBasis',...
          'the global bas.basis matrix is retired; use bas.basis{n} and bas.offsets from basis().');
end
if isfield(spin_system.bas,'irrep')
    error('Spinach:basis:retiredIrrep',...
          'bas.irrep is retired; use bas.sym_fact(n).irr_projectors and irr_dimensions.');
end
if ~ischar(mode)||(~isempty(mode)&&~strcmp(mode,'report'))
    error('Spinach:kinetics:mode','the optional kinetics mode must be ''report''.');
end
if ~isempty(spin_system.chem.reactions)&&strcmp(spin_system.bas.formalism,'zeeman-wavef')
    error('Spinach:kinetics:wavefunction',...
          'chemical reactions are not supported in zeeman-wavef formalism.');
end
for n=1:numel(spin_system.chem.reactions)
    reaction=spin_system.chem.reactions{n};
    if strcmp(spin_system.bas.formalism,'zeeman-hilb')&&numel(reaction.reactants)>1
        error('Spinach:kinetics:hilbertMassAction',...
              'mass-action chemistry is not supported in zeeman-hilb formalism.');
    end
    if isfield(reaction,'selector')&&ischar(reaction.selector{1})&&...
       any(spin_system.comp.mults(reaction.selector{2})~=2)
        error('Spinach:kinetics:selectorSpin','named singlet/triplet selectors require two spin-half electrons.');
    end
end
end

% The egoist is the absolute sense is not the man who sacrifices others. He
% is the man who stands above the need of using others in any manner. He do-
% es not function through them. He is not concerned with them in any primary
% matter. Not in his aim, not in his motive, not in his thinking, not in his
% desires, not in the source of his energy. He does not exist for any other 
% man - and he asks no other man to exist for him. This is the only form of
% brotherhood and mutual respect possible between men.
%
% Ayn Rand, "The Fountainhead"

% #NGRUM

