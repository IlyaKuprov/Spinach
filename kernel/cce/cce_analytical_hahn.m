% Evaluate the complete secular Hahn CCE-2 family without materialising
% its graph or cluster records. Internal clexpn crystallite evaluator.
%
% Syntax:
%
%   reduced=cce_analytical_hahn(spin_system,parameters,pulse_sequence,assumptions)
%
% Parameters:
%
%   spin_system    - create output with sys.tols.prox_cutoff=0; point
%                    dipoles are taken from inter.coordinates here
%   parameters     - resolved clexpn request for one orientation
%   pulse_sequence - registered @hahn_echo handle
%   assumptions    - 'labframe'
%
% Outputs:
%
%   reduced        - normalised Hahn echo, orders, and screening bounds
%
% The observed-only and singleton echoes are exactly one in the ideal-pulse
% limit, so the irreducible pair factors are the raw two-spin echoes.
%
% ilya.kuprov@weizmann.ac.il

function reduced=cce_analytical_hahn(spin_system,parameters,pulse_sequence,assumptions)

% Check that the requested experiment has the analytical pair limit
grumble(spin_system,parameters,pulse_sequence,assumptions);

% Construct the per-orientation output, including the trivial orders
order=parameters.cce.order;
shape=[numel(parameters.times) 1 order+1];
reduced.by_order=ones(shape);
reduced.by_order_error=zeros(shape);
reduced.by_order_valid=true(shape);
reduced.by_order_reason=repmat({'ok'},shape);
if order==2

    % Project central hyperfine interactions along the field in the molecular frame
    R=euler2dcm(parameters.orientation);
    axis=R(3,:).';
    xyz=vertcat(spin_system.inter.coordinates{2:end});
    bath=spin_system.comp.isotopes{2};
    gamma=spin(bath);
    scale=1e30*spin_system.tols.mu0/(4*pi);
    dipolar=scale*spin_system.tols.hbar*gamma^2;
    hf=zeros(size(xyz,1),1);
    if ~isempty(spin_system.inter.coordinates{1})
        xyz=xyz-spin_system.inter.coordinates{1};
        dist=vecnorm(xyz,2,2);
        g=spin_system.tols.freeg*spin_system.inter.zeeman.ddscal{1};
        hf=-scale*spin_system.tols.muB*gamma*((axis.'*g*axis)-...
           3*((xyz*g.'*axis).*(xyz*axis))./dist.^2)./dist.^3;
    end

    % Add the longitudinal part of user-supplied electron-nuclear tensors
    pairs=spin_system.inter.coupling.matrix;
    keys=sort(reshape(cell2mat(pairs(:,1:2)),[],2),2)-1;
    central=(keys(:,1)==0);
    zz=cellfun(@(tensor)axis.'*tensor*axis,pairs(central,3));
    hf=hf+accumarray(keys(central,2),zz,size(hf));

    % Isotropic bath exchange J adds J/2 to the flip-flop matrix element
    J=cellfun(@(tensor)tensor(1,1),pairs(~central,3));
    exchange=sparse(keys(~central,2),keys(~central,1),J,size(xyz,1),size(xyz,1));

    % Stream the complete pair family and retain only its product and screening bound
    numerics=parameters.cce.numerics;
    [reduced.by_order(:,:,3),omitted]=cce_pair_product(xyz,axis,hf,dipolar,...
                                                     exchange,...
                                                     parameters.times.',numerics.atol,...
                                                     numerics.max_memory_bytes);
    reduced.by_order_error(:,:,3)=omitted;
    reduced.by_order_error(parameters.times==0,1,3)=0;

end

% Match the standard reduced-output interface without raw cluster records
reduced.value=reduced.by_order(:,:,end);
reduced.abs_error_bound=reduced.by_order_error(:,:,end);
reduced.valid=reduced.by_order_valid(:,:,end);
reduced.reason=reduced.by_order_reason(:,:,end);
reduced.samples={};

end

% Enforce the physical and algorithmic domain of the approximation
function grumble(spin_system,p,sequence,assumptions)
c=p.cce;
comp=spin_system.comp;
if (comp.nspins<2)||~strcmp(comp.isotopes{1},'E')||~isnucleus(comp.isotopes{2})||...
   ~all(strcmp(comp.isotopes(2:end),comp.isotopes{2}))||(comp.mults(2)~=2)||...
   any(cellfun(@isempty,spin_system.inter.coordinates(2:end)))
    error('clexpn:analyticalHahn','Analytical Hahn pairs require one electron followed by a located bath of identical spin-1/2 nuclei.');
end
if ~isempty(spin_system.inter.coordinates{1})&&...
   any(all(vertcat(spin_system.inter.coordinates{2:end})==spin_system.inter.coordinates{1},2))
    error('clexpn:analyticalHahn','Analytical Hahn pairs require every bath nucleus to be away from the electron.');
end
if ~isfinite(spin_system.inter.magnet)||(spin_system.inter.magnet<=0)||...
   ~isempty(spin_system.rlx.theories)||~isscalar(spin_system.chem.parts)
    error('clexpn:analyticalHahn','Analytical Hahn pairs require a positive finite field, no relaxation theories, and one chemical species.');
end
shifts=spin_system.inter.zeeman.matrix(2:end);
if ~all(cellfun(@(z)isequal(z,shifts{1}),shifts))
    error('clexpn:analyticalHahn','Secular pairs require a common nuclear Zeeman tensor.');
end
keys=reshape(cell2mat(spin_system.inter.coupling.matrix(:,1:2)),[],2);
bath=spin_system.inter.coupling.matrix(all(keys>1,2),3);
if (spin_system.tols.prox_cutoff~=0)||any(keys(:,1)==keys(:,2))||...
   any(cellfun(@(tensor)~isequal(tensor,tensor(1,1)*eye(3)),bath))
    error('clexpn:analyticalHahn','Analytical Hahn pairs take point dipoles from coordinates: set sys.tols.prox_cutoff=0 and supply only hyperfine tensors and isotropic bath exchange.');
end
if ~isequal(sequence,@hahn_echo)||~strcmp(assumptions,'labframe')||...
   ~isequal(p.observed,1)||~isequal(c.bath_spins,2:spin_system.comp.nspins)||...
   ~isfield(p,'transition')||~isfield(p.transition,'ms')||...
   numel(fieldnames(p.transition))~=2||~isequal(p.transition.ms,[1/2 -1/2])||...
   ~strcmp(c.approximation,'projected')
    error('clexpn:analyticalHahn','Analytical Hahn pairs require a projected electron Hahn echo with the full ordered bath.');
end
if c.order>2||~strcmp(c.graph.mode,'complete')||numel(fieldnames(c.graph))~=1||...
   ~strcmp(c.expansion,'standard')||~strcmp(c.bath_state.mode,'mixed_product')||...
   ~strcmp(c.meanfield,'none')||~strcmp(c.averaging,'after_expansion')||...
   c.second_order.enabled||c.retain_raw
    error('clexpn:analyticalHahn','Analytical Hahn pairs support only the complete mixed-state secular expansion through order two.');
end
if any(~ismember(fieldnames(p),{'observed','times','orientation','transition','cce'}))||...
   any(isfield(c,{'basis','parallel','convergence','checkpoint','groups','partition',...
                  'bath_radius','interaction_cutoff','interlaced'}))||...
   isfield(c.numerics,'propagator')||~strcmp(c.reducer,'powers')
    error('clexpn:analyticalHahn','The analytical ideal Hahn evaluator does not support the requested experiment controls.');
end
end


