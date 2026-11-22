% Evaluate the complete secular Hahn CCE-2 family without materialising
% its graph or cluster records. Internal clexpn crystallite evaluator.
%
% Syntax:
%
%   reduced=cce_analytical_hahn(spin_system,parameters,pulse_sequence,assumptions)
%
% Parameters:
%
%   spin_system    - create output with compact dipolar geometry
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
    model=spin_system.inter.coupling;
    R=euler2dcm(parameters.orientation);
    axis=R(3,:).';
    xyz=vertcat(spin_system.inter.coordinates{2:end})-spin_system.inter.coordinates{1};
    if isscalar(model.hyperfine)
        dist=vecnorm(xyz,2,2);
        hf=model.hyperfine*(1-3*((xyz*axis)./dist).^2)./dist.^3;
    else
        hf=(kron(axis,axis).'*reshape(model.hyperfine,9,[])).';
    end

    % Stream the complete pair family and retain only its product and screening bound
    numerics=parameters.cce.numerics;
    [reduced.by_order(:,:,3),omitted]=cce_pair_product(xyz,axis,hf,model.nuclear,...
                                                     spin_system.inter.tunnel_splittings,...
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
if ~isfield(spin_system.inter.coupling,'hyperfine')||...
   ~isequal(sequence,@hahn_echo)||~strcmp(assumptions,'labframe')||...
   ~isequal(p.observed,1)||~isequal(c.bath_spins,2:spin_system.comp.nspins)||...
   ~isfield(p,'transition')||~isfield(p.transition,'ms')||...
   numel(fieldnames(p.transition))~=2||~isequal(p.transition.ms,[1/2 -1/2])||...
   ~strcmp(c.approximation,'projected')
    error('clexpn:analyticalHahn','Compact dipolar geometry requires a projected electron Hahn echo with the full ordered bath.');
end
if c.order>2||~strcmp(c.graph.mode,'complete')||numel(fieldnames(c.graph))~=1||...
   ~strcmp(c.expansion,'standard')||~strcmp(c.bath_state.mode,'mixed_product')||...
   ~strcmp(c.meanfield,'none')||~strcmp(c.averaging,'after_expansion')||...
   c.second_order.enabled||c.retain_raw||~strcmp(c.backend,'auto')
    error('clexpn:analyticalHahn','Compact Hahn supports only the complete mixed-state secular expansion through order two.');
end
if any(~ismember(fieldnames(p),{'observed','times','orientation','transition','cce'}))||...
   any(isfield(c,{'basis','parallel','convergence','checkpoint','groups','partition',...
                  'bath_radius','interaction_cutoff','interlaced'}))||...
   isfield(c.numerics,'propagator')||~strcmp(c.reducer,'powers')
    error('clexpn:analyticalHahn','The compact ideal Hahn evaluator does not support the requested experiment controls.');
end
end


