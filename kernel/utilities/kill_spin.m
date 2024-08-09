% Removes the specified spins from the spin_system structure and
% updates it accordingly. Syntax:
%
%            spin_system=kill_spin(spin_system,hit_list)
%
% Parameters:
%
%     spin_system   - primary Spinach data structure
%
%     hit_list      - a vector of integers or a logical
%                     vector giving the numbers of spins
%                     to be removed from the system
%
% Outputs:
%
%     spin_system   - the data structure with the indica-
%                     ted spins and all dependent infor-
%                     mation (basis, assumptions) removed
%
% Notes: basis and assumption information is destroyed by this
%        function; you would need to call basis.m and assume.m
%        again.
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=kill_spin.m>

function spin_system=kill_spin(spin_system,hit_list)

% Check consistency
grumble(spin_system,hit_list)

% Catch logical indexing
if any(hit_list==0), hit_list=find(hit_list); end

% Inform the user
report(spin_system,['removing ' num2str(numel(hit_list)) ...
                    ' spins from the system...']);

% Update isotopes list and its hash
spin_system.comp.isotopes(hit_list)=[];
if ismember('op_cache',spin_system.sys.enable)
    spin_system.comp.iso_hash=md5_hash(spin_system.comp.isotopes);
end

% Update spin numbers
spin_system.comp.nspins=spin_system.comp.nspins-numel(hit_list);

% Update labels list
spin_system.comp.labels(hit_list)=[];

% Update multiplicities and magnetogyric ratios
spin_system.comp.mults(hit_list)=[];
spin_system.inter.gammas(hit_list)=[];

% Update base frequencies
spin_system.inter.basefrqs(hit_list)=[];

% Update Zeeman tensor array
spin_system.inter.zeeman.matrix(hit_list)=[];

% Update DD scaling multipliers
spin_system.inter.zeeman.ddscal(hit_list)=[];

% Update giant spin Hamiltonian terms
spin_system.inter.giant.coeff(hit_list)=[];

% Update coupling tensor array
spin_system.inter.coupling.matrix(hit_list,:)=[];
spin_system.inter.coupling.matrix(:,hit_list)=[];

% Update coordinates
spin_system.inter.coordinates(hit_list)=[];

% Update proximity matrix
spin_system.inter.proxmatrix(hit_list,:)=[];
spin_system.inter.proxmatrix(:,hit_list)=[];

% Update relaxation parameters
if ~isempty(spin_system.rlx.r1_rates)
    spin_system.rlx.r1_rates(hit_list)=[];
end
if ~isempty(spin_system.rlx.r2_rates)
    spin_system.rlx.r2_rates(hit_list)=[];
end
if ~isempty(spin_system.rlx.lind_r1_rates)
    spin_system.rlx.lind_r1_rates(hit_list)=[];
end
if ~isempty(spin_system.rlx.lind_r2_rates)
    spin_system.rlx.lind_r2_rates(hit_list)=[];
end
if ~isempty(spin_system.rlx.srfk_mdepth)
    spin_system.rlx.srfk_mdepth(hit_list,:)=[];
    spin_system.rlx.srfk_mdepth(:,hit_list)=[];
end
if ~isempty(spin_system.rlx.weiz_r1d)
    spin_system.rlx.weiz_r1d(hit_list,:)=[];
    spin_system.rlx.weiz_r1d(:,hit_list)=[];
end
if ~isempty(spin_system.rlx.weiz_r2d)
    spin_system.rlx.weiz_r2d(hit_list,:)=[];
    spin_system.rlx.weiz_r2d(:,hit_list)=[];
end

% Update kinetics parameters
for n=1:numel(spin_system.chem.parts)
    subsystem_idx=false(1,spin_system.comp.nspins);
    subsystem_idx(spin_system.chem.parts{n})=true();
    subsystem_idx(hit_list)=[];
    spin_system.chem.parts{n}=find(subsystem_idx);
end
if ~isempty(spin_system.chem.flux_rate)
    spin_system.chem.flux_rate(hit_list,:)=[];
    spin_system.chem.flux_rate(:,hit_list)=[];
end

% Update radical recombination parameters
reacting_spins=zeros(1,spin_system.comp.nspins+numel(hit_list));
reacting_spins(spin_system.chem.rp_electrons)=1; reacting_spins(hit_list)=[];
spin_system.chem.rp_electrons=find(reacting_spins);
if (~isempty(spin_system.chem.rp_rates))&&(numel(spin_system.chem.rp_electrons)<2)
    error('cannot destroy an essential electron in a radical pair system.');
end

% If any basis set information is found, destroy it
if isfield(spin_system,'bas')
    spin_system=rmfield(spin_system,'bas');
    report(spin_system,'WARNING - basis set information must be re-created.');
end

% If any assumption information is found, destroy it
if isfield(spin_system.inter.zeeman,'strength')
    spin_system.inter.zeeman=rmfield(spin_system.inter.zeeman,'strength');
    report(spin_system,'WARNING - assumption information must be re-created.');
end
if isfield(spin_system.inter.coupling,'strength')
    spin_system.inter.coupling=rmfield(spin_system.inter.coupling,'strength');
    report(spin_system,'WARNING - assumption information must be re-created.');
end

end

% Consistency enforcement
function grumble(spin_system,hit_list)
if islogical(hit_list)
    if(numel(hit_list)~=spin_system.comp.nspins)
        error('the size of the hit mask does not match the number of spins in the system.');
    end
else
    if (~isnumeric(hit_list))||any(hit_list<1)
        error('hit_list must be a logical matrix or an array of positive numbers.');
    end
    if any(hit_list>spin_system.comp.nspins)
        error('at least one number in hit_list exceeds the number of spins.');
    end
end
end

% I do not see that the sex of the candidate is an argument against
% her admission as privatdozent. After all, we are a university, not
% a bath house.
%
% David Hilbert, about Emmy Noether, in 1915.

