% Combinatorial PCS fitting function. Takes into account potential am-
% biguities in diamagnetic and paramagnetic NMR assignments. Syntax:
%
%     [d_shifts,p_shifts,pcs_theo,...
%               pcs_expt,chi,total_theo]=pcs_combi_fit(parameters)
%
% Parameters:
%
%     parameters.hfcs   - cell array of 3x3 hyperfine tensors,
%                         in Gauss, usually out of gparse() or
%                         something similar
%
%     parameters.nel    - number of unpaired electrons involved
%
%     parameters.isotopes - cell array of isotope specificati-
%                           ons, e.g. {'1H','1H'}
%
%     parameters.spin_groups - cell array of integer vectors
%                              specifying the numbers of spins
%                              that have each of the chemical 
%                              shifts specified, e.g. 
%                              {[28 22 30]; [25 33 81]}
%
%     parameters.d_shifts - a vector of unique diamagnetic che-
%                           mical shifts, in ppm
%
%     parameters.p_shifts - a vector of unique paramagnetic 
%                           chemical shifts, in ppm
%
%     parameters.d_ambig - a cell array of integer vectors spe-
%                          cifying the spins for which the dia-
%                          magnetic assignment can potentially
%                          be swapped around.
%
%     parameters.p_ambig - a cell array of integer vectors spe-
%                          cifying the spins for which the para-
%                          magnetic assignment can potentially
%                          be swapped around.
%
% Outputs:
%
%     d_shifts - diamagnetic chemical shifts, optimally permuted
%
%     p_shifts - paramagnetic chemical shifts, optimally permuted
%
%     pcs_theo - theoretical pseudocontact shifts
%
%     pcs_expt - experimental pseudocontact shifts from optimally
%                permuted assignments
%
%     chi      - rank 2 part of the magnetic susceptibility tensor,
%                in cubic Angstrom
%
%     total_theo - theoretical total NMR chemcial shifts, computed
%                  as a sum of d_shifts and pcs_theo
%
% ilya.kuprov@weizmann.ac.il
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=pcs_combi_fit.m>

function [d_shifts,p_shifts,pcs_theo,pcs_expt,chi,total_theo]=pcs_combi_fit(parameters)

% Check consistency
grumble(parameters);

% Build diamagnetic shift permutation table
perm_table_d=perms(parameters.d_ambig{1});
for n=2:numel(parameters.d_ambig)
    local_perms=perms(parameters.d_ambig{n});
    perm_table_d=[kron(perm_table_d,ones(size(local_perms,1),1)) kron(ones(size(perm_table_d,1),1),local_perms)];
end
disp([num2str(size(perm_table_d,1)) ' diamagnetic shift permutations']);
diamag_shift_perms=zeros(size(perm_table_d,1),numel(parameters.d_shifts));
for n=1:size(perm_table_d,1)
    diamag_shift_perms(n,:)=parameters.d_shifts;
    diamag_shift_perms(n,cell2mat(parameters.d_ambig))=diamag_shift_perms(n,perm_table_d(n,:));
end

% Build paramagnetic shift permutation table
perm_table_p=perms(parameters.p_ambig{1});
for n=2:numel(parameters.p_ambig)
    local_perms=perms(parameters.p_ambig{n});
    perm_table_p=[kron(perm_table_p,ones(size(local_perms,1),1)) kron(ones(size(perm_table_p,1),1),local_perms)];
end
disp([num2str(size(perm_table_p,1)) ' paramagnetic shift permutations']);
paramag_shift_perms=zeros(size(perm_table_p,1),numel(parameters.p_shifts));
for n=1:size(perm_table_p,1)
    paramag_shift_perms(n,:)=parameters.p_shifts;
    paramag_shift_perms(n,cell2mat(parameters.p_ambig))=paramag_shift_perms(n,perm_table_p(n,:));
end

% Compute permutation scores
N=size(diamag_shift_perms,1); M=size(paramag_shift_perms,1); score=zeros(N*M,1);
disp(['Evaluating a total of ' num2str(N*M) ' permutations...']);
parfor n=1:N*M %#ok<*PFBNS>
    
    % Disentangle indices
    [d_index,p_index]=ind2sub([N M],n);
    
    % Compute pseudocontact shifts
    pcs=paramag_shift_perms(p_index,:)-diamag_shift_perms(d_index,:); 
    
    % Build the input arrays
    hfcs=cell(numel(cell2mat(parameters.spin_groups)),1); p=1;
    pcs_expt=zeros(numel(cell2mat(parameters.spin_groups)),1); 
    for k=1:numel(parameters.spin_groups)
        for m=1:numel(parameters.spin_groups{k})
            hfcs{p}=parameters.hfcs{parameters.spin_groups{k}(m)};
            pcs_expt(p)=pcs(k); p=p+1;
        end
    end
    
    % Get the assignment score
    [~,score(n)]=pcs2chi(hfcs,pcs_expt,parameters.isotopes,parameters.nel);
    
end

% Find the shift permutations that led to the best score
[min_score,n]=min(score);
[d_index,p_index]=ind2sub([N M],n);
d_shifts=diamag_shift_perms(d_index,:);
p_shifts=paramag_shift_perms(p_index,:);
disp(['Minimum score: ' num2str(min_score)]);

% Build the input arrays
hfcs=cell(numel(cell2mat(parameters.spin_groups)),1); p=1;
pcs_expt=zeros(numel(cell2mat(parameters.spin_groups)),1);
for k=1:numel(parameters.spin_groups)
    for m=1:numel(parameters.spin_groups{k})
        hfcs{p}=parameters.hfcs{parameters.spin_groups{k}(m)};
        pcs_expt(p)=p_shifts(k)-d_shifts(k); p=p+1;
    end
end

% Recover the susceptibility
chi=pcs2chi(hfcs,pcs_expt,parameters.isotopes,parameters.nel);

% Compute the predicted pseudocontact shifts
pcs_theo=zeros(numel(hfcs),1);
for k=1:numel(hfcs)
    pcs_theo(k)=hfc2pcs(hfcs{k},chi,parameters.isotopes{k},parameters.nel);
end

% Compute the predicted total shifts
total_theo=pcs_theo+kron(d_shifts',[1; 1; 1]);

% Report to the user
disp('Best permuted diamagnetic shifts:'); disp(d_shifts');
disp('Best permuted paramagnetic shifts:'); disp(p_shifts');
disp('PCS (theory, experiment):'); disp([pcs_theo pcs_expt]);
disp('Susceptibility tensor anisotropy (cubic Angstroms):'); disp(chi);
disp('Expected experimental chemical shifts:'); disp(total_theo);

end

% Consistency enforcement
function grumble(parameters)
if ~isfield(parameters,'hfcs')
    error('parameters.hfcs field is missing.');
end
if ~isfield(parameters,'isotopes')
    error('parameters.isotopes field is missing.');
end
if ~isfield(parameters,'spin_groups')
    error('parameters.spin_groups field is missing.');
end
if ~isfield(parameters,'d_shifts')
    error('parameters.d_shifts field is missing.');
end
if ~isfield(parameters,'p_shifts')
    error('parameters.p_shifts field is missing.');
end
if ~isfield(parameters,'d_ambig')
    error('parameters.d_ambig field is missing.');
end
if ~isfield(parameters,'p_ambig')
    error('parameters.p_ambig field is missing.');
end
end

% Shall I refuse my dinner because I do not fully understand the 
% process of digestion?
%
% Oliver Heaviside

