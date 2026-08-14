% Unbiased stochastic compression of a state vector. Reduces the
% number of non-zeroes to a user-specified target in such a way
% that the expectation value of every element of the result is
% exactly equal to the corresponding element of the input. Ele-
% ments with magnitudes above a calibrated threshold are kept
% exactly, and the remainder are retained with probabilities pro-
% portional to their magnitudes using systematic resampling, the
% survivors being rescaled by the reciprocal of their retention
% probability. In contrast to the deterministic amplitude thres-
% holding performed by clean_up.m, the error introduced here has
% zero mean, and therefore averages out across an ensemble of in-
% dependent trajectories, allowing the truncation error to be es-
% timated from the scatter of the ensemble. Syntax:
%
%                   v=stoch_chop(v,nnz_target)
%
% Parameters:
%
%    v          - a state vector, a complex column vector
%
%    nnz_target - the number of non-zeroes to keep, a posi-
%                 tive integer; an input having this many
%                 non-zeroes or fewer is returned unchanged
%
% Outputs:
%
%    v          - a column vector of the same dimension and
%                 storage class as the input, having at most
%                 nnz_target non-zeroes
%
% Note: the threshold is the unique solution of the requirement
%       that the retention probabilities of the sub-threshold
%       elements must sum to the number of slots left over by
%       the elements kept exactly, and must not exceed unity.
%
% Note: the random number stream is not reset, call rng() before
%       this function when a reproducible sequence is required.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=stoch_chop.m>

function v=stoch_chop(v,nnz_target)

% Check consistency
grumble(v,nnz_target);

% Locate the non-zeroes
nz_idx=find(v);

% Bail out when the vector is already sparse enough
if numel(nz_idx)<=nnz_target, return; end

% Rank the non-zeroes by decreasing magnitude
[mags,sort_idx]=sort(abs(nonzeros(v)),'descend');
ranked_idx=nz_idx(sort_idx);

% Sums of the magnitudes left over by each retention count
tail_sums=sum(mags)-[0; cumsum(mags(1:(nnz_target-1)))];

% Thresholds implied by each retention count
thresholds=tail_sums./((nnz_target:-1:1)');

% Smallest retention count at which no probability exceeds unity
n_kept=find(mags(1:nnz_target)<=thresholds,1)-1;

% Retention probabilities of the sub-threshold elements
retn_prob=mags((n_kept+1):end)/thresholds(n_kept+1);

% Cumulative probability axis, closed at the exact total
bin_edges=[0; cumsum(retn_prob)];
bin_edges(end)=nnz_target-n_kept;

% Systematically placed sample points, one per available slot
sample_pts=rand()+(0:(nnz_target-n_kept-1))';

% Sub-threshold elements selected by the sample points
picked=discretize(sample_pts,bin_edges);

% Rescale the survivors to keep the expectation value exact
v(ranked_idx(n_kept+picked))=v(ranked_idx(n_kept+picked))./retn_prob(picked);

% Discard everything that was neither kept nor picked
dropped=true(numel(nz_idx),1);
dropped(1:n_kept)=false;
dropped(n_kept+picked)=false;
v(ranked_idx(dropped))=0;

end

% Consistency enforcement
function grumble(v,nnz_target)
if (~isnumeric(v))||(~ismatrix(v))||(size(v,2)~=1)
    error('v must be a numeric column vector.');
end
if any(~isfinite(nonzeros(v)))
    error('v must not have NaN or Inf elements.');
end
if (~isnumeric(nnz_target))||(~isreal(nnz_target))||(~isscalar(nnz_target))||...
   (mod(nnz_target,1)~=0)||(nnz_target<1)
    error('nnz_target must be a positive real integer.');
end
end

% Anyone who considers arithmetical methods of producing random
% digits is, of course, in a state of sin.
%
% John von Neumann

