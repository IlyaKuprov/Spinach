% Computes trajectory similarity scores. Returns a function representing 
% "similarity" of the two state space trajectories at different points in
% time. See http://dx.doi.org/10.1016/j.jmr.2013.02.012 for further infor-
% mation. Syntax:
%
%         trajsimil(spin_system,trajectory_1,trajectory_2,method)
%
% Parameters:
%
%   trajectory_1,2 - spin system trajectories, supplied as nstates
%                    x nsteps matrices.
%
%   method         - similarity scoring method; possibilities are:
%
%                     'RSP'  - running scalar product. Computes
%                              scalar products between the cor-
%                              responding vectors of the trajec-
%                              tories.
%
%                     'RDN'  - running difference norm. The two
%                              trajectories are subtracted and 
%                              difference 2-norms returned.
%
%                     'SG-'  - prefix that turns on state grou-
%                              ping. T(l,m) and T(l,-m) states
%                              of each spin (standalone or in 
%                              direct products with other ope-
%                              rators)will be considered equva-
%                              lent.
%
%                     'BSG-' - prefix that turns on broad state
%                              grouping. All states of a given
%                              spin (standalone or in direct 
%                              products with other operators)
%                              will be considered equivalent.
%
%                    The possible combinations are: 'RSP','RDN',
%                    'SG-RSP','SG-RDN','BSG-RSP','BSG-RDN'.
%
% State grouping consists in summing the absolute squares of the coeffi-
% cients to be grouped and taking the square root. The trajectories would
% usually come out of the evolution.m or krylov.m run from a given star-
% ting point under a given Liouvillian.
%
% Output:
%
%     this function writes into the current figure
%
% Note: SG and BSG options require sphten-liouv formalism.
%
% ilya.kuprov@weizmann.ac.il
% kpervushin@ntu.edu.sg
%
% <https://spindynamics.org/wiki/index.php?title=trajsimil.m>

function score=trajsimil(spin_system,trajectory_1,trajectory_2,scorefcn)

% Check consistency
grumble(spin_system,trajectory_1,trajectory_2,scorefcn);

% Run state grouping
if strcmp(scorefcn(1:3),'SG-')||strcmp(scorefcn(1:3),'BSG')
    
    % Grab the description of the current basis
    state_list=spin_system.bas.basis;
    
    % Tell the user we're started
    report(spin_system,'collapsing equivalent subspaces...');
    
    % Decide how to proceed
    if strcmp(scorefcn(1:3),'SG-')
        
        % Rename all T(l,-m) states into T(l,m) states
        [L,M]=lin2lm(state_list); state_list=lm2lin(L,abs(M));
        
        % Update the method variable
        scorefcn=scorefcn(4:end);
        
    elseif strcmp(scorefcn(1:3),'BSG')
        
        % Rename all non-identity states into Lz
        state_list(state_list~=0)=2;
        
        % Update the method variable
        scorefcn=scorefcn(5:end);
        
    end
    
    % Index all unique and repeated states on the modified state list
    [grouped_state_list,index_forward,index_backward]=unique(state_list,'rows');
    
    % Preallocate state-grouped trajectories
    grouped_trajectory_1=zeros(length(index_forward),size(trajectory_1,2));
    grouped_trajectory_2=zeros(length(index_forward),size(trajectory_2,2));
    
    % Group trajectory tracks corresponding to the states that are
    % flagged as identical in the indices (root-sum-square)
    for n=1:length(index_forward)
        grouped_trajectory_1(n,:)=sqrt(sum(abs(trajectory_1(index_backward==n,:)).^2,1));
        grouped_trajectory_2(n,:)=sqrt(sum(abs(trajectory_2(index_backward==n,:)).^2,1));
    end
    
    % Update the variables
    trajectory_1=grouped_trajectory_1;
    trajectory_2=grouped_trajectory_2;
    
    % Tell the user we're done
    report(spin_system,[num2str(size(state_list,1)) ' states collected into ' num2str(size(grouped_state_list,1)) ' groups.']);
    
end

% Preallocate the result array
trajectory_length=size(trajectory_1,2);
score=zeros(1,trajectory_length);

% Compute the score function
switch scorefcn
    
    case 'RSP'
        
        % Compute running scalar product scores
        for n=1:trajectory_length
            score(n)=dot(trajectory_1(:,n),trajectory_2(:,n));
        end
        
    case 'RDN'
        
        % Compute running difference norm scores
        for n=1:trajectory_length
            score(n)=1-norm(trajectory_1(:,n)-trajectory_2(:,n),2)/2;
        end
       
    otherwise
        
        % Complain and bomb out
        error('unknown similarity score function.');
        
end

end

% Consistency enforcement 
function grumble(spin_system,trajectory_1,trajectory_2,scorefcn)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('this function is only available for Liouville space formalisms.');
end
if (strcmp(scorefcn(1:3),'SG-')||strcmp(scorefcn(1:3),'BSG'))&&...
   (~ismember(spin_system.bas.formalism,{'sphten-liouv'}))
    error('state grouping is only available for sphten-liouv formalism.');
end
if any(size(trajectory_1)~=size(trajectory_2))
    error('matrix dimensions of the two trajectories should match.');
end
if (size(trajectory_1,1)~=size(spin_system.bas.basis,1))||...
   (size(trajectory_2,1)~=size(spin_system.bas.basis,1))
    error('trajectory dimension should be equal to the basis set dimension.');
end
if (~isnumeric(trajectory_1))||(~isnumeric(trajectory_2))
    error('both trajectories should be arrays of doubles.');
end
if ~ismember(scorefcn,{'RSP','RDN','SG-RSP','SG-RDN','BSG-RSP','BSG-RDN'})
    error('available score functions are RSP, RDN, SG-RSP, SG-RDN, BSG-RSP, and BSG-RDN.');
end
end

% How would we look for a new law? In general we look for a new law by the
% following process. First, we guess it. Then we... don't laugh. That's the
% damned truth. Then we compute the consequences of the guess... to see if
% this is right, to see if this law we guessed is right, to see what it
% would imply. And then we compare those computation results to nature. Or
% we say to compare it to experiment, or to experience. Compare it directly
% with observations to see if it works. If it disagrees with experiment,
% it's wrong. In that simple statement is the key to science. It doesn't
% make a difference how beautiful your guess is. It doesn't make a diffe-
% rence how smart you are, who made the guess or what his name is... If it
% disagrees with experiment, it's wrong. That's all there is to it.
%
% Richard Feynman

