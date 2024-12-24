% Equilibrates linear chemical kinetics and returns a vector of
% equilibrium concentrations. Syntax:
%
%                      c=equilibrate(K,c0)
%
% Parameters:
%
%     K  - reaction rate matrix corresponding to 
%          dc/dt=K*c, where c is the concentration
%          vector
%
%     c0 - vector of initial concentrations
%
% Outputs:
%
%     c  - vector of equilibrium concentrations
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=equilibrate.m>

function c=equilibrate(K,c0)

% Check consistency
grumble(K,c0);

% Shortcut for zero concentrations
if norm(c0,2)==0, c=c0; return; end

% Recursive calls for independent reactions
indep_rxn_idx=scomponents(logical(K));
n_indep_rxns=numel(unique(indep_rxn_idx));
if n_indep_rxns>1
    c=zeros(size(c0));
    for n=1:n_indep_rxns
        P=eye(size(K)); P=P(:,indep_rxn_idx==n);
        c=c+P*equilibrate(P'*K*P,P'*c0);
    end
    return;
end

% Assemble the steady state system
A=vertcat(ones(1,size(K,2)),K);
b=vertcat(sum(c0),zeros(size(K,1),1));

% Check condition number
if cond(A)>1/sqrt(eps('double'))
    error('kinetic matrix does not constrain the equilibrium.');
end

% Solve the system
c=A\b;

end

% Consistency enforcement
function grumble(K,c0)
if (~isnumeric(K))||(~isreal(K))||(size(K,1)~=size(K,2))
    error('K must be a real square matrix.');
end
if ~all(abs(sum(K,1))<10*eps('double'))
    error('K violates conservation of matter: column sums must be zero.');
end
if (~isnumeric(c0))||(~isreal(c0))||(~iscolumn(c0))||any(c0<0)
    error('c0 must be a column vector of non-negative real numbers.');
end
if numel(c0)~=size(K,2)
    error('dimension mismatch between K and c0');
end
end

% In 2005, psychologists Belinda Board and Katarina Fritzon at the
% University of Surrey, UK, interviewed and gave personality tests
% to high-level British executives and compared their profiles with
% those of criminal psychiatric patients at Broadmoor Hospital in 
% the UK. They found that three out of eleven personality disorders
% were actually more common in executives than in the disturbed cri-
% minals. They were:
%
%    Histrionic personality disorder: including superficial charm, 
%    insincerity, egocentricity and manipulation.
%
%    Narcissistic personality disorder: including grandiosity,
%    self-focused lack of empathy for others, exploitativeness
%    and independence.
%
%    Obsessive-compulsive personality disorder: including perfec-
%    tionism, excessive devotion to work, rigidity, stubbornness
%    and dictatorial tendencies.
%
% They described these business people as successful psychopaths
% and the criminals as unsuccessful psychopaths.
%
% Wikipedia

