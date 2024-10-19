% Plots grid integration quality as a function of spherical rank. The
% quality is defined as the norm of the residual of spherical harmon-
% ics or Wigner functions integrated using the grid provided. Syntax:
%
%  grid_profile=grid_test(alphas,betas,gammas,weights,max_rank,sfun)
%
% Parameters:
%
%      alphas - alpha Euler angles of the grid, in radians,
%               zeros for single-angle grids
%
%       betas - beta Euler angles of the grid, in radians
%
%      gammas - gamma Euler angles of the grid, in radians,
%               zeros for two-angle grids
%
%     weights - point weights of the grid
%
%       ranks - spherical ranks to consider
%
%        sfun - spherical function type: for three-angle
%               grids use 'D_lmn', for two-angle grids use
%               'Y_lm', for single-angle grids use 'Y_l0'.
%
% The output is a vector of residual norms in each spherical rank.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=grid_test.m>

function grid_profile=grid_test(alphas,betas,gammas,weights,ranks,sfun)

% Check consistency
grumble(alphas,betas,gammas,weights,ranks,sfun);

% Preallocate the answer
grid_profile=zeros(size(ranks));

% Loop over spherical ranks
for k=1:numel(ranks)
    
    % Preallocate Wigner matrix
    D=zeros(2*ranks(k)+1,'like',1i);
    
    % Loop over grid points
    parfor n=1:numel(alphas)
        D=D+weights(n)*wigner(ranks(k),alphas(n),...
                              betas(n),gammas(n)); %#ok<PFBNS>
    end
    
    % Update grid profile
    if strcmp(sfun,'D_lmn')
        grid_profile(k)=norm(D,2)-krondelta(0,ranks(k));
    elseif strcmp(sfun,'Y_lm')
        grid_profile(k)=norm(D(ranks(k)+1,:),2)-krondelta(0,ranks(k));
    elseif strcmp(sfun,'Y_l0')
        grid_profile(k)=norm(D(ranks(k)+1,ranks(k)+1),2)-krondelta(0,ranks(k));
    else
        error('unknown diagnostics function.');
    end
    
    % Update the user
    disp(['Spherical rank ' num2str(ranks(k)) ...
          ', residual ' sfun ' norm: ' num2str(grid_profile(k))]);
    
end
    
% Do the plotting
if nargout==0
    figure(); plot(ranks,grid_profile);
    kxlabel('spherical rank'); kgrid;
    kylabel('integration residual');
end

end

% Consistency enforcement
function grumble(alphas,betas,gammas,weights,ranks,sfun)
if (~isnumeric(alphas))||(~isreal(alphas))||...
   any(~isfinite(alphas))||(size(alphas,2)~=1)
    error('alphas must be a column vector of real numbers.');
end
if (~isnumeric(betas))||(~isreal(betas))||...
   any(~isfinite(betas))||(size(betas,2)~=1)
    error('betas must be a column vector of real numbers.');
end
if (~isnumeric(gammas))||(~isreal(gammas))||...
   any(~isfinite(gammas))||(size(gammas,2)~=1)
    error('gammas must be a column vector of real numbers.');
end
if (~isnumeric(weights))||(~isreal(weights))||...
   any(~isfinite(weights))||(size(weights,2)~=1)||any(weights<=0)
    error('weights must be a column vector of positive real numbers.');
end
if numel(alphas)~=numel(betas)
    error('the number of elements in alphas and betas must be the same.');
end
if numel(betas)~=numel(gammas)
    error('the number of elements in betas and gammas must be the same.');
end
if numel(gammas)~=numel(weights)
    error('the number of elements in gammas and weights must be the same.');
end
if (~isnumeric(ranks))||(~isreal(ranks))||(~isvector(ranks))||...
   any(~isfinite(ranks))||any(ranks<0)||any(mod(ranks,1))
    error('ranks must be a vector of non-negative integers.');
end
if ~ismember(sfun,{'D_lmn','Y_lm','Y_l0'})
    error('sfun argument can be ''D_lmn'', ''Y_lm'' or ''Y_l0''');
end
end

% Art is what you can get away with.
%
% Andy Warhol

