% Computes SHREWD weights for a given two- or three-angle spherical
% grid. See the paper by Eden and Levitt for details on now the al-
% gorithm works: http://dx.doi.org/10.1006/jmre.1998.1427 Syntax:
%
%      weights=shrewd(alphas,betas,gammas,max_rank,max_error)
%
% Parameters:
%
%      alphas - alpha Euler angles (ZYZ active) of the 
%               grid, in radians, set to all-zeros for
%               two-angle grids
%
%       betas - beta Euler angles (ZYZ active) of the 
%               grid, in radians
%
%      gammas - gamma Euler angles (ZYZ active) of the
%               grid, in radians
%
%    max_rank - maximum spherical rank to take into consi-
%               deration when minimizing residuals
%
%   max_error - maximum residual absolute error per spheri-
%               cal function
%
% Outputs
%
%     weights - a vector of grid weights for each 
%               [alpha beta gamma] point supplied.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=shrewd.m>

function weights=shrewd(alphas,betas,gammas,max_rank,max_error)

% Check consistency
grumble(alphas,betas,gammas,max_rank,max_error);

% Decide the grid type
if all(alphas==0)

    % Preallocate spherical harmonic matrix
    H=zeros(lm2lin(max_rank,-max_rank)+1,numel(alphas),'like',1i);

    % Fill spherical harmonic matrix
    for k=1:numel(alphas)
        for l=0:max_rank
            D=wigner(l,alphas(k),betas(k),gammas(k));
            for m=l:-1:-l
                H(lm2lin(l,m)+1,k)=D(l+1,l+m+1);
            end
        end
    end

    % Get the right hand side vector
    v=max_error*ones(lm2lin(max_rank,-max_rank)+1,1);
    v(1)=1-max_error;
    
else
    
    % Preallocate Wigner function matrix
    H=zeros(lmn2lin(max_rank,-max_rank,-max_rank),numel(alphas));

    % Fill Wigner function matrix
    for k=1:numel(alphas)
        for l=0:max_rank
            D=wigner(l,alphas(k),betas(k),gammas(k));
            for m=l:-1:-l
                for n=l:-1:-l
                    H(lmn2lin(l,m,n),k)=D(l+m+1,l+n+1);
                end
            end
        end
    end

    % Get the right hand side vector
    v=max_error*ones(lmn2lin(max_rank,-max_rank,-max_rank),1);
    v(1)=1-max_error;
    
end

% Compute the weights
weights=real(H\v); weights=weights/sum(weights);

% Run some diagnostics
if any(weights==0)
    error('zero weights detected - increase the maximum rank.');
end
if any(weights<0)
    error('negative weights detected - reduce the accuracy threshold.');
end

end

% Consistency enforcement
function grumble(alphas,betas,gammas,max_rank,max_error)
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
if numel(alphas)~=numel(betas)
    error('the number of elements in alphas and betas must be the same.');
end
if numel(betas)~=numel(gammas)
    error('the number of elements in betas and gammas must be the same.');
end
if (~isnumeric(max_rank))||(~isreal(max_rank))||(~isfinite(max_rank))||...
   (numel(max_rank)~=1)||(max_rank<1)||mod(max_rank,1)
    error('max_rank must be a positive real integer.');
end
if (~isnumeric(max_error))||(~isreal(max_error))||(~isfinite(max_error))||...
   (numel(max_error)~=1)||(max_error<0)
    error('max_error must be a non-negative real number.');
end
end

% A good theory explains a lot, but postulates little.
%
% Richard Dawkins

