% Performs TT-orthogonalisation for a tensor train (or for each tensor 
% train in a buffered sum). Syntax:
%
%                  [tt,lognrm]=ttort(tt,direct)
%
% Parameters:
%
%    direct=+1 - {default} gives you left-to-right orthogonality,
%    direct=-1 - gives right-to-left orthogonality
%
%    tt     - tensor train object, possibly with buffered sums
%
% Outputs:
%
%    tt     - tensor train object with all terms in the buffe-
%             red sum has all of them orthogonalised in the 
%             direction requested
%
%    lognrm - if this output is present, all buffered trains 
%             are also normalized, and natural logs of their
%             norms returned in the vector lognrm. Use this
%             option if the tensor norm is likely to exceed
%             realmax()=1.7977e+308.
%
% Note: normally you should not call this subroutine directly.
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/ttort.m>

function [tt,lognrm]=ttort(tt,direct)

% Read tensor ranks and dimensions
sz=tt.sizes; rnk=tt.ranks;
d=tt.ncores; N=tt.ntrains;

% Set the default dir if it is not present
if ~exist('dir','var'), direct=+1; end

% Compute logs of norms
if (nargout>1)
    lognrm=log(tt.coeff);
    tt.coeff=ones(1,N);
end

% Decide the direction
switch direct
    
    % Left-to-right direction
    case +1
        
        % Loop over buffered trains
        for n=1:N
            
            % Get the ranks of the current train
            r=rnk(:,n);
            
            % TT-othogonalise left-to-right
            for k=1:d-1
                
                % Shape cores into matrices
                B=reshape(tt.cores{k,n}, [r(k)*sz(k,1)*sz(k,2), r(k+1)]);
                C=reshape(tt.cores{k+1,n}, [r(k+1), sz(k+1,1)*sz(k+1,2)*r(k+2)]);
                
                % Run orthogonal-triangular decomposition
                [Q,R]=qr(B,0); nrm=norm(R,2);
                
                % Take care of the norm
                if (nrm>0)
                    R=R/nrm;
                    if (nargout>1)
                        lognrm(1,n)=lognrm(1,n)+log(nrm);
                    else
                        tt.coeff(1,n)=tt.coeff(1,n)*nrm;
                    end
                end
                
                % Multiply triangular part forward
                C=R*C; rnew=size(Q,2);
                
                % Fold the cores back
                tt.cores{k,n}=reshape(Q, [r(k), sz(k,1), sz(k,2), rnew]);
                tt.cores{k+1,n}=reshape(C, [rnew, sz(k+1,1), sz(k+1,2), r(k+2)]);
                
                % Update the ranks
                r(k+1)=rnew;
                
            end
            
            % Take care of the last norm
            nrm=norm(tt.cores{d,n}(:),2);
            if (nrm>0)
                tt.cores{d,n}=tt.cores{d,n}/nrm;
                if (nargout>1)
                    lognrm(1,n)=lognrm(1,n)+log(nrm);
                else
                    tt.coeff(1,n)=tt.coeff(1,n)*nrm;
                end
            end
            
        end
        
    % Right-to-left direction
    case -1
        
        % Loop over buffered trains
        for n=1:N
            
            % Get the ranks of the current train
            r=rnk(:,n);
            
            % TT-othogonalise right-to-left
            for k=d:-1:2
                
                % Shape cores into matrices
                B=reshape(tt.cores{k,n}, [r(k), sz(k,1)*sz(k,2)*r(k+1)]);
                C=reshape(tt.cores{k-1,n}, [r(k-1)*sz(k-1,1)*sz(k-1,2), r(k)]);
                
                % Run orthogonal-triangular decomposition
                [Q,R]=qr(B.',0); nrm=norm(R,2);
                
                % Take care of the norm
                if (nrm>0)
                    R=R/nrm;
                    if (nargout>1)
                        lognrm(1,n)=lognrm(1,n)+log(nrm);
                    else
                        tt.coeff(1,n)=tt.coeff(1,n)*nrm;
                    end
                end
                
                % Multiply triangular part backward
                C=C*R.'; rnew=size(Q,2);
                
                % Fold the cores back
                tt.cores{k,n}=reshape(Q.', [rnew, sz(k,1), sz(k,2), r(k+1)]);
                tt.cores{k-1,n}=reshape(C, [r(k-1), sz(k-1,1), sz(k-1,2), rnew]);
                
                % Update the ranks
                r(k)=rnew;
                
            end
            
            % Take care of the last norm
            nrm=norm(tt.cores{1,n}(:),2);
            if (nrm>0)
                tt.cores{1,n}=tt.cores{1,n}/nrm;
                if (nargout>1)
                    lognrm(1,n)=lognrm(1,n)+log(nrm);
                else
                    tt.coeff(1,n)=tt.coeff(1,n)*nrm;
                end
            end
            
        end
        
    otherwise
        
        % Complain and bomb out
        error('unrecognized direction.');
        
end

end

% The order of authorship was determined from a twenty-five-game 
% croquet series held at Imperial College Field Station during 
% summer 1973.
%
% A footnote in http://dx.doi.org/10.2307/3384

