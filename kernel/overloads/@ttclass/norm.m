% Computes the norm of the matrix represented by a tensor train. Syntax:
%
%                    ttnorm=norm(ttrain,norm_type)
%
% Parameters:
%
%    ttrain - a tensor train representation of a matrix
%
%    norm_type:
%
%       norm_type=1         not available for ttclass
%       norm_type=inf       not available for ttclass
%       norm_type=2         not available for ttclass
%       norm_type='fro'     returns the Frobenius norm
%
% Outputs:
%
%    ttnorm - a positive real number
%
% Note: only Frobenius norm is currently available for tensor trains;
%       other norm types raise errors.
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/norm.m>

function ttnorm=norm(ttrain,norm_type) %#NORMOK

% Compute the norm
switch norm_type
    
    case 'fro'
        
        % Frobenius norm
        ttrain=pack(ttrain); ttrain=ttort(ttrain,-1);
        ttnorm=abs(ttrain.coeff)*norm(ttrain.cores{1,1}(:),2);
        
    case 1
        
        % Maximum absolute column sum
        error('1-norm is not available for ttclass');
        
    case inf
        
        % Maximum absolute row sum
        error('inf-norm is not available for ttclass');
        
    case 2
        
        % Maximum absolute eigenvalue
        error('2-norm is not available for ttclass');
        
    otherwise
        
        % Complain and bomb out
        error('unrecognized norm type.');
        
end

end

% The most exciting phrase to hear in science, the one that heralds new
% discoveries, is not "eureka!" but rather "hmm....that's funny..."
%
% Isaac Asimov
%
%
% Contrary to what Asimov says, the most exciting phrase in science, the
% one that heralds new discoveries, is not "eureka!" or "that's funny...",
% it's "your research grant has been approved".
%
% John Alejandro King

