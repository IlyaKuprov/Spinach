% Computes the norm of the matrix represented by a tensor train. Syntax:
%
% #NORMOK               ttnorm=norm(ttrain,norm_type)
%
%       norm_type=1         returns the 1-norm
%       norm_type=inf       returns the inf-norm
%       norm_type=2         returns the 2-norm
%       norm_type='fro'     returns the Frobenius norm
%
% Note: norms other than Frobenius norm are expensive for tensor trains.
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

