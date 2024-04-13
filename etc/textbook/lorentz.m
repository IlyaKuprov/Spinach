% The (L,0)(+)(0,L) irreducible matrix representation of the
% Lorentz group with inversion. Syntax:
%
%                    [J,K,Kil]=lorentz(L)
%
% Parameters:
%
%     L    -   irreducible representation 
%              rank, e.g. 1/2
%
% Outputs:
%
%     J    -   three rotation generators
%
%     K    -   three boost generators
%
%     Kil  -   Killing form (expensive)
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=lorentz.m>

function [J,K,Kil]=lorentz(L)

% Check consistency
grumble(L);

% Dimension and Pauli blocks
D=2*L+1; s=pauli(D);

% Rotation and boost generators
J.x=full([s.x zeros(D,D); zeros(D,D) s.x]);
J.y=full([s.y zeros(D,D); zeros(D,D) s.y]);
J.z=full([s.z zeros(D,D); zeros(D,D) s.z]);
K.x=full([1i*s.x zeros(D,D); zeros(D,D) -1i*s.x]);
K.y=full([1i*s.y zeros(D,D); zeros(D,D) -1i*s.y]);
K.z=full([1i*s.z zeros(D,D); zeros(D,D) -1i*s.z]);

% Collect generators
gens={J.x J.y J.z K.x K.y K.z};

% Killing form is expensive
if nargout>2
    
    % Preallocate Killing form
    Kil=zeros(numel(gens),numel(gens),'like',1i);
    
    % Loop over generators
    for n=1:numel(gens)
        
        % Get the first operator adjoint
        AdA=kron(eye(D,D),gens{n})-kron(transpose(gens{n}),eye(D,D));
        
        % Loop over generators
        for k=1:numel(gens)
            
            % Get the second operator adjoint
            AdB=kron(eye(D,D),gens{k})-kron(transpose(gens{k}),eye(D,D));
            
            % Get the Killing form element
            Kil(n,k)=trace(AdA*AdB);
            
        end
        
    end
    
end
    
end

% Consistency enforcement
function grumble(L)
if (~isnumeric(L))||(~isscalar(L))||(~isreal(L))||...
   (mod(2*L,1)~=0)||(L<0.5)
    error('L must be a positive integer or half-integer.');
end
end

% The FBI is the only organization on Earth complaining
% that computer security is too good.
%
% Matt Blaze

