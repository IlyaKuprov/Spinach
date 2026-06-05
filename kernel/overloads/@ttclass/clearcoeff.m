% Absorbs physical coefficients into tensor train cores without
% changing the value represented by the tensor train. Syntax:
% 
%                     tt=clearcoeff(tt)
%
% Parameters:
%
%    tt - tensor train object
%
% Outputs:
%
%    tt - tensor train object with each coefficient distributed
%         into its cores and the coefficient array set to one
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/clearcoeff.m>

function tt=clearcoeff(tt)

% Get the number of cores and trains
[ncores,ntrains]=size(tt.cores);

% Loop over the trains in the buffer
for n=1:ntrains
    
    % Scale the coefficient
    A=tt.coeff(1,n)^(1/ncores);
    
    % Apply it to cores
    for k=1:ncores
        tt.cores{k,n}=A*tt.cores{k,n};
    end
    
    % Erase the coefficient
    tt.coeff(1,n)=1;
    
end

end

% "Moral outrage is a middle-class luxury."
%
% Anonymous philosophy professor

% #NGRUM

