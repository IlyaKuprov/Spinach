% Wigner matrix element correlation function under isotropic, axial,
% and rhombic rotational diffusion. Syntax:
%
%       [weights,rates,states]=corrfun(spin_system,n,k,m,p,q)
%
% Parameters:
%
%      spin_system - the output of [[create.m]] to which ro-
%                    tational correlation time should have 
%                    been supplied. For a single correlation
%                    time, the isotropic rotational diffusi-
%                    on model is used; a vector with two
%                    correlation times is assumed to be cor-
%                    relation times for rotation around and
%                    perpendicularly to the main axis res-
%                    pectively); a vector with three corre-
%                    lation times is assumed to be the cor-
%                    relation times for the rotation around
%                    the XX, YY and ZZ direction respecti-
%                    vely of the rotational diffusion tensor.
%
%      n,k,m,p,q   - the five indices found in the ensemble-
%                    averaged Wigner function product:
%
%                          <D{n}{k,m}(0)*D{n}{p,q}(t)'>
% 
% Outputs:
%
%      weights     - a cell array (one element for each che-
%                    mical species) of vectors listing the 
%                    weights of the exponential components
%                    of the decays
%
%      rates       - a cell array (one element for each che-
%                    mical species) of vectors listing the 
%                    decay rates (negative numbers) of the
%                    exponential components of the decays
%
%      states      - a cell array (one element for each che-
%                    mical species) of logical vectors indi-
%                    cating which states in the basis set
%                    belong to which chemical species
%
% Note: Wigner function indices are sorted in descending order, that is,
%       k=[1 2 3 4 5] in the input represents [2 1 0 -1 -2] for n=2.
%
% Note: second rank rotational correlation times (as per Spinach input) 
%       will be updated automatically if other ranks are specified.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=corrfun.m>

function [weights,rates,states]=corrfun(spin_system,n,k,m,p,q)

% Check consistency
grumble(spin_system,n,k,m,p,q);

% Preallocate outputs
weights=cell(1,numel(spin_system.chem.parts));
rates=cell(1,numel(spin_system.chem.parts));
states=cell(1,numel(spin_system.chem.parts));

% Index basis states for different chemical species
for s=1:numel(spin_system.chem.parts)
    states{s}=(sum(spin_system.bas.basis(:,spin_system.chem.parts{s}),2)>0);
end

% Loop over chemical species
for s=1:numel(spin_system.chem.parts)
    
    % Select rotational diffusion model
    switch numel(spin_system.rlx.tau_c{s})
        
        case 1
            
            % Assume second rank correlation time is supplied
            D=1/(6*spin_system.rlx.tau_c{s});
            
            % Use isotropic rotational diffusion model
            weights{s}=(1/(2*n+1))*krondelta(k,p)*krondelta(m,q);
            rates{s}=-n*(n+1)*D;
            
        case 2
            
            % Assume second rank correlation time is supplied
            D_ax=1/(6*spin_system.rlx.tau_c{s}(1));
            D_eq=1/(6*spin_system.rlx.tau_c{s}(2));
            
            % Use axial rotational diffusion model
            weights{s}=(1/(2*n+1))*krondelta(k,p)*krondelta(m,q);
            rates{s}=-(n*(n+1)*D_eq+((n-m+1)^2)*(D_ax-D_eq));
            
        case 3
            
            % Use anisotropic rotational diffusion model, only L=2 permitted
            Dxx=1/(6*spin_system.rlx.tau_c{s}(1));
            Dyy=1/(6*spin_system.rlx.tau_c{s}(2));
            Dzz=1/(6*spin_system.rlx.tau_c{s}(3));
            
            % Refuse to process degenerate cases and ranks other than 2
            if (abs(Dxx-Dyy)<1e-6*mean([Dxx Dyy Dzz]))||...
               (abs(Dyy-Dzz)<1e-6*mean([Dxx Dyy Dzz]))||...
               (abs(Dzz-Dxx)<1e-6*mean([Dxx Dyy Dzz]))
                error('the three rotational correlation times must be different.');
            end
            if n~=2
                error('rhombic rotational diffusion only implemented for n=2.');
            end
            
            % Compute decay rates
            delta=sqrt(Dxx^2+Dyy^2+Dzz^2-Dxx*Dyy-Dxx*Dzz-Dyy*Dzz);
            rates{s}(1)=-(4*Dxx+Dyy+Dzz);
            rates{s}(2)=-(Dxx+4*Dyy+Dzz);
            rates{s}(3)=-(Dxx+Dyy+4*Dzz);
            rates{s}(4)=-(2*Dxx+2*Dyy+2*Dzz-2*delta);
            rates{s}(5)=-(2*Dxx+2*Dyy+2*Dzz+2*delta);
            
            % Compute coefficients
            lambda_p=sqrt(2/3)*(Dxx+Dyy-2*Dzz+2*delta)/(Dxx-Dyy);
            lambda_m=sqrt(2/3)*(Dxx+Dyy-2*Dzz-2*delta)/(Dxx-Dyy);
            h(1,2)=1/sqrt(2); h(1,4)=1/sqrt(2); h(2,2)=-1/sqrt(2);
            h(2,4)=1/sqrt(2); h(3,5)=1/sqrt(2); h(3,1)=-1/sqrt(2);
            h(4,1)=1/sqrt(2+lambda_m^2); h(4,3)=lambda_m/sqrt(2+lambda_m^2); h(4,5)=1/sqrt(2+lambda_m^2);
            h(5,1)=1/sqrt(2+lambda_p^2); h(5,3)=lambda_p/sqrt(2+lambda_p^2); h(5,5)=1/sqrt(2+lambda_p^2);
            
            % Compute weights
            for j=1:5, weights{s}(j)=(1/5)*krondelta(k,p)*h(j,m)*h(j,q); end
            
    end
    
end

end

% Consistency enforcement
function grumble(spin_system,n,k,m,p,q)
if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
    error('Redfield relaxation theory is only available for sphten-liouv formalism.');
end
if (~isnumeric(n))||(~isnumeric(k))||(~isnumeric(m))||...
   (~isnumeric(p))||(~isnumeric(q))
    error('all indices must be numeric.');
end
if ~isreal([n k m p q])
    error('all indices must be real.');
end
if any(mod([n k m p q],1)~=0)
    error('all indices must be integers.');
end
if n<0
    error('Wigner function rank must be greater than zero.');
end
if any([k m p q]>(2*n+1))||any([k m p q]<1)
    error('k,m,p,q must be from [1,2*n+1] interval.');
end
end

% "But they are useless. They can only give you answers."
%
% Pablo Picasso, about computers.

