% Average Hamiltonian theories under Zeeman interaction rotating frame
% transformations. Syntax:
%
%             H=average(spin_system,Hp,H0,Hm,omega,theory)
%
% Parameters:
%
%     Hp  -  the part of the rotating frame Hamiltonian that has positive
%            frequency +omega under the rotating frame transformation
%
%     H0  -  the part of the rotating frame Hamiltonian that has zero
%            frequency under the rotating frame transformation 
%
%     Hm  -  the part of the rotating frame Hamiltonian that has negative
%            frequency -omega under the rotating frame transformation
%
%  omega  -  the frequency of the rotating frame transformation, rad/s
%
%  theory -  the level of the average Hamiltonian theory:
%
%            'ah_first_order'  - first order in Waugh theory
%
%            'ah_second_order' - second order in Waugh theory
%
%            'ah_third_order'  - third order in Waugh theory
%
%            'matrix_log'      - exact algorithm (very expensive,
%                                uses dense matrix algebra)
%
%            'kb_first_order'  - first order in Krylov-Bogolyubov
%                                theory (DNP experiments only)
%
%            'kb_second_order' - second order in Krylov-Bogolyubov
%                                theory (DNP experiments only)
%
%            'kb_second_order' - third order in Krylov-Bogolyubov
%                                theory (DNP experiments only)
%
% Outputs:
%
%      H  -  average Hamiltonian
%
% Note: Krylov-Bogolyubov averging theory as applied to DNP systems is
%       described in detail here:
%
%                http://dx.doi.org/10.1039/C2CP23233B
%                http://dx.doi.org/10.1007/s00723-012-0367-0
%
% i.kuprov@soton.ac.uk
% alexander.karabanov@nottingham.ac.uk
% walter.kockenberger@nottingham.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=average.m>

function H=average(spin_system,Hp,H0,Hm,omega,theory)

% Check consistency
grumble(Hp,H0,Hm,omega,theory);

% Print diagnostics
report(spin_system,['H+ 1-norm: ' num2str(norm(Hp,1))]);
report(spin_system,['H0 1-norm: ' num2str(norm(H0,1))]);
report(spin_system,['H- 1-norm: ' num2str(norm(Hm,1))]);
report(spin_system,['frequency denominator ' num2str(omega/(2*pi)) ' Hz']);

% Run the averaging
switch theory
    
    case 'kb_first_order'
        
        H1=+(1/omega^1)*(Hp*Hm-Hm*Hp); H=H0+H1;
        report(spin_system,['Krylov-Bogolyubov 0 order 1-norm: ' num2str(norm(H0,1))]);
        report(spin_system,['Krylov-Bogolyubov 1 order 1-norm: ' num2str(norm(H1,1))]);
        
    case 'kb_second_order'
        
        H1=+(1/omega^1)*(Hp*Hm-Hm*Hp);
        H2=-(1/omega^2)*(Hp*(Hm*H0-H0*Hm)+Hm*(Hp*H0-H0*Hp));
        H=H0+H1+H2;
        report(spin_system,['Krylov-Bogolyubov 0 order 1-norm: ' num2str(norm(H0,1))]);
        report(spin_system,['Krylov-Bogolyubov 1 order 1-norm: ' num2str(norm(H1,1))]);
        report(spin_system,['Krylov-Bogolyubov 2 order 1-norm: ' num2str(norm(H2,1))]);
         
    case 'kb_third_order'
        
        H1=+(1/omega^1)*(Hp*Hm-Hm*Hp);
        H2=-(1/omega^2)*(Hp*(Hm*H0-H0*Hm)+Hm*(Hp*H0-H0*Hp));
        H3=+(1/omega^3)*(0.5*(Hm*Hm*Hp*Hp-Hp*Hp*Hm*Hm)+...
                             (Hp*Hm+Hm*Hp)*(Hp*Hm-Hm*Hp)+...
                              Hp*(H0*(Hm*H0-H0*Hm)-(Hm*H0-H0*Hm)*H0)-...
                              Hm*(H0*(Hp*H0-H0*Hp)-(Hp*H0-H0*Hp)*H0));
        H=H0+H1+H2+H3;
        report(spin_system,['Krylov-Bogolyubov 0 order 1-norm: ' num2str(norm(H0,1))]);
        report(spin_system,['Krylov-Bogolyubov 1 order 1-norm: ' num2str(norm(H1,1))]);
        report(spin_system,['Krylov-Bogolyubov 2 order 1-norm: ' num2str(norm(H2,1))]);
        report(spin_system,['Krylov-Bogolyubov 3 order 1-norm: ' num2str(norm(H3,1))]);
                                       
    case 'ah_first_order'
        
        H1=-(1/omega)*(H0*Hm-H0*Hp-Hm*H0+Hm*Hp+Hp*H0-Hp*Hm); 
        H=H0+H1;
        report(spin_system,['average Hamiltonian 0 order 1-norm: ' num2str(norm(H0,1))]);
        report(spin_system,['average Hamiltonian 1 order 1-norm: ' num2str(norm(H1,1))]);
        
    case 'ah_second_order'
        
        H1=-(1/omega^1)*(H0*Hm-H0*Hp-Hm*H0+Hm*Hp+Hp*H0-Hp*Hm);
        H2=-(1/omega^2)*(2*H0*H0*Hm + 2*H0*H0*Hp - 4*H0*Hm*H0 - 1*H0*Hm*Hm + 2*H0*Hm*Hp - ... 
                         4*H0*Hp*H0 + 2*H0*Hp*Hm - 1*H0*Hp*Hp + 2*Hm*H0*H0 + 2*Hm*H0*Hm - ... 
                         4*Hm*H0*Hp - 1*Hm*Hm*H0 + 2*Hm*Hm*Hp + 2*Hm*Hp*H0 - 4*Hm*Hp*Hm + ... 
                         2*Hm*Hp*Hp + 2*Hp*H0*H0 - 4*Hp*H0*Hm + 2*Hp*H0*Hp + 2*Hp*Hm*H0 + ... 
                         2*Hp*Hm*Hm - 4*Hp*Hm*Hp - 1*Hp*Hp*H0 + 2*Hp*Hp*Hm)/2;
        H=H0+H1+H2;
        report(spin_system,['average Hamiltonian 0 order 1-norm: ' num2str(norm(H0,1))]);
        report(spin_system,['average Hamiltonian 1 order 1-norm: ' num2str(norm(H1,1))]);
        report(spin_system,['average Hamiltonian 2 order 1-norm: ' num2str(norm(H2,1))]);
        
    case 'ah_third_order'
        
        H1=-(1/omega^1)*(H0*Hm-H0*Hp-Hm*H0+Hm*Hp+Hp*H0-Hp*Hm);
        H2=-(1/omega^2)*(2*H0*H0*Hm + 2*H0*H0*Hp - 4*H0*Hm*H0 - 1*H0*Hm*Hm + 2*H0*Hm*Hp - ... 
                         4*H0*Hp*H0 + 2*H0*Hp*Hm - 1*H0*Hp*Hp + 2*Hm*H0*H0 + 2*Hm*H0*Hm - ... 
                         4*Hm*H0*Hp - 1*Hm*Hm*H0 + 2*Hm*Hm*Hp + 2*Hm*Hp*H0 - 4*Hm*Hp*Hm + ... 
                         2*Hm*Hp*Hp + 2*Hp*H0*H0 - 4*Hp*H0*Hm + 2*Hp*H0*Hp + 2*Hp*Hm*H0 + ... 
                         2*Hp*Hm*Hm - 4*Hp*Hm*Hp - 1*Hp*Hp*H0 + 2*Hp*Hp*Hm)/2;
        H3=+(1/omega^3)*(12*H0*H0*H0*Hm - 12*H0*H0*H0*Hp - 36*H0*H0*Hm*H0 -  9*H0*H0*Hm*Hm + ...
                         12*H0*H0*Hm*Hp + 36*H0*H0*Hp*H0 - 12*H0*H0*Hp*Hm +  9*H0*H0*Hp*Hp + ...
                         36*H0*Hm*H0*H0 + 18*H0*Hm*H0*Hm - 36*H0*Hm*H0*Hp +  2*H0*Hm*Hm*Hm + ...
                         18*H0*Hm*Hm*Hp + 12*H0*Hm*Hp*H0 - 36*H0*Hm*Hp*Hm - 36*H0*Hp*H0*H0 + ...
                         36*H0*Hp*H0*Hm - 18*H0*Hp*H0*Hp - 12*H0*Hp*Hm*H0 + 36*H0*Hp*Hm*Hp - ...
                         18*H0*Hp*Hp*Hm -  2*H0*Hp*Hp*Hp - 12*Hm*H0*H0*H0 + 36*Hm*H0*H0*Hp - ...
                         18*Hm*H0*Hm*H0 -  6*Hm*H0*Hm*Hm - 36*Hm*H0*Hp*H0 + 36*Hm*H0*Hp*Hm - ...
                         18*Hm*H0*Hp*Hp +  9*Hm*Hm*H0*H0 +  6*Hm*Hm*H0*Hm - 18*Hm*Hm*H0*Hp - ...
                          2*Hm*Hm*Hm*H0 +  6*Hm*Hm*Hm*Hp - 18*Hm*Hm*Hp*Hm + 18*Hm*Hm*Hp*Hp + ...
                         12*Hm*Hp*H0*H0 - 36*Hm*Hp*H0*Hm + 36*Hm*Hp*Hm*H0 + 18*Hm*Hp*Hm*Hm - ...
                         36*Hm*Hp*Hm*Hp + 18*Hm*Hp*Hp*H0 +  6*Hm*Hp*Hp*Hp + 12*Hp*H0*H0*H0 - ...
                         36*Hp*H0*H0*Hm + 36*Hp*H0*Hm*H0 + 18*Hp*H0*Hm*Hm - 36*Hp*H0*Hm*Hp + ...
                         18*Hp*H0*Hp*H0 +  6*Hp*H0*Hp*Hp - 12*Hp*Hm*H0*H0 + 36*Hp*Hm*H0*Hp - ...
                         18*Hp*Hm*Hm*H0 -  6*Hp*Hm*Hm*Hm - 36*Hp*Hm*Hp*H0 + 36*Hp*Hm*Hp*Hm - ...
                         18*Hp*Hm*Hp*Hp -  9*Hp*Hp*H0*H0 + 18*Hp*Hp*H0*Hm -  6*Hp*Hp*H0*Hp - ...
                         18*Hp*Hp*Hm*Hm + 18*Hp*Hp*Hm*Hp +  2*Hp*Hp*Hp*H0 -  6*Hp*Hp*Hp*Hm)/12;
        H=H0+H1+H2+H3;
        report(spin_system,['average Hamiltonian 0 order 1-norm: ' num2str(norm(H0,1))]);
        report(spin_system,['average Hamiltonian 1 order 1-norm: ' num2str(norm(H1,1))]);
        report(spin_system,['average Hamiltonian 2 order 1-norm: ' num2str(norm(H2,1))]);
        report(spin_system,['average Hamiltonian 3 order 1-norm: ' num2str(norm(H3,1))]);
        
    case 'matrix_log'
        
        % Discretise the period of the rotating frame, 16 intervals are enough
        nslices=16; time_grid=linspace(0,2*pi/omega,nslices+1); time_step=time_grid(2);
        
        % Compute the product integral using 4th order Lie quadrature
        P=eye(size(H0)); report(spin_system,'computing product integral...');
        parfor n=1:nslices
            report(spin_system,['interval ' num2str(n) ' out of ' num2str(nslices) '...']);
            HL=H0+exp(+1i*omega*(time_grid(n)+0.0*time_step))*Hp+exp(-1i*omega*(time_grid(n)+0.0*time_step))*Hm;
            HM=H0+exp(+1i*omega*(time_grid(n)+0.5*time_step))*Hp+exp(-1i*omega*(time_grid(n)+0.5*time_step))*Hm;
            HR=H0+exp(+1i*omega*(time_grid(n)+1.0*time_step))*Hp+exp(-1i*omega*(time_grid(n)+1.0*time_step))*Hm;
            P=propagator(spin_system,isergen(HL,HM,HR,time_step),time_step)*P;
        end
        
        % Compute matrix logarithm (slow, needs work)
        report(spin_system,'computing matrix logarithm...');
        H=1i*(omega/(2*pi))*logm(P);
        
    otherwise
        
        % Complain and bomb out
        error('unknown theory specification.');
        
end

% Clean up the result and force the matrix into sparse format
H=clean_up(spin_system,H,spin_system.tols.liouv_zero); H=sparse(H);

% Report to the user
report(spin_system,['average Hamiltonian dimension ' num2str(size(H,1))...
                    ', nnz ' num2str(nnz(H)) ', density ' num2str(100*nnz(H)/numel(H))...
                    '%, 1-norm ' num2str(norm(H,1)) ', sparsity ' num2str(issparse(H))]); 

end

% Consistency enforcement
function grumble(Hp,H0,Hm,omega,theory)
if (~isnumeric(Hp))||(~isnumeric(H0))||(~isnumeric(Hm))||...
   (~ismatrix(Hp))||(~ismatrix(H0))||(~ismatrix(Hm))
    error('the three Hamiltonian components must be matrices.');
end
if any(size(Hp)~=size(H0))||any(size(H0)~=size(Hm))
    error('dimensions of the three Hamiltonian components must be the same.');
end
if (~isnumeric(omega))||(~isreal(omega))||(~isscalar(omega))
    error('omega parameter must be a real number.');
end
if ~ischar(theory)
    error('theory parameter must be a character string.');
end
end

% I save about twenty drafts - that's ten meg of disc space - and the last
% one contains all the final alterations. Once it has been printed out and
% received by the publishers, there's a cry here of "Tough shit, literary
% researchers of the future, try getting a proper job!" and the rest are
% wiped.
%
% Terry Pratchett

