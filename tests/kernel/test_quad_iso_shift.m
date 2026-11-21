% Test: kernel/quad_iso_shift | Quadrupolar powder shifts
% Checks second-order powder shifts with exact oriented energies.
% Syntax:
%
%                 result=test_quad_iso_shift()
%
% Outputs:
%
%     result - regression test result with explanatory messages
%
% talos@spindynamics.org

function result=test_quad_iso_shift()

% Set up the signed Zeeman coefficient and a non-axial tensor
result=new_test_result('kernel/quad_iso_shift','Quadrupolar powder shifts',...
                       'central and satellite corrections approach exact high-field energies.');
cq_hz=2e6;
eta=0.64;
fields=[80e6 320e6];
m=[1.5 0.5 -0.5 -1.5];
lz=diag(m);
lp=diag([sqrt(3) 2 sqrt(3)],1);
operators={(lp+lp')/2,(lp-lp')/(2i),lz};
pairs=[2 3; 1 2];

% Allow exact integration of fourth-degree polynomials in cos(beta)
nodes=[-0.861136311594053 -0.339981043584856 0.339981043584856 0.861136311594053];
weights=[0.347854845137454 0.652145154862546 0.652145154862546 0.347854845137454]/2;
errors=zeros(2,2);
for n=1:numel(fields)

    % Form an oriented full Hamiltonian, then track levels by Zeeman ordering
    nu_hz=fields(n);
    shifts=zeros(1,4);
    for alpha=(0:7)*pi/4
        for k=1:numel(nodes)
            for gamma=(0:7)*pi/4
                Q=eeqq2nqi(cq_hz,eta,3/2,[alpha acos(nodes(k)) gamma]);
                H=zeros(4);
                for a=1:3
                    for b=1:3
                        H=H+Q(a,b)*operators{a}*operators{b};
                    end
                end
                energies=sort(eig(nu_hz*lz+(H+H')/2),'descend');
                shifts=shifts+weights(k)*(energies.'-nu_hz*m)/64;
            end
        end
    end

    % Compare both transitions to the exact powder energy differences
    for k=1:2
        levels=pairs(k,:);
        observed=quad_iso_shift(3/2,cq_hz,eta,nu_hz,m(levels(1)),m(levels(2)));
        expected=diff(shifts(levels));
        errors(n,k)=abs((observed-expected)/observed);
        result=test_true(result,['signed transition ' int2str(k) ' field ' int2str(n)],...
                         sign(observed)==sign(expected)&&errors(n,k)<0.005,...
                         'full-Hamiltonian powder shifts agree in Hz and sign with perturbation theory');
    end

end

% Check asymptotic convergence, signed frequency reversal, and zero coupling
result=test_true(result,'high-field convergence',all(errors(2,:)<errors(1,:)),...
                 'the relative error decreases at higher field for both transitions');
plus=quad_iso_shift(3/2,cq_hz,eta,fields(1),0.5,-0.5);
minus=quad_iso_shift(3/2,cq_hz,eta,-fields(1),0.5,-0.5);
result=test_close(result,'signed field reversal',minus,-plus,0,0,...
                  'the second-order correction reverses with the signed Zeeman coefficient');
result=test_close(result,'zero coupling',quad_iso_shift(3/2,0,eta,fields(1),0.5,-0.5),...
                  0,0,0,'without quadrupolar interactions the shift vanishes');

end


