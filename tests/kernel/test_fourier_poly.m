% Tests spectral polyadics and shaped RF pulse integration. Syntax:
%
%                       result=test_fourier_poly()
%
% Outputs:
%
%    result - differentiation and pulse propagation regression checks
%
% ilya.kuprov@weizmann.ac.il

function result=test_fourier_poly()

% Describe the independently materialised references
result=new_test_result('kernel/fourier_poly','Spectral polyadic derivatives',...
                       'Fourier actions and RF pulse observables must retain explicit-path results.');

% Compare arbitrary complex blocks on odd and even periodic grids
spin_system.sys.enable={'polyadic'}; explicit_system.sys.enable={};
for npoints=[1 2 3 8 9 32]
    for order=1:4
        [~,D]=fourdif(explicit_system,npoints,order);
        D=(2*pi/3)^order*D;
        [~,implicit]=fourdif(spin_system,npoints,order);
        implicit=(2*pi/3)^order*implicit;
        rhs=reshape(sin(1:(3*npoints))+1i*cos(1:(3*npoints)),npoints,3);
        label=[num2str(npoints) '/' num2str(order)];
        result=test_close(result,['action ' label],implicit*rhs,D*rhs,...
                          1e-9,1e-11,'complex block action matches fourdif');
        result=test_close(result,['adjoint ' label],implicit'*rhs,D'*rhs,...
                          1e-9,1e-11,'transform adjoints have the correct normalisation');
        result=test_close(result,['transpose ' label],implicit.'*rhs,D.'*rhs,...
                          1e-9,1e-11,'non-conjugating transpose preserves complex states');
    end
end

% Build a physical one-spin Liouville system for RF pulse comparisons
sys.magnet=1; sys.isotopes={'1H'}; inter.zeeman.scalar={0};
bas.formalism='zeeman-liouv'; bas.approximation={'none'};
spin_system=test_spin_system(sys,inter,bas);
Lx=operator(spin_system,'Lx',1); Ly=operator(spin_system,'Ly',1);
L0=2*pi*37*operator(spin_system,'Lz',1);
rho=[state(spin_system,'Lz',1) state(spin_system,'L+',1)];
methods={'expv','evolution','expm'};

% Exercise signed frequency slices, phase, rank, and stacked states
for rank=[1 3]
    for n=1:numel(methods)
        spin_system.sys.enable={};
        [reference,ref_traj]=shaped_pulse_af(spin_system,L0,Lx,Ly,rho,...
                                            [130 -70 0],[500 230 420],...
                                            [2e-4 3e-4 1e-4],pi/7,rank,methods{n});
        spin_system.sys.enable={'polyadic'};
        [observed,obs_traj]=shaped_pulse_af(spin_system,L0,Lx,Ly,rho,...
                                           [130 -70 0],[500 230 420],...
                                           [2e-4 3e-4 1e-4],pi/7,rank,methods{n});
        label=[methods{n} '/' num2str(rank)];
        result=test_close(result,['pulse ' label],observed,reference,...
                          1e-9,1e-9,'the complete RF pulse retains its final states');
        result=test_close(result,['trajectory ' label],cell2mat(obs_traj),cell2mat(ref_traj),...
                          1e-9,1e-9,'every folded trajectory point agrees');
    end
end

% Verify that effective propagators remain explicitly available
spin_system.sys.enable={};
[~,~,reference]=shaped_pulse_af(spin_system,L0,Lx,Ly,rho,...
                              [130 -70],[500 230],[2e-4 3e-4],pi/7,2,'expm');
spin_system.sys.enable={'polyadic'};
[observed,~,P]=shaped_pulse_af(spin_system,L0,Lx,Ly,rho,...
                             [130 -70],[500 230],[2e-4 3e-4],pi/7,2,'expm');
result=test_close(result,'effective propagator',P,reference,1e-12,1e-12,...
                  'expm retains its explicit propagator even when polyadics are enabled');
result=test_close(result,'effective action',observed,P*rho,1e-9,1e-9,...
                  'the effective propagator reproduces the returned state stack');

end


