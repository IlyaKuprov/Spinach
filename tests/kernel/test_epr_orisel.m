% Focused orientation-selection regression tests. Syntax:
%
%                   result=test_epr_orisel()
%
% Output: standard Spinach regression result; failures throw errors.
%
% Checks Gaussian FWHM, field scaling, selectors, removed options,
% quadrature pruning, and visible true-colour plotting on a small grid.
%
% Talos

function result=test_epr_orisel()

% Check the fixed no-input test interface
grumble();
result=new_test_result('kernel/epr_orisel','EPR orientation selection',...
                       'Validate linewidth, field scaling, quadrature, and plotting.');

% Build an isotropic electron at a deliberately different reference field
sys.isotopes={'E'}; sys.magnet=0.7; sys.output='hush';
sys.parallel={'local',1}; sys.disable={'hygiene'};
inter.zeeman.eigs={[2 2 2]}; inter.zeeman.euler={[0 0 0]};
inter.temperature=20;
bas.formalism='zeeman-hilb'; bas.approximation={'none'};
s=create(sys,inter); s=basis(s,bas);
p.grid='leb_2ang_rank_5'; p.selection_spin='E';
p.exclude_from_EPR=[]; p.magnetic=1.21; p.mw_frq=34e9;
p.bandwidth_EPR=100e6; p.rspt_order=inf; p.intensity_cutoff=0;
[g,t]=epr_orisel_transitions(s,p);

% Check an absolute S=1/2 intensity using a pure ground state
h=hamiltonian(assume(s,'labframe'));
[v,e]=eig(full(h),'vector'); [~,ground]=min(e);
pure=p; pure.rho0=v(:,ground)*v(:,ground)';
[~,pure_trans]=epr_orisel_transitions(s,pure);
strengths=cellfun(@(x)x(1),pure_trans.intensities);
fprintf('ABSOLUTE INTENSITY actual %.12g, expected 0.25\n',strengths(1));
assert(max(abs(strengths-0.25))<1e-12);

% Check unequal orientation-dependent strengths by direct diagonalisation
inter_a=inter; inter_a.zeeman.eigs={[2.15 1.91 2.03]};
inter_a.zeeman.euler={[0.2 0.4 0.7]};
sa=create(sys,inter_a); sa=basis(sa,bas);
[ga,ta]=epr_orisel_transitions(sa,pure);
[iz,qz]=hamiltonian(assume(sa,'labframe','zeeman'));
sx=operator(sa,'Lx',1); reference=zeros(size(ga.weights));
for n=1:numel(ga.weights)
    hz=iz+orientation(qz,[ga.alphas(n) ga.betas(n) ga.gammas(n)]);
    [va,ea]=eig(full(hz*p.magnetic/sys.magnet),'vector');
    [~,order]=sort(ea); va=va(:,order);
    populations=real(diag(va'*pure.rho0*va));
    reference(n)=abs(va(:,1)'*sx*va(:,2))^2*(populations(1)-populations(2));
    assert(abs(ta.intensities{n}(1)-reference(n))<1e-12);
end
assert(max(reference)>min(reference));
fprintf('ANISOTROPIC INTENSITY range %.12g to %.12g\n',...
        min(reference),max(reference));

% Verify the Zeeman frequency and its field derivative independently
expected=abs(spin('E'))*(2/s.tols.freeg)*p.magnetic/(2*pi);
assert(max(abs(cellfun(@(x)x(1),t.frequencies)-expected))<1e-4);
assert(max(abs(cellfun(@(x)x(1),t.dfdb_hz_per_t)-expected/p.magnetic))<1e-4);
fprintf('PASS reference-field scaling: %.9f GHz at %.4f T\n',expected/1e9,p.magnetic);

% Verify area normalisation and half-height points through production
p.mw_frq=t.frequencies{1}(1);
[~,peak]=epr_orisel_transitions(s,p);
peak_height=sqrt(4*log(2)/pi)/p.bandwidth_EPR;
expected=peak_height*cellfun(@(x)x(1),peak.intensities);
assert(max(abs(peak.selection_weights./expected-1))<1e-12);
for sign=[-1 1]
    half=p; half.mw_frq=p.mw_frq+sign*p.bandwidth_EPR/2;
    [~,t]=epr_orisel_transitions(s,half);
    ratio=t.selection_weights./peak.selection_weights;
    assert(max(abs(ratio-.5))<1e-10);
    fprintf('PASS Gaussian half-height at %+g MHz: %.12f\n',sign*50,ratio(1));
end

% Verify that obsolete fields no longer affect physics or row selection
old=p; old.selection_exclude=1; old.only_ENDOR_list=1; old.epr_transition=99;
[~,t]=epr_orisel_transitions(s,old);
assert(isequal(t.selection_weights,peak.selection_weights));
old.grid_cutoff=0;
[~,sel]=epr_orisel_grid(s,old);
assert(isequal(sel.transition_indices,1));
assert(isequal(sel.selection_weights,sel.weights));
assert(isequal(sel.frequencies,peak.frequencies)&&isequal(sel.grid,g));
export_file=[tempname '.mat'];
cleanup=onCleanup(@()delete(export_file));
old.export_grid=export_file; sg=epr_orisel_grid(s,old);
assert(isequal(load(export_file),sg));
clear cleanup
fprintf('PASS removed options have no effect\n');

% Check weighted pruning at an exact cutoff boundary
q.selection_weights=zeros(1,numel(g.weights));
q.selection_weights(1:3)=[.5 .3 .2]./g.weights(1:3)';
q.grid_cutoff=.5;
sg=epr_orisel_grid(g,q);
assert(isscalar(sg.weights)&&sg.weights==1);
q.grid_cutoff=0; sg=epr_orisel_grid(g,q);
assert(numel(sg.weights)==3&&abs(sum(sg.weights)-1)<eps);
assert(max(abs(sg.weights-[.5;.3;.2]))<1e-14);
q.selection_weights(:)=0;
assert_fails(@()epr_orisel_grid(g,q));
q.selection_weights(1)=-1;
assert_fails(@()epr_orisel_grid(g,q));
q.selection_weights(1)=1; q.grid_cutoff=1;
assert_fails(@()epr_orisel_grid(g,q));

% Reject nonfinite, complex, empty and mismatched supplied responses
q.grid_cutoff=0;
for bad=[NaN inf 1i]
    q.selection_weights(1)=bad;
    assert_fails(@()epr_orisel_grid(g,q));
end
q.selection_weights=[]; assert_fails(@()epr_orisel_grid(g,q));
q.selection_weights=1; assert_fails(@()epr_orisel_grid(g,q));
q.selection_weights=ones(1,numel(g.weights));
for bad=[-eps NaN 1i]
    q.grid_cutoff=bad;
    assert_fails(@()epr_orisel_grid(g,q));
end
q.grid_cutoff=1-eps; sg=epr_orisel_grid(g,q);
assert(isscalar(sg.weights)&&sg.weights==1);
fprintf('PASS quadrature weighting, cutoff equality, zero and invalid inputs\n');

% Check retained exclusion selectors on an electron-nitrogen system
sys.isotopes={'E','14N'}; sys.labels={'radical','nitrogen'};
inter.zeeman.eigs={[2 2 2],[0 0 0]};
inter.zeeman.euler={[0 0 0],[0 0 0]};
inter.coupling.matrix=cell(2);
inter.coupling.matrix{1,2}=diag([14 14 95])*1e6;
sn=create(sys,inter); sn=basis(sn,bas);
selectors={2,[false true],'14N',{'nitrogen'}};
for k=1:numel(selectors)
    ex=p; ex.exclude_from_EPR=selectors{k};
    [~,t]=epr_orisel_transitions(sn,ex);
    assert(max(abs(t.selection_weights-peak.selection_weights))<1e-14);
end
ex=p; ex.exclude_from_EPR=1;
assert_fails(@()epr_orisel_transitions(sn,ex));
ex=p; ex.rspt_order=1.5;
assert_fails(@()epr_orisel_transitions(sn,ex));
ex=p; ex.selection_spin='14N';
assert_fails(@()epr_orisel_transitions(sn,ex));
ex=p; ex.selection_spin=2;
assert_fails(@()epr_orisel_transitions(sn,ex));

% Check nonfinite scalar controls and the preserved selector entry point
for field={'mw_frq','magnetic','bandwidth_EPR'}
    for bad=[0 NaN inf]
        ex=p; ex.(field{1})=bad;
        assert_fails(@()epr_orisel_transitions(sn,ex));
    end
end
selector.spin_selector={'E',2};
assert(isequal(epr_orisel_transitions(sn,selector),[1 2]));
selector.spin_selector=[1;2];
assert(isequal(epr_orisel_transitions(sn,selector),[1 2]));
ex=rmfield(p,'selection_spin'); [~,t]=epr_orisel_transitions(s,ex);
assert(isequal(t.selection_weights,peak.selection_weights));
fprintf('PASS exclusion selectors and unsupported inputs\n');

% Check both adjacent manifolds of a spin-one electron
sys.isotopes={'E3'}; sys=rmfield(sys,'labels');
inter=rmfield(inter,'coupling');
inter.zeeman.eigs={[2 2 2]}; inter.zeeman.euler={[0 0 0]};
se=create(sys,inter); se=basis(se,bas);
pe=p; pe.selection_spin='E3'; pe.selection_transition=1;
[~,first]=epr_orisel_grid(se,pe);
pe.selection_transition=2; [~,second]=epr_orisel_grid(se,pe);
assert(size(first.selection_weights,1)==2);
assert(isequal(first.selected_weights,first.selection_weights(1,:)));
assert(isequal(second.selected_weights,second.selection_weights(2,:)));
assert(all(first.selected_weights>0)&&all(second.selected_weights>0));
pe.selection_transition=3;
assert_fails(@()epr_orisel_grid(se,pe));
fprintf('PASS spin-one manifold row selection\n');

% Check true colours, finite zero-weight plotting, ticks, and hold state
f=figure('Visible','off'); cleanup=onCleanup(@()close(f));
w=linspace(0,1,numel(g.weights)); plot_orisel(g,w,1);
a=gca; patches=findobj(a,'Type','patch'); cmap=parula(256);
assert(numel(patches)==numel(w));
for k=1:numel(patches)
    rgb=patches(k).FaceColor;
    assert(size(rgb,2)==3&&any(all(abs(cmap-rgb)<1e-14,2)));
end
assert(~isempty(a.XTick)&&~isempty(a.YTick)&&~isempty(a.ZTick));
assert(~isequal(a.XColor,'none')&&~isequal(a.YColor,'none')&&~isequal(a.ZColor,'none'));
assert(~ishold); cla; hold on; plot_orisel(g,zeros(size(w)),1); assert(ishold);
clear cleanup
fprintf('PASS RGB plotting, zero weights, visible axes, and hold state\n');
result=test_true(result,'Orientation selection',true,...
                 'All numerical, input, and graphics assertions passed.');

end

% Require a clear rejection for unsupported input
function assert_fails(call)

% Record the production error without substituting a successful result
failed=false;
try
    call();
catch
    failed=true;
end
assert(failed,'unsupported input was not rejected.');

end

% Check the fixed no-input interface
function grumble()
end


