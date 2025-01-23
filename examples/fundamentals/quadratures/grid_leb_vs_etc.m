% Heuristic vs Lebedev spherical quadrature bake-off, illus-
% trating the fact that, well... heuristic grids suck.
%
% ilya.kuprov@weizmann.ac.il

function grid_leb_vs_etc()

% Evaluate Lebedev grid
load('../../../kernel/grids/leb_2ang_rank_29.mat',...
     'alphas','betas','gammas','weights');
perf_leb=grid_test(alphas,betas,gammas,weights,4:2:60,'Y_lm');
 
% Evaluate Repulsion grid with Voronoi weights
[alphas,betas,gammas]=repulsion(numel(alphas),3,10000);
[~,~,~,weights]=voronoisphere([sin(betas').*cos(gammas');
                               sin(betas').*sin(gammas');
                               cos(betas')]);
weights=weights/(4*pi);
perf_rep=grid_test(alphas,betas,gammas,weights,4:2:60,'Y_lm');

% Evaluate ZCWn grid with Voronoi weights
[alphas,betas,gammas,weights]=grid_fibon('zcwn',302);
perf_zcw=grid_test(alphas,betas,gammas,weights,4:2:60,'Y_lm');

% Evaluate Igloo grid with Voronoi weights
[alphas,betas,gammas,weights]=grid_igloo(17);
perf_igl=grid_test(alphas,betas,gammas,weights,4:2:60,'Y_lm');

% Evaluate Stoll grid with Voronoi weights
[alphas,betas,gammas,weights]=grid_trian('stoll',9);
perf_esp=grid_test(alphas,betas,gammas,weights,4:2:60,'Y_lm');

% Evaluate ASG grid with Voronoi weights
[alphas,betas,gammas,weights]=grid_trian('asg',9);
perf_asg=grid_test(alphas,betas,gammas,weights,4:2:60,'Y_lm');

% Evaluate SOPHE grid with Voronoi weights
[alphas,betas,gammas,weights]=grid_trian('sophe',9);
perf_sop=grid_test(alphas,betas,gammas,weights,4:2:60,'Y_lm');

% Plot the profiles
figure(); hold on;
plot(4:2:60,perf_rep);
plot(4:2:60,perf_esp);
plot(4:2:60,perf_leb,'k-');
plot(4:2:60,perf_sop);
plot(4:2:60,perf_zcw);
plot(4:2:60,perf_igl);
plot(4:2:60,perf_asg);
kgrid; ylim([1e-16 1]);
set(gca,'YScale','log');
xlim([0 64]); box on;

% Residual cosmetics
kxlabel('spherical rank');
kylabel('integration error');
legend({'Repulsion, 302 pts',...
        'EasySpin, 326 pts',...
        'Lebedev, 302 pts',...
        'SOPHE, 326 pts',...
        'ZCWn, 302 pts',...
        'Igloo, 328 pts',...
        'ASG, 326 pts'},...
        'Location','southeast');
scale_figure([1.0 0.75]);

end

