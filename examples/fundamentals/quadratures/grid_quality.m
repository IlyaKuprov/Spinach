% Performance analysis for the spherical and SO(3) integration
% grids supplied with Spinach kernel.
%
% i.kuprov@soton.ac.uk

function grid_quality()

%% Two-angle REPULSION grids

% Create a figure
figure(); hold on; kgrid; box on;
set(gca,'YScale','log'); legend_txt={};

% Loop over two-angle REPULSION grids
for n=[100 200 400 800 1600 3200 6400 12800]
    
    % Load the grid
    load(['../../../kernel/grids/rep_2ang_' num2str(n) 'pts_sph.mat'],...
          'alphas','betas','gammas','weights');
      
    % Evaluate the grid 
    grid_profile=grid_test(alphas,betas,gammas,weights,0:80,'Y_lm');
    
    % Plot the evaluation
    plot(grid_profile(2:end)); axis tight; drawnow;
    legend_txt{end+1}=['REPULSION ' num2str(n) ' pts']; %#ok<AGROW>
    
end

% Residual cosmetics
xlabel('spherical rank');
ylabel('integration error');
legend(legend_txt,'Location','southeast');

%% Three-angle REPULSION grids

% Create a figure
figure(); hold on; kgrid; box on;
set(gca,'YScale','log'); legend_txt={};

% Loop over three-angle REPULSION grids
for n=[100 200 400 800 1600 3200 6400 12800]
    
    % Load the grid
    load(['../../../kernel/grids/rep_3ang_' num2str(n) 'pts.mat'],...
          'alphas','betas','gammas','weights');
      
    % Evaluate the grid 
    grid_profile=grid_test(alphas,betas,gammas,weights,0:30,'D_lmn');
    
    % Plot the evaluation
    plot(grid_profile(2:end)); axis tight; drawnow;
    legend_txt{end+1}=['REPULSION ' num2str(n) ' pts']; %#ok<AGROW>
    
end

% Residual cosmetics
xlabel('spherical rank');
ylabel('integration error');
legend(legend_txt,'Location','southeast');

%% Two-angle Lebedev grids
 
% Create a figure
figure(); hold on; kgrid; box on;
set(gca,'YScale','log'); legend_txt={};

% Loop over two-angle Lebedev grids
for n=[5 17 29 41 53]
    
    % Load the grid
    load(['../../../kernel/grids/leb_2ang_rank_' num2str(n) '.mat'],...
          'alphas','betas','gammas','weights');
      
    % Set even ranks
    ranks=2:2:100;
      
    % Evaluate the grid 
    grid_profile=grid_test(alphas,betas,gammas,weights,ranks,'Y_lm');
    
    % Plot the evaluation
    plot(ranks(2:end),grid_profile(2:end)); axis tight; drawnow;
    legend_txt{end+1}=['Lebedev rank ' num2str(n)]; %#ok<AGROW>
    
end

% Residual cosmetics
xlabel('spherical rank');
ylabel('integration error');
legend(legend_txt,'Location','southeast');

end

