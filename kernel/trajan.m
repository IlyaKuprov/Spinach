% Trajectory analysis function. Plots the time dependence of the densi-
% ty matrix norm, partitioned into user-specified property classes. See
%        
%              http://dx.doi.org/10.1016/j.jmr.2013.02.012
%
% for further information. Syntax:
%
%                trajan(spin_system,traj,property,time_axis)
%
% Arguments:
%
%     traj       - a stack of state vectors of any length. The
%                  number of rows in the trajectory array must
%                  match the number of states in the basis.
%
%     property   - if set to 'correlation_order', returns the
%                  time dependence of the total populations of
%                  one-spin, two-spin, three-spin, etc. corre-
%                  lations in the system.
%
%                  if set to 'coherence_order', returns the ti-
%                  me dependence of different orders of coheren-
%                  ce in the system, where a coherence order is
%                  defined as the sum of projection quantum num-
%                  bers in the spherical tensor representation
%                  of each state.
%
%                  if set to 'total_each_spin', returns the ti-
%                  me dependence of total state space populati-
%                  on that involves each individual spin in the
%                  system in any way (all local populations and
%                  coherences of the spin as well as all of its 
%                  correlations to all third party spins).
%
%                  if set to 'local_each_spin', returns the ti-
%                  me dependence of the population of the sub-
%                  space of states that are local to each indi-
%                  vidual spin and do not involve any correla-
%                  tions to other spins in the system.
%
%                  if set to 'level_populations', returns the
%                  populations of the Zeeman energy levels.
%
%      time_axis - (optional) user specified time axis, a row 
%                  vector of time positions of each state vec-
%                  tor inthe trajectory array.
%
% The trajectory would usually come out of the evolution.m run from a
% given starting point under a given Liouvillian.
%
% Output:
%
%     this function writes into the current figure
%
% Note: this function is only applicable to the trajectories recorded
%       in sphten-liouv formalism.
%
% Note: unit state population is ignored.
%
% gareth.charnock@oerc.ox.ac.uk
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=trajan.m>

function trajan(spin_system,traj,property,time_axis)

% Check the input
grumble(spin_system,traj);

% Set the defaults
if ~exist('property','var'), property='correlation_order'; end

% Project out unit state
if ~strcmp(property,'level_populations')
    unit=unit_state(spin_system);
    traj=traj-(unit*unit')*traj;
end

% Determine how to proceed
switch property
    
    case 'correlation_order'
        
        % Determine the correlation order of each state
        correlation_orders=sum(logical(spin_system.bas.basis),2);
        
        % Find out which correlation orders are present
        unique_correlation_orders=unique(correlation_orders);
        
        % Eliminate zero-spin order
        unique_correlation_orders=setdiff(unique_correlation_orders,0);
        
        % Preallocate the norm trajectory array
        result=zeros(numel(unique_correlation_orders),size(traj,2));
        
        % Loop over the unique correlation orders that are present
        for n=1:numel(unique_correlation_orders)
            
            % Find the subspace of the given correlation order
            subspace_mask=(correlation_orders==unique_correlation_orders(n));
            
            % Get the part of the trajectory belonging to the subspace
            subspace_trajectory=traj(subspace_mask,:);
            
            % Get the norm of the trajectory
            result(n,:)=sqrt(sum(subspace_trajectory.*conj(subspace_trajectory),1));
            
        end
        
        % Create a legend
        legend_text=cell(numel(unique_correlation_orders),1);
        for n=1:numel(unique_correlation_orders)
            legend_text{n}=[num2str(unique_correlation_orders(n)) '-spin'];
        end
        
        % Create labels
        label_text='correlation order amplitude';
        title_text='correlation orders';

        % Compute Y axis extents
        max_val=max(result,[],'all');
        y_axis_extents=[-0.05*max_val 1.05*max_val];
        
    case 'coherence_order'
        
        % Determine projection quantum numbers of the basis
        [~,M]=lin2lm(spin_system.bas.basis);
        
        % Determine the coherence order of each state
        coherence_orders=sum(M,2);
        
        % Find out which coherence orders are present
        unique_coherence_orders=unique(coherence_orders);
        
        % Preallocate the norm trajectory array
        result=zeros(numel(unique_coherence_orders),size(traj,2));
        
        % Loop over the unique coherence orders that are present
        for n=1:numel(unique_coherence_orders)
            
            % Find the subspace of the given coherence order
            subspace_mask=(coherence_orders==unique_coherence_orders(n));
            
            % Get the part of the trajectory belonging to the subspace
            subspace_trajectory=traj(subspace_mask,:);
            
            % Get the norm of the trajectory
            result(n,:)=sqrt(sum(subspace_trajectory.*conj(subspace_trajectory),1));
            
        end
        
        % Create a legend
        legend_text=cell(numel(unique_coherence_orders),1);
        for n=1:numel(unique_coherence_orders)
            legend_text{n}=num2str(unique_coherence_orders(n));
        end
        
        % Create labels
        label_text='coherence order amplitude';
        title_text='coherence orders';

        % Compute Y axis extents
        max_val=max(result,[],'all');
        y_axis_extents=[-0.05*max_val 1.05*max_val];
        
    case 'total_each_spin'
        
        % Preallocate the norm trajectory array
        result=zeros(spin_system.comp.nspins,size(traj,2));
        
        % Loop over spins in the system
        for n=1:spin_system.comp.nspins
            
            % Find the subspace of states that involve the current spin
            subspace_mask=(spin_system.bas.basis(:,n)~=0);
            
            % Get the part of the trajectory belonging to the subspace
            subspace_trajectory=traj(subspace_mask,:);
            
            % Get the norm of the trajectory
            result(n,:)=sqrt(sum(subspace_trajectory.*conj(subspace_trajectory),1));
            
        end
        
        % Create a legend
        legend_text=cell(spin_system.comp.nspins,1);
        for n=1:spin_system.comp.nspins
            if ~isempty(spin_system.comp.labels{n})

                % User-specified labels if available
                legend_text{n}=spin_system.comp.labels{n};

            else

                % Otherwise isotopes and numbers
                legend_text{n}=[spin_system.comp.isotopes{n} ' (' num2str(n) ')'];

            end
        end
        
        % Create labels
        label_text='density touching each spin';
        title_text='spin populations';

        % Compute Y axis extents
        max_val=max(result,[],'all');
        y_axis_extents=[-0.05*max_val 1.05*max_val];
        
    case 'local_each_spin'
        
        % Preallocate the norm trajectory array
        result=zeros(spin_system.comp.nspins,size(traj,2));
        
        % Loop over spins in the system
        for n=1:spin_system.comp.nspins
            
            % Find the subspace of states that are local to current spin
            subspace_mask=(spin_system.bas.basis(:,n)~=0)&...
                          (sum(spin_system.bas.basis,2)==spin_system.bas.basis(:,n));
            
            % Get the part of the trajectory belonging to the subspace
            subspace_trajectory=traj(subspace_mask,:);
            
            % Get the norm of the trajectory
            result(n,:)=sqrt(sum(subspace_trajectory.*conj(subspace_trajectory),1));
            
        end
        
        % Create a legend
        legend_text=cell(spin_system.comp.nspins,1);
        for n=1:spin_system.comp.nspins
            if ~isempty(spin_system.comp.labels{n})

                % User-specified labels if available
                legend_text{n}=spin_system.comp.labels{n};

            else

                % Otherwise isotopes and numbers
                legend_text{n}=[spin_system.comp.isotopes{n} ' (' num2str(n) ')'];
                
            end
        end
        
        % Create labels
        label_text='density local to each spin';
        title_text='spin populations';

        % Compute Y axis extents
        max_val=max(result,[],'all');
        y_axis_extents=[-0.05*max_val 1.05*max_val];
        
    case 'level_populations'
        
        % Move trajectory into the Zeeman basis set
        traj=sphten2zeeman(spin_system)*traj;
        
        % Find out the number of energy levels
        nlevels=sqrt(size(traj,1));
        
        % Preallocate population dynamics array
        result=zeros(nlevels,size(traj,2));
        
        % Extract the populations
        for n=1:size(traj,2)
            result(:,n)=real(diag(reshape(traj(:,n),[nlevels nlevels])));
        end
        
        % Create labels
        label_text='energy level populations';
        title_text='level populations';

        % Compute Y axis extents (log scale is used later)
        y_axis_extents=[spin_system.tols.subs_drop 1.05*max(result,[],'all')];

    otherwise
        
        % Complain and bomb out
        error('unknown property.');
        
end

% Set plot option
if exist('time_axis','var')
    
    % Plot with time axis supplied
    p=plot(time_axis,result');  
      
else
    
    % Plot with point axis
    p=plot(result'); kxlabel('trajectory point');   
    
end

% Plot with predictable colours
for n=1:numel(p) 
    p(n).Color=hsv2rgb([n/numel(p) 0.75 0.75]);  
end

% Log scale for level populations  
if ~ismember(property,'level_populations')
    set(gca,'yscale','log'); 
end

% Do not redraw the legend (expensive)
if exist('legend_text','var')&&isempty(gca().Legend)  
    leg_obj=legend(gca(),legend_text,'Location','NorthEast','AutoUpdate','off');  
    set(leg_obj.BoxFace,'ColorType','truecoloralpha',...
                        'ColorData',uint8([200 200 200 64]')); 
end

% Set Y axis label and axis scaling
kylabel(label_text); ylim(y_axis_extents); xlim tight;

% Set the title and the grid
ktitle(title_text); kgrid;

end

% Consistency enforcement
function grumble(spin_system,trajectory)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('trajectory analysis is only available for sphten-liouv formalism.');
end
if ~isnumeric(trajectory)
    error('trajectory should be an array of doubles.');
end
if size(trajectory,1)~=size(spin_system.bas.basis,1)
    error('trajectory dimension should match basis dimension.');
end
end

% And this is the whole shabby secret: to some men, the sight of
% an achievement is a reproach, a reminder that their own lives 
% are irrational, and that there is no loophole - no escape from
% reason and reality.
%
% Ayn Rand

