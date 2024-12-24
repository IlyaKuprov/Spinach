% A simple 1D bin packing algorithm. Collects the list of  numbers
% supplied into sublists that sum to the number that is smaller or
% equal to the number specified. The algorithm is not optimal, but
% it does the job. Syntax:
%
%                  bins=binpack(box_sizes,bin_size)
%
% Parameters: 
%
%    box_sizes - a row vector of box sizes
%
%    bin_size  - an integer specifying the bin size
%
% Outputs:
%
%    bins      - a cell array of index vectors specifying
%                boxes allocated into each bin
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=binpack.m>

function bins=binpack(box_sizes,bin_size)

% Check consistency
grumble(box_sizes,bin_size);

% Number the boxes
box_index=(1:numel(box_sizes))';

% Find boxes that are bigger than bins
big_boxes=(box_sizes>bin_size); 
bins=num2cell(box_index(big_boxes));
box_sizes(big_boxes)=[];
box_index(big_boxes)=[];

% Pack the rest of the boxes
while numel(box_index)>0
    
    % Get enough boxes to fill the bin
    current_boxes=(cumsum(box_sizes)<=bin_size);

    % Stuff them into the bin
    bins{end+1}=box_index(current_boxes); %#ok<AGROW>
    box_sizes(current_boxes)=[];
    box_index(current_boxes)=[];
    
end

end

% Consistency enforcement
function grumble(box_sizes,bin_size)
if (~isnumeric(box_sizes))||(~isreal(box_sizes))||(size(box_sizes,1)~=1)||...
     any(box_sizes<1)||any(~isfinite(box_sizes))||any(mod(box_sizes,1))
    error('box_sizes must be a row vector of positive real integers.');
end
if (~isnumeric(bin_size))||(~isreal(bin_size))||(~isfinite(bin_size))||...
   (numel(bin_size)~=1)||(bin_size<1)||mod(bin_size,1)
    error('bin_size must be a positive real integer.');
end
end

% "Ford!" he said, "there's an infinite number of monkeys
% outside who want to talk to us about this script for 
% Hamlet they've worked out."
%
% Douglas Adams, "The Hitchhiker's Guide to the Galaxy"

