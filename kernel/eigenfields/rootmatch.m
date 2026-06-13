% Global order-preserving root matching between three magnetic field root
% lists. The routine returns indices into the input lists that identify
% matching roots with minimum field-continuation cost. Syntax:
%
%           [idx1,idx2,idx3]=rootmatch(field1,field2,field3,...
%                                      edge12,edge23,edge31)
%
% Parameters:
%
%     field1 - real vector of roots at the first triangle vertex
%
%     field2 - real vector of roots at the second triangle vertex
%
%     field3 - real vector of roots at the third triangle vertex
%
%     edge12 - positive distance between the first and the second
%              triangle vertices
%
%     edge23 - positive distance between the second and the third
%              triangle vertices
%
%     edge31 - positive distance between the third and the first
%              triangle vertices
%
% Outputs:
%
%     idx1   - indices of matched roots in field1
%
%     idx2   - indices of matched roots in field2
%
%     idx3   - indices of matched roots in field3
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=rootmatch.m>

function [idx1,idx2,idx3]=rootmatch(field1,field2,field3,...
                                    edge12,edge23,edge31)

% Check consistency
grumble(field1,field2,field3,edge12,edge23,edge31);

% Initialise empty output
idx1=[]; idx2=[]; idx3=[];

% Nothing to match
if isempty(field1)||isempty(field2)||isempty(field3)
    return;
end

% Sort roots by field
[field1,ord1]=sort(field1(:));
[field2,ord2]=sort(field2(:));
[field3,ord3]=sort(field3(:));

% Equal lists match by root order
if (numel(field1)==numel(field2))&&(numel(field2)==numel(field3))
    idx1=ord1.';
    idx2=ord2.';
    idx3=ord3.';
    return;
end

% Get list dimensions
n1=numel(field1); n2=numel(field2); n3=numel(field3);

% Initialise dynamic programme
match_count=zeros(n1+1,n2+1,n3+1);
match_cost=inf(n1+1,n2+1,n3+1);
prev_move=zeros(n1+1,n2+1,n3+1,'uint8');
match_cost(1,1,1)=0;

% Dynamic programme over sorted root lists
for a=1:(n1+1)
    for b=1:(n2+1)
        for c=1:(n3+1)

            % Skip the initial state
            if (a==1)&&(b==1)&&(c==1), continue; end

            % Initialise the local optimum
            best_count=-Inf;
            best_cost=Inf;
            best_move=uint8(0);

            % Skip a root from the first list
            if (a>1)&&isfinite(match_cost(a-1,b,c))
                test_count=match_count(a-1,b,c);
                test_cost=match_cost(a-1,b,c);
                if (test_count>best_count)||...
                   ((test_count==best_count)&&(test_cost<best_cost))
                    best_count=test_count;
                    best_cost=test_cost;
                    best_move=uint8(1);
                end
            end

            % Skip a root from the second list
            if (b>1)&&isfinite(match_cost(a,b-1,c))
                test_count=match_count(a,b-1,c);
                test_cost=match_cost(a,b-1,c);
                if (test_count>best_count)||...
                   ((test_count==best_count)&&(test_cost<best_cost))
                    best_count=test_count;
                    best_cost=test_cost;
                    best_move=uint8(2);
                end
            end

            % Skip a root from the third list
            if (c>1)&&isfinite(match_cost(a,b,c-1))
                test_count=match_count(a,b,c-1);
                test_cost=match_cost(a,b,c-1);
                if (test_count>best_count)||...
                   ((test_count==best_count)&&(test_cost<best_cost))
                    best_count=test_count;
                    best_cost=test_cost;
                    best_move=uint8(3);
                end
            end

            % Match one root from each list
            if (a>1)&&(b>1)&&(c>1)&&...
               isfinite(match_cost(a-1,b-1,c-1))
                sheet_cost=((field1(a-1)-field2(b-1))/edge12)^2+...
                           ((field2(b-1)-field3(c-1))/edge23)^2+...
                           ((field3(c-1)-field1(a-1))/edge31)^2;
                test_count=match_count(a-1,b-1,c-1)+1;
                test_cost=match_cost(a-1,b-1,c-1)+sheet_cost;
                if (test_count>best_count)||...
                   ((test_count==best_count)&&(test_cost<best_cost))
                    best_count=test_count;
                    best_cost=test_cost;
                    best_move=uint8(4);
                end
            end

            % Store the local optimum
            match_count(a,b,c)=best_count;
            match_cost(a,b,c)=best_cost;
            prev_move(a,b,c)=best_move;

        end
    end
end

% Start traceback
a=n1+1; b=n2+1; c=n3+1;

% Trace the optimal path
while (a>1)||(b>1)||(c>1)

    % Follow the stored move
    switch prev_move(a,b,c)
        case 1
            a=a-1;
        case 2
            b=b-1;
        case 3
            c=c-1;
        case 4
            idx1(end+1)=ord1(a-1); %#ok<AGROW>
            idx2(end+1)=ord2(b-1); %#ok<AGROW>
            idx3(end+1)=ord3(c-1); %#ok<AGROW>
            a=a-1; b=b-1; c=c-1;
        otherwise
            break;
    end

end

% Restore forward order
idx1=fliplr(idx1);
idx2=fliplr(idx2);
idx3=fliplr(idx3);

end

% Consistency enforcement
function grumble(field1,field2,field3,edge12,edge23,edge31)
if (~isnumeric(field1))||(~isreal(field1))||...
   ((~isempty(field1))&&(~isvector(field1)))||...
   any(~isfinite(field1(:)))
    error('field1 must be a finite real vector.');
end
if (~isnumeric(field2))||(~isreal(field2))||...
   ((~isempty(field2))&&(~isvector(field2)))||...
   any(~isfinite(field2(:)))
    error('field2 must be a finite real vector.');
end
if (~isnumeric(field3))||(~isreal(field3))||...
   ((~isempty(field3))&&(~isvector(field3)))||...
   any(~isfinite(field3(:)))
    error('field3 must be a finite real vector.');
end
if (~isnumeric(edge12))||(~isreal(edge12))||...
   (~isscalar(edge12))||(~isfinite(edge12))||(edge12<=0)
    error('edge12 must be a finite positive real scalar.');
end
if (~isnumeric(edge23))||(~isreal(edge23))||...
   (~isscalar(edge23))||(~isfinite(edge23))||(edge23<=0)
    error('edge23 must be a finite positive real scalar.');
end
if (~isnumeric(edge31))||(~isreal(edge31))||...
   (~isscalar(edge31))||(~isfinite(edge31))||(edge31<=0)
    error('edge31 must be a finite positive real scalar.');
end
end

% Of all that is written, I love only what a man 
% has written with his own blood.
%
% Friedrich Nietzsche

