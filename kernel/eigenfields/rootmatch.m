% Global order-preserving root matching
function [idx1,idx2,idx3]=rootmatch(field1,field2,field3,...
                                    edge12,edge23,edge31)

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

% Of all that is written, I love only what a man 
% has written with his own blood.
%
% Friedrich Nietzsche

