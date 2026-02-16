function [product_table_left_ortho,product_table_right_ortho]= bos_product_table_ortho(nlevels)

% Check consistency
grumble(nlevels);

% Generate cache record name
own_path=mfilename('fullpath');
own_path=own_path(1:(end-23));
table_file=[own_path 'bos_product_table_ortho_' ...
            num2str(nlevels) '.mat'];

% Check the cache
if exist(table_file,'file')
    
    % Lift data from the cache if the file is already available
    load(table_file,'product_table_left_ortho','product_table_right_ortho');
    
else
    
    % Get orthogonal bosonic monomials
    B = boson_mono_ortho(nlevels);

    % Precompute norms
    norms=zeros(nlevels^2,1);
    for n=1:nlevels^2
        norms(n)=sqrt(hdot(B{n},B{n}));
    end
    
    % Preallocate IST product tables
    product_table_left_ortho= zeros(nlevels^2,nlevels^2,nlevels^2);
    product_table_right_ortho=zeros(nlevels^2,nlevels^2,nlevels^2);

    % Populate product tables
    for k=1:nlevels^2
        for m=1:nlevels^2
            for n=1:nlevels^2

                % Left product action: carefully tiptoe around extreme norms
                product_table_left_ortho(n,m,k)= norms(n)*hdot((B{k}/norms(k)),...
                                                        (B{n}/norms(n))*...
                                                        (B{m}/norms(m)));

                % Right product action: carefully tiptoe around extreme norms
                product_table_right_ortho(n,m,k)=norms(n)*hdot((B{k}/norms(k)),...
                                                        (B{m}/norms(m))*...
                                                        (B{n}/norms(n)));

            end
        end
    end

    try % Try to save a cache record, but do not insist
        save(table_file,'product_table_left_ortho',...
                        'product_table_right_ortho','-v7.3'); drawnow;
    catch
        warning('Spinach directory appears to be write-protected');
    end

end

end

% Consistency enforcement
function grumble(nlevels)
if (numel(nlevels)~=1)||(~isnumeric(nlevels))||...
   (~isreal(nlevels))||(nlevels<2)||(mod(nlevels,1)~=0)
    error('nlevels must be a real integer greater than 1.');
end
end