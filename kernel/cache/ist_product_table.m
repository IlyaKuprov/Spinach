% Structure coefficient tables for the associative envelopes of su(mult)
% algebras. Syntax:
%
%    [product_table_left,product_table_right]=ist_product_table(mult)
%
% The input parameter is the multiplicity of the spin in question. Output
% contains the structure coefficients in the following conventions:
%
%  product_table_left(n,m,k)=trace(T{n}*T{m}*T{k}')/...
%                            sqrt(trace(T{k}*T{k}')*trace(T{m}*T{m}'))
%    
%  product_table_right(n,m,k)=trace(T{m}*T{n}*T{k}')/...
%                             sqrt(trace(T{k}*T{k}')*trace(T{m}*T{m}'));
%
% corresponding to Eq 7.18 of the first edition of IK's book (normalisati-
% on is missing in the book, that's a typo).
%
% These are expensive tables: disk cache is created and used automatically.
%
% hannah.hogben@chem.ox.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ist_product_table.m>

function [product_table_left,product_table_right]=ist_product_table(mult)

% Check consistency
grumble(mult);

% Hard-coded for spin-1/2 to avoid expensive disk hit
if mult==2
    product_table_left=zeros(4,4,4);
    product_table_left(:,:,1)= [1.0     0     0     0
                                  0     0     0  -0.5
                                  0     0   0.5     0
                                  0  -0.5     0     0];
    product_table_left(:,:,2)= [  0   1.0     0     0
                                0.5     0  -0.5     0
                                  0   0.5     0     0
                                  0     0     0     0];
    product_table_left(:,:,3)= [  0     0   1.0     0
                                  0     0     0  -0.5
                                0.5     0     0     0
                                  0   0.5     0     0];
    product_table_left(:,:,4)= [  0     0     0   1.0
                                  0     0     0     0
                                  0     0     0  -0.5
                                0.5     0   0.5     0];
    product_table_right=zeros(4,4,4);
    product_table_right(:,:,1)=[1.0     0     0     0
                                  0     0     0  -0.5
                                  0     0   0.5     0
                                  0  -0.5     0     0];
    product_table_right(:,:,2)=[  0   1.0     0     0
                                0.5     0   0.5     0
                                  0  -0.5     0     0
                                  0     0     0     0];
    product_table_right(:,:,3)=[  0     0   1.0     0
                                  0     0     0   0.5
                                0.5     0     0     0
                                  0  -0.5     0     0];
    product_table_right(:,:,4)=[  0     0     0   1.0
                                  0     0     0     0
                                  0     0     0   0.5
                                0.5     0  -0.5     0]; return; 
end

% Generate cache record name
own_path=mfilename('fullpath');
own_path=own_path(1:(end-17));
table_file=[own_path 'ist_product_table_' ...
            num2str(mult) '.mat'];

% Check the cache
if exist(table_file,'file')
    
    % Lift data from the cache if the file is already available
    load(table_file,'product_table_left','product_table_right');
    
else
    
    % Get the irreducible spherical tensors
    T=irr_sph_ten(mult);
    
    % Preallocate the arrays
    product_table_left=zeros(mult^2,mult^2,mult^2);
    product_table_right=zeros(mult^2,mult^2,mult^2);
    
    % Get the structure coefficients
    for k=1:mult^2
        for m=1:mult^2
            normalization=sqrt(trace(T{k}*T{k}')*trace(T{m}*T{m}'));
            for n=1:mult^2
                product_table_left(n,m,k)=trace(T{n}*T{m}*T{k}')/normalization;
                product_table_right(n,m,k)=trace(T{m}*T{n}*T{k}')/normalization;
            end
        end
    end

    try % Try to save a cache record, but don't insist
        save(table_file,'product_table_left',...
                        'product_table_right'); drawnow;
    catch
        warning('Spinach installation appears to be write-protected');
    end

end

end

% Consistency enforcement
function grumble(mult)
if (numel(mult)~=1)||(~isnumeric(mult))||(~isreal(mult))||...
   (mult<2)||(mod(mult,1)~=0)
    error('mult must be a real integer greater than 1.');
end
end

% According to Oxford Chemistry folklore, Peter Atkins has once asked the
% following question at an interview for a Lecturer post: "What is it that
% you have done that a technician would not do?". The candidate produced a
% reasonable answer to the effect that a technician would not have the re-
% quired skills. A better answer was suggested by a member of the Inter-
% view Board some time later: "And what is it that you, Atkins, have done
% that a journalist would not do?"

