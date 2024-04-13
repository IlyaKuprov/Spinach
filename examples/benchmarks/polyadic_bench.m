% A benchmark for the polyadic object.
%
% i.kuprov@soton.ac.uk

function polyadic_bench()

% Statistics parameters
n_mats=3;       % number of matrices in the polyadic
max_dim=32;     % max dimension of each matrix
n_stats=100;    % statistics sample size
nnz_per_col=5;  % nnz per column for sparse matrices

%% Full matrix benchmark

% Result array
runtimes_full_poly=zeros(n_stats,1);
runtimes_full_flat=zeros(n_stats,1);

% Full matrix statistics loop
for n=1:n_stats
    
    % Update the user
    disp(['stat loop instance ' num2str(n) '/' num2str(n_stats) '...']);
    
    % Get random full complex matrices
    terms=cell(1,n_mats);
    for k=1:n_mats
        dim=randi(max_dim,1);
        terms{k}=1i*randn(dim,dim)+...
                    randn(dim,dim);
    end
    
    % Form a polyadic
    P=polyadic({terms});
    
    % Get a random full complex vector
    v=1i*randn(size(P,2),1)+...
         randn(size(P,2),1);
     
    % Time polyadic multiplication
    tic; P*v; runtimes_full_poly(n)=1000*toc; %#ok<VUNUS>
    
    % Inflate the polyadic
    P=full(P);
    
    % Time flat multiplication
    tic; P*v; runtimes_full_flat(n)=1000*toc; %#ok<VUNUS>
 
end

%% Sparse matrix benchmark

% Result array
runtimes_sparse_poly=zeros(n_stats,1);
runtimes_sparse_flat=zeros(n_stats,1);

% Sparse matrix statistics loop
for n=1:n_stats
    
    % Update the user
    disp(['stat loop instance ' num2str(n) '/' num2str(n_stats) '...']);
    
    % Get random sparse complex matrices
    terms=cell(1,n_mats);
    for k=1:n_mats
        dim=randi(max_dim,1);
        terms{k}=exp(1i/3)*sprandn(dim,dim,nnz_per_col/dim);
    end
    
    % Form a polyadic
    P=polyadic({terms});
    
    % Get a random full complex vector
    v=1i*randn(size(P,2),1)+...
         randn(size(P,2),1);
     
    % Time polyadic multiplication
    tic; P*v; runtimes_sparse_poly(n)=1000*toc; %#ok<VUNUS>
    
    % Inflate the polyadic
    P=inflate(P);
    
    % Time flat multiplication
    tic; P*v; runtimes_sparse_flat(n)=1000*toc; %#ok<VUNUS>
        
end

%% Print statistics

disp(['Average polyadic, full: ' num2str(mean(runtimes_full_poly)), ...
      ', stdev: ' num2str(std(runtimes_full_poly)/sqrt(n_stats))]);
disp(['Average standard, full: ' num2str(mean(runtimes_full_flat)), ...
      ', stdev: ' num2str(std(runtimes_full_flat)/sqrt(n_stats))]);
disp(['Average polyadic, sparse: ' num2str(mean(runtimes_sparse_poly)), ...
      ', stdev: ' num2str(std(runtimes_sparse_poly)/sqrt(n_stats))]);
disp(['Average standard, sparse: ' num2str(mean(runtimes_sparse_flat)), ...
      ', stdev: ' num2str(std(runtimes_sparse_flat)/sqrt(n_stats))]);
        
end

