function [] = generate_data_BIQ()
rng('default');
probtype = 'BIQ';
nk = [5000 10000 15000 20000];
for k = 1:length(nk)
    %% tuning par
    n = nk(k);
    m = 0;
    %% fixed par
    bidx = 1:n;
    Q = (randi(200,n,n)-100).*(rand(n)<0.1);
    Q = triu(Q)+triu(Q,1)';
    Q = -Q;
    % Q = sparse(Q);
    c = zeros(n,1);
    A = zeros(0,n);
    b = zeros(0,1); 
    %% record data
    data.problemtype = probtype;
    data.n = n;
    data.m = m;
    data.bidx = bidx;
    data.Q = Q;
    data.c = c;
    data.A = A;
    data.b = b;
    %% save data
    probname = strcat('BIQ-','n',string(n),'.mat');
    save(probname,'data','-mat');
end
end