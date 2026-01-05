function [] = generate_data_Theta()
rng('default');
probtype = 'Theta';
fname = problems_Theta()';
    
%% Gset
for i = 1:100
    %%
    dataname = fname{i};
    G = readmatrix([dataname,'.txt']);
    n = G(1,1);
    edgenum = G(1,2);
    G = G(2:end,1:2);
    Q = -speye(n);
    m = 0;
    %%
    M = sparse(G(1:end,1),G(1:end,2),ones(size(G,1),1),n,n);
    M = M + M';
    Mp = [zeros(1,n+1);zeros(n,1),M];
    zeroidvec = find(Mp);
    %% fixed par
    bidx = 1:n;
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
    data.edge = edgenum;
    data.zeroidvec = zeroidvec;
    %% save data
    probname = strcat('Theta-',string(dataname),'.mat');
    save(probname,'data','-mat');
end
end