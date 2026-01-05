function [] = generate_data_QKP()
rng('default');
probtype = 'QKP';
nk = [10000];%[500 1000 2000 5000];
beta = 0.9;
pl = [0.1 0.5 0.9];
for k = 1:length(nk)
    for l = 1:length(pl)
        %% tuning par
        n = nk(k);
        p = pl(l);
        m = 1;
        %% fixed par
        bidx = 1:n;
        Q = randi(100,n,n).*(rand(n)<p);
        Q = triu(Q)+triu(Q,1)';
        Q = -Q;
        % Q = sparse(Q);
        c = zeros(n,1);
        A = randi(50,m,n);
        b = max(50,floor(beta*sum(A,'all'))); %one constraint
        %% record data
        data.problemtype = probtype;
        data.n = n;
        data.m = m;
        data.bidx = bidx;
        data.Q = Q;
        data.c = c;
        data.A = A;
        data.b = b;
        %
        data.bscale = beta;
        data.Qdensity = p;
        %% save data
        probname = strcat('QKP-','n',string(n),'-p0',string(p*10),'-beta0',string(beta*10),'.mat');
        save(probname,'data','-mat');
    end
end
end