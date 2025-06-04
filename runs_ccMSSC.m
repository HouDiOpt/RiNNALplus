%% preparation
close all;
clear
lastwarn('');
warning off;
rng('default');
setup_path;

%% solver
useRiNNAL = 1;
useSDPNAL = 0; %%set to 1 if wants to run SDPNALplus

%% pars
probtype = 'ccMSSC';
MBQPtype = 1; % 0. no strengthen 1. (<=) 2. (=)slack 3. (=) comp
tol = 1e-6;
record = 1;
fname = feval(strcat(['problems_',probtype]))';

%% record results
rrALM   = {"data","n","m","binary","rank","pobj","dobj","pfeas","dfeas","comp","pdgap",...
    "iter","itersub","time","norg","korg","strengthen","time_PG","PGiter","SSNiter","SSNCGiter","aveCG"};
rrNAL   = {"data","n","m","binary","rank","pobj","dobj","pfeas","dfeas","comp","pdgap","iter","itersub","iterA","time","norg","korg","strengthen"};

%% start test
for i = [5]; [1:length(fname)];
    %% data
    TRAIN_data = load(strcat(fname{i},'_TRAIN.tsv'));
    TEST_data  = load(strcat(fname{i},'_TEST.tsv'));
    all_data   = [TRAIN_data;TEST_data];
    [Q,c,A,b,B,d,bidx,zid,m,n,l,norg,korg,c_vec,neq,niq] = formulate_cluster(all_data,MBQPtype);
    fprintf('\n ===> Start testing %s problems: %s \n',probtype,fname{i});
    %% RNNAL
    if useRiNNAL
        %% parameters
        clear par
        par.tol = tol;
        r = [];
        %%
        par.projtype = 0;
        par.beta = 1e-2;
        par.Amap1 = @(X)AmapGEP(X,norg,korg);
        par.Atmap1 = @(lam)AtmapGEP(lam,norg,korg);
        if korg == 2
            bidx = bidx(1:length(bidx)/2); % avoid nonsmoothness
        end
        %%
        [obj_ALM,X_ALM,info_ALM] = RiNNAL_plus(Q,c,A,b,B,d,bidx,r,zid,par);
        if record
            nlALM = {fname{i},n,m,length(bidx),size(X_ALM.R,2),obj_ALM(1),obj_ALM(2),...
                info_ALM.pfeas,info_ALM.dfeas,info_ALM.comp,info_ALM.pdgap,...
                info_ALM.ALMite,info_ALM.BBite,info_ALM.ttime,norg,korg,MBQPtype,...
                info_ALM.ttime_PG,info_ALM.PGiterTotal,info_ALM.SSNTotal,info_ALM.CGiterTotal,info_ALM.aveCGiter};
            rrALM = [rrALM;nlALM];
            cname = strcat(pwd,"/results/",probtype,'/',probtype,"-ALM-",string(datetime('now','Format','yyyy-MM-dd--HH-mm-ss')),".mat");
            save(cname,"rrALM");
        end
    end

    %% SDPNAL
    if useSDPNAL && n <= 1500
        [blk,At,Bt,C,bb,dd,flag] = MBQP_to_SDPNAL(Q,c,A,b,B,d,bidx,zid);
        LL = 0;
        OPTIONS.tol = tol;
        OPTIONS.maxtime = 3600;
        OPTIONS.maxiter = 1e5;
        if size(At{1},2)+size(Bt{1},2) < 1000
            OPTIONS.AATsolve.method = 'direct';
        else
            OPTIONS.AATsolve.method = 'iterative';
        end
        OPTIONS.stopoption = 0;
        tic;
        [obj_NAL,X_NAL,~,~,~,~,~,~,info_NAL,~] = ...
            sdpnalplus(blk,At,C,bb,LL,[],Bt,dd,[],OPTIONS);
        ttime_NAL = toc;
        fprintf('\n SDPNAL+ rank  = %3.0d\n',rank(X_NAL{1},1e-6));
        if record
            nlNAL = {fname{i},n,m,length(bidx),rank(X_NAL{1},1e-6),obj_NAL(1),obj_NAL(2),...
                max([info_NAL.etaRp,info_NAL.etaK1,info_NAL.etaK2]),...
                info_NAL.etaRd,max([info_NAL.etaC1,info_NAL.etaC1]),...
                info_NAL.relgap,info_NAL.iterSSN,info_NAL.iterSSNsub,info_NAL.iterADM,ttime_NAL,norg,korg,MBQPtype};
            rrNAL = [rrNAL;nlNAL];
            cname = strcat(pwd,"/results/",probtype,'/',probtype,"-NAL-",string(datetime('now','Format','yyyy-MM-dd--HH-mm-ss')),".mat");
            save(cname,"rrNAL");
        end
    end
end

%%======================================================================
%% formulate SparseQP
%%======================================================================
function [Q,c,A,b,B,d,bidx,zid,m,n,l,norg,korg,c_vec,neq,niq] = formulate_cluster(raw_data,type)

label = raw_data(:,1);
unique_integers = unique(label); % Get unique values
d = size(raw_data,2)-1;
n = size(raw_data,1);
k = numel(unique_integers); % Count the unique values
c_vec = histc(label, unique_integers);

if n*k <= 5000

    hasNaN = any(isnan(raw_data), 'all');  % 'all' checks the entire matrix
    if hasNaN
        disp('The matrix contains missing (NaN) values. fill missing by interpolation.');
        raw_data = fillmissing(raw_data,'linear',2,'EndValues','nearest');
    else
        % disp('The matrix does not contain any missing (NaN) values.');
    end

    A = generate_ot_matrix(n,k,1); % 1 for sparse, (col sum) + (row sum-1)
    b = [c_vec;ones(n-1,1)];
    RD = raw_data(:,2:end);
    D = sum(RD.^2, 2) + sum(RD.^2, 2)' - 2 * (RD * RD');
    Q = kron(spdiags(1./c_vec,0,k,k),D);

    norg = n;
    korg = k;
    nk = n*k;

    if type == 0
        n = n*k;
        c = zeros(n,1);
        bidx = 1:n;
        B = [];
        d = [];
        zid = [];
    elseif type == 1
        n = n*k;
        c = zeros(n,1);
        bidx = 1:n;
        B = speye(n);
        d = ones(n,1);
        zid = [];
    elseif type == 2
        Q = [Q,zeros(n*k,n*k);zeros(n*k,2*n*k)];
        A = [A,zeros(n+k-1,n*k);speye(n*k),speye(n*k)];
        b = [b;ones(n*k,1)];
        B = [];
        d = [];
        n = n*k+n*k;
        c = zeros(n,1);
        bidx = 1:(norg*k);
        zid = [];
    elseif type == 3
        Q = [Q,zeros(n*k,n*k);zeros(n*k,2*n*k)];
        A = [A,zeros(n+k-1,n*k);speye(n*k),speye(n*k)];
        b = [b;ones(n*k,1)];
        B = [];
        d = [];
        %
        G = zeros(2*nk,2*nk);
        G(1:nk,nk+1:2*nk) = eye(nk);
        G(nk+1:2*nk,1:nk) = eye(nk);
        M = [zeros(1,2*nk+1);zeros(2*nk,1),G];
        zid = find(M);
        %
        n = n*k+n*k;
        c = zeros(n,1);
        bidx = [];
    end

    A = sparse(A);
    B = sparse(B);
    Q = sparse(Q);
    m = size(A,1);
    l = length(bidx);
    nc = m+m*n+l+1;
else
    Q = [];
    c = [];
    A = [];
    B = [];
    b = [];
    d = [];
    bidx = [];
    zid = [];
    m = 0;
    n = 0;
    l = 0;
    nc = 0;
    norg = 0;
    korg = 0;
    c_vec = [];
end

p = size(B,1);
neq = m+m*n+length(bidx)+1+length(zid)/2;
niq = p*n+p*(p+1)/2;

end


