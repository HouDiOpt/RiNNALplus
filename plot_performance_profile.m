%% =======================================================================
%% Generate BIQ/QKP performance profile plots and optionally recompute .mat results.
%% Set regenerate=true to recompute BIQ_all/QKP_all before plotting.
%% =======================================================================

clear 
setup_path;
regenerate = false;
baseDir = fullfile('results','performance_profile');

if regenerate
    runs_QKP_all();
    runs_BIQ_all();
end

% BIQ
load(fullfile(baseDir,'BIQ_all.mat'));
figure;
perf(fullfile(baseDir,'table_BIQ_pp'),BIQ_all, 1);
title('Performance Profile: BIQ');

% QKP
load(fullfile(baseDir,'QKP_all.mat'));
figure;
perf(fullfile(baseDir,'table_QKP_pp'),QKP_all,1);
title('Performance Profile: QKP');

%% ========================= Local functions =========================
function runs_BIQ_all()
% Generate BIQ_all for performance profiles (RLT/DNN/COMP + SDPNAL+)
baseDir = fullfile('results','performance_profile');
useRiNNAL = 1;
useSDPNAL = 1;

probtype = 'BIQ';
MBQPtypes = [1 2 3]; % 1: RLT, 2: slack (DNN), 3: comp
tol = 1e-6;
fname = feval(strcat(['problems_' probtype]))';

idxList = 31:61;
nprob = numel(idxList);
BIQ_all = NaN(nprob,4); % [RLT, DNN, COMP, SDPNAL+]

for ii = 1:nprob
    i = idxList(ii);
    load(strcat(fname{i},'.mat'),'data');
    Q = data.Q;
    n0 = size(Q,1);
    c = data.c;
    A = data.A;
    b = data.b;
    B = [];
    d = [];
    bidx = data.bidx;
    zid  = [];

    for k = 1:numel(MBQPtypes)
        MBQPtype = MBQPtypes(k);
        [Qf,cf,Af,bf,Bf,df,bidxf,zidf] = formulate_MBQP(Q,c,A,b,B,d,bidx,zid,MBQPtype);
        if useRiNNAL
            try
                clear par
                par.tol = tol;
                r = [];
                [~,~,info_ALM] = RiNNAL_plus(Qf,cf,Af,bf,Bf,df,bidxf,r,zidf,par);
                BIQ_all(ii,k) = info_ALM.ttime;
                if BIQ_all(ii,k) > 3600
                    BIQ_all(ii,k) = NaN;
                end
            catch
                BIQ_all(ii,k) = NaN;
            end
        end
    end

    if useSDPNAL
        MBQPtype = 1;
        [Qf,cf,Af,bf,Bf,df,bidxf,zidf] = formulate_MBQP(Q,c,A,b,B,d,bidx,zid,MBQPtype);
        try
            [blk,At,Bt,C,bb,dd] = MBQP_to_SDPNAL(Qf,cf,Af,bf,Bf,df,bidxf,zidf);
            LL = 0;
            OPTIONS.tol = tol;
            OPTIONS.maxtime = 3600;
            OPTIONS.maxiter = 2e6;
            if size(At{1},2)+size(Bt{1},2) < 1000
                OPTIONS.AATsolve.method = 'direct';
            else
                OPTIONS.AATsolve.method = 'iterative';
            end
            OPTIONS.stopoption = 0;
            tic;
            sdpnalplus(blk,At,C,bb,LL,[],Bt,dd,[],OPTIONS);
            BIQ_all(ii,4) = toc;
            if BIQ_all(ii,4) > 3600
                BIQ_all(ii,4) = NaN;
            end
        catch
            BIQ_all(ii,4) = NaN;
        end
    end
end

if ~exist(baseDir,'dir')
    mkdir(baseDir);
end
save(fullfile(baseDir,'BIQ_all.mat'),'BIQ_all');
end

function runs_QKP_all()
% Generate QKP_all for performance profiles (3 formulations)
baseDir = fullfile('results','performance_profile');
useRiNNAL = 1;

probtype = 'QKP';
tol = 1e-6;
fname = feval(strcat(['problems_' probtype]))';

idxList = 1:12;
nprob = numel(idxList);
QKP_all = NaN(nprob,3); % [RLT, DNN, COMP]

for ii = 1:nprob
    i = idxList(ii);
    load(strcat(fname{i},'.mat'),'data');
    Q0 = data.Q;
    n0 = size(Q0,1);
    c0 = data.c;
    A0 = data.A;
    b0 = data.b;
    bidx0 = data.bidx;
    zid0 = [];

    % Formulation 1: Ax=b, x<=1 (SDP-RLT)
    Q = Q0; c = c0; A = A0; b = b0;
    B = []; d = [];
    bidx = bidx0; zid = zid0;
    try
        clear par
        par.tol = tol;
        r = [];
        [Qf,cf,Af,bf,Bf,df,bidxf,zidf] = formulate_MBQP(Q,c,A,b,B,d,bidx,zid,1);
        [~,~,info_ALM] = RiNNAL_plus(Qf,cf,Af,bf,Bf,df,bidxf,r,zidf,par);
        QKP_all(ii,1) = info_ALM.ttime;
        if QKP_all(ii,1) > 3600
            QKP_all(ii,1) = NaN;
        end
    catch
        QKP_all(ii,1) = NaN;
    end

    % Formulation 2: Bx+s=d (DNN slack)
    Q = Q0; c = c0; A = A0; b = b0;
    B = []; d = [];
    bidx = bidx0; zid = zid0;
    try
        clear par
        par.tol = tol;
        r = [];
        [Qf,cf,Af,bf,Bf,df,bidxf,zidf] = formulate_MBQP(Q,c,A,b,B,d,bidx,zid,2);
        [~,~,info_ALM] = RiNNAL_plus(Qf,cf,Af,bf,Bf,df,bidxf,r,zidf,par);
        QKP_all(ii,2) = info_ALM.ttime;
        if QKP_all(ii,2) > 3600
            QKP_all(ii,2) = NaN;
        end
    catch
        QKP_all(ii,2) = NaN;
    end

    % Formulation 3: Bx<=d (COMP-style)
    Q = Q0; c = c0; A = A0; b = b0;
    B = []; d = [];
    bidx = bidx0; zid = zid0;
    try
        clear par
        par.tol = tol;
        r = [];
        [Qf,cf,Af,bf,Bf,df,bidxf,zidf] = formulate_MBQP(Q,c,A,b,B,d,bidx,zid,3);
        [~,~,info_ALM] = RiNNAL_plus(Qf,cf,Af,bf,Bf,df,bidxf,r,zidf,par);
        QKP_all(ii,3) = info_ALM.ttime;
        if QKP_all(ii,3) > 3600
            QKP_all(ii,3) = NaN;
        end
    catch
        QKP_all(ii,3) = NaN;
    end
end

if ~exist(baseDir,'dir')
    mkdir(baseDir);
end
save(fullfile(baseDir,'QKP_all.mat'),'QKP_all');
end

function [Q,c,A,b,B,d,bidx,zid] = formulate_MBQP(Q,c,A,b,B,d,bidx,zid,type)
% Shared formulation helper (BIQ/QKP)
m1 = size(A,1);
m2 = size(B,1);
n  = size(Q,1);
l  = length(bidx);

if l == 0
    error('\n No binary constraint, rewrite formulation function!')
end

eyeid = eye(n);
eyeid = eyeid(bidx,:);

if type == 1
    % Add redundant bound x_B <= e
    B = eyeid;
    d = ones(l,1);
elseif type == 2
    if m2 == 0
        B = eyeid;
        d = ones(l,1);
        m2 = size(B,1);
    else
        B = [B; eyeid];
        d = [d; ones(l,1)];
        m2 = size(B,1);
    end
    Q = [Q,zeros(n,m2);zeros(m2,n+m2)];
    c = [c;zeros(m2,1)];
    A = [A,zeros(m1,m2);B,eye(m2)];
    b = [b;d];
    B = zeros(0,n+m2);
    d = zeros(0,1);
    G = zeros(n+1,n+1);
    G(zid) = 1;
    G = [G,zeros(n+1,m2);zeros(m2,n+m2+1)];
    zid = find(G);
elseif type == 3
    m2orig = m2;
    if m2 == 0
        B = eyeid;
        d = ones(l,1);
    else
        B = [B; eyeid];
        d = [d; ones(l,1)];
    end
    m2 = size(B,1);
    Q = [Q,zeros(n,m2);zeros(m2,n+m2)];
    c = [c;zeros(m2,1)];
    A = [A,zeros(m1,m2);B,eye(m2)];
    b = [b;d];
    B = zeros(0,n+m2);
    d = zeros(0,1);
    bidx = [];
    G = zeros(n+1,n+1);
    G(zid) = 1;
    addCols = zeros(n+1,m2);
    addCols(2:end,(m2orig+1):m2) = eyeid';
    bottom = zeros(m2,1+n+m2);
    bottom((m2orig+1):m2,2:(n+1)) = eyeid;
    G = [G,addCols;bottom];
    zid = find(G);
end

if nnz(Q)/numel(Q) <= 0.5
    Q = sparse(Q);
end
if nnz(A)/numel(A) <= 0.5
    A = sparse(A);
end
if nnz(B)/numel(B) <= 0.5
    B = sparse(B);
end
end
