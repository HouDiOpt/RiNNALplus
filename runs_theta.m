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
probtype = 'Theta';
MBQPtype = 1; % 0. no strengthen 1. (<=) 2. (=)slack 3. (=) comp
tol = 1e-6;
record = 1;
fname = feval(strcat(['problems_' probtype]))';

%% record results
rrALM = {"data","n","eq","binary","rank","pobj","dobj","pfeas","dfeas","comp","pdgap","iter","itersub","time","iq",...
    "time_PG","PGiter","SSNiter","SSNCGiter","aveCG"};
rrNAL = {"data","n","eq","binary","rank","pobj","dobj","pfeas","dfeas","comp","pdgap","iter","itersub","iterA","time","iq"};

%% start test
for i = [1]; [1:55 105:106 111:112 117:118 122:123 125];
    % Gset  : 1:55    (low    rank)
    % Coding: 101:132 (medium rank)

    %% data
    load(strcat('Theta-',fname{i},'.mat'),'data');
    Q = data.Q;
    c = data.c;
    A = data.A;
    b = data.b;
    B = [];
    d = [];
    bidx = data.bidx;
    zid = data.zeroidvec;
    [Q,c,A,b,B,d,bidx,zid,n,l,m1,m2,neq,niq] = formulate_MBQP(Q,c,A,b,B,d,bidx,zid,MBQPtype);
    fprintf('\n ===> Start testing %s problems: %s \n',probtype,fname{i});

    %% RiNNAL+
    if useRiNNAL
        %%
        clear par
        par.tol = tol;
        r = [];
        %%
        [obj_ALM,X_ALM,info_ALM] = RiNNAL_plus(Q,c,A,b,B,d,bidx,r,zid,par);
        if record
            nlALM = {fname{i},n,m1,length(bidx),size(X_ALM.R,2),obj_ALM(1),obj_ALM(2),...
                info_ALM.pfeas,info_ALM.dfeas,info_ALM.comp,info_ALM.pdgap,...
                info_ALM.ALMite,info_ALM.BBite,info_ALM.ttime,m2,...
                info_ALM.ttime_PG,info_ALM.PGiterTotal,info_ALM.SSNTotal,info_ALM.CGiterTotal,info_ALM.aveCGiter};
            rrALM = [rrALM;nlALM];
            cname = strcat(pwd,"/results/",probtype,'/',probtype,"-ALM-",string(datetime('now','Format','yyyy-MM-dd--HH-mm-ss')),".mat");
            save(cname,"rrALM");
        end
    end

    %% SDPNAL
    if useSDPNAL
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
            nlNAL = {fname{i},n,m1,length(bidx),rank(X_NAL{1},1e-6),obj_NAL(1),obj_NAL(2),...
                max([info_NAL.etaRp,info_NAL.etaK1,info_NAL.etaK2]),...
                info_NAL.etaRd,max([info_NAL.etaC1,info_NAL.etaC1]),...
                info_NAL.relgap,info_NAL.iterSSN,info_NAL.iterSSNsub,info_NAL.iterADM,ttime_NAL,m2};
            rrNAL = [rrNAL;nlNAL];
            cname = strcat(pwd,"/results/",probtype,'/',probtype,"-NAL-",string(datetime('now','Format','yyyy-MM-dd--HH-mm-ss')),".mat");
            save(cname,"rrNAL");
        end
    end
end

%%
function [Q,c,A,b,B,d,bidx,zid,n,l,m1,m2,neq,niq] = formulate_MBQP(Q,c,A,b,B,d,bidx,zid,type)

m1 = size(A,1);
m2 = size(B,1);
n  = size(Q,1);
l  = length(bidx);

if m2 ~= 0
    error('\n Have inequality constraint, rewrite formulation function!');
elseif l == 0
    error('\n No binary constraint, rewrite formulation function!')
end

eyeid = eye(n);
eyeid = eyeid(bidx,:);

if type == 0

elseif type == 1
    B = eyeid;
    d = [ones(l,1)];
elseif type == 2
    Q = [Q,zeros(n,l);zeros(l,n+l)];
    c = [c;zeros(l,1)];
    A = [A,zeros(m1,l);eyeid,eye(l)];
    b = [b;ones(l,1)];
    G = zeros(n+1,n+1);
    G(zid) = 1;
    G = [G,zeros(n+1,l);zeros(l,n+l+1)];
    zid = find(G);
elseif type == 3
    Q = [Q,zeros(n,l);zeros(l,n+l)];
    c = [c;zeros(l,1)];
    A = [A,zeros(m1,l);eyeid,eye(l)];
    b = [b;ones(l,1)];
    bidx = [];
    G = zeros(n+1,n+1);
    G(zid) = 1;
    G = [G,[zeros(1,l);eyeid'];zeros(l,1),eyeid,zeros(l,l)];
    zid = find(G);
end

m1 = size(A,1);
m2 = size(B,1);
n = size(Q,1);
l = length(bidx);
neq = m1+m1*n+l+1;
niq = m2*(2*n+m2+1)/2;

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
