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
probtype = 'QKP';
tol = 1e-6;
record = 1;
fname = feval(strcat(['problems_' probtype]))';

%% record results
rrALM = {"data","bscale","Qdensity","n","m","binary","rank","pobj","dobj","pfeas","dfeas","comp","pdgap","iter","itersub","time",...
    "time_PG","PGiter","SSNiter","SSNCGiter","aveCG"};
rrNAL = {"data","bscale","Qdensity","n","m","binary","rank","pobj","dobj","pfeas","dfeas","comp","pdgap","iter","itersub","iterA","time"};

%% start test
for i = [1]; [1:12];

    %% data
    load(strcat(fname{i},'.mat'),'data');
    m = data.m;
    n = data.n;
    Q = data.Q;
    c = data.c;
    A = data.A;
    b = data.b;
    B = [];
    d = [];
    %% x<=1
    B = speye(n);
    d = ones(n,1);
    %%
    bidx = data.bidx;
    zid  = [];
    bs = data.bscale;
    Qd = data.Qdensity;
    l = length(bidx);
    nc = m+m*n+l+1;
    fprintf('\n Start testing %s problems: n=%2d, m=%2d, binary=%2d\n',probtype,n,m,l);

    %% RiNNAL+
    if useRiNNAL
        %%
        clear par
        par.tol = tol;
        r = [];
        %%
        [obj_ALM,X_ALM,info_ALM] = RiNNAL_plus(Q,c,A,b,B,d,bidx,r,zid,par);
        if record
            nlALM = {fname{i},bs,Qd,n,m,length(bidx),size(X_ALM.R,2),obj_ALM(1),obj_ALM(2),...
                info_ALM.pfeas,info_ALM.dfeas,info_ALM.comp,info_ALM.pdgap,...
                info_ALM.ALMite,info_ALM.BBite,info_ALM.ttime,...
                info_ALM.ttime_PG,info_ALM.PGiterTotal,info_ALM.SSNTotal,info_ALM.CGiterTotal,info_ALM.aveCGiter};
            rrALM = [rrALM;nlALM];
            cname = strcat(pwd,"/results/QKP/",probtype,"-ALM-",string(datetime('now','Format','yyyy-MM-dd--HH-mm-ss')),".mat");
            save(cname,"rrALM");
        end
    end

    %% SDPNAL
    if useSDPNAL && n<2000
        [blk,At,Bt,C,bb,dd] = MBQP_to_SDPNAL(Q,c,A,b,B,d,bidx,zid);
        LL = 0;
        OPTIONS.tol = tol;
        OPTIONS.maxtime = 3600;
        OPTIONS.maxiter = 2e6;
        if nc < 1000
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
            nlNAL = {fname{i},bs,Qd,n,m,length(bidx),rank(X_NAL{1},1e-6),obj_NAL(1),obj_NAL(2),...
                max([info_NAL.etaRp,info_NAL.etaK1,info_NAL.etaK2]),...
                info_NAL.etaRd,max([info_NAL.etaC1,info_NAL.etaC1]),...
                info_NAL.relgap,info_NAL.iterSSN,info_NAL.iterSSNsub,info_NAL.iterADM,ttime_NAL};
            rrNAL = [rrNAL;nlNAL];
            cname = strcat(pwd,"/results/QKP/",probtype,"-NAL-",string(datetime('now','Format','yyyy-MM-dd--HH-mm-ss')),".mat");
            save(cname,"rrNAL");
        end
    end
end
