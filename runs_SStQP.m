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
probtype = 'SStQP';
MBQPtype = 4;
%% constraint type
% 1: Big-M
% 2: Complementarity
% 3. Partial Comp (without upperbound)
% 4. Ineq (Big-M)
% 5. Ineq (Comp)
% 6. Ineq (Big-M) with e'u<= rho
rhorate = 1/2;
tol = 1e-6;
record = 1;
fname = feval(strcat(['problems_' probtype]))';

%% record results
rrALM   = {"data","n","eq","binary","rank","pobj","dobj","pfeas","dfeas","comp","pdgap",...
    "iter","itersub","time","rhorate","MBQPtype","time_PG","PGiter","SSNiter","SSNCGiter","aveCG"};
rrNAL   = {"data","n","eq","binary","rank","pobj","dobj","pfeas","dfeas","comp","pdgap","iter","itersub","iterA","time","rhorate","MBQPtype"};

%% start test
for i = [8]; [8 50 92 24 66 108 32 74];
    %% data
    %  n  |   25    |   50    |  100   |  200    |  500    |
    % ------------------------------------------------------
    % COP | 26:30   | 38:42   | 8:13   | 20:25   | 32:37   |
    % PSD | 68:72   | 80:84   | 50:55  | 62:67   | 74:79   |
    % SPN | 110:114 | 122:126 | 92:97  | 104:109 | 116:121 |

    %% real data
    Qtmp = load(strcat(fname{i},'.txt'));
    ntmp = size(Qtmp,1);
    ctmp = zeros(ntmp,1);
    rho0_match = regexp(fname{i}, 'rho_(\d+)', 'tokens');
    rho0_value = str2double(rho0_match{1}{1});
    rho = round(rhorate*rho0_value);
    [Q,c,A,b,B,d,bidx,zid,m,n,l,nc] = formulate_SparseStQP(Qtmp,ctmp,MBQPtype,rho);
    fprintf('\n ===> Start testing %s problems: %s \n',probtype,fname{i});

    %% RiNNAL+
    if useRiNNAL
        %%
        clear par
        par.tol = tol;
        r = [];
        %%
        par.projtype = 0;
        if size(Q,1)>= 500
            par.beta = 1e-3;
        end
        %%
        [obj_ALM,X_ALM,info_ALM] = RiNNAL_plus(Q,c,A,b,B,d,bidx,r,zid,par);
        if record
            nlALM = {fname{i},n,m,length(bidx),size(X_ALM.R,2),obj_ALM(1),obj_ALM(2),...
                info_ALM.pfeas,info_ALM.dfeas,info_ALM.comp,info_ALM.pdgap,...
                info_ALM.ALMite,info_ALM.BBite,info_ALM.ttime,rhorate,MBQPtype,...
                info_ALM.ttime_PG,info_ALM.PGiterTotal,info_ALM.SSNTotal,info_ALM.CGiterTotal,info_ALM.aveCGiter};
            rrALM = [rrALM;nlALM];
            cname = strcat(pwd,"/results/",probtype,'/',probtype,"-ALM-",string(datetime('now','Format','yyyy-MM-dd--HH-mm-ss')),".mat");
            save(cname,"rrALM");
        end
    end

    %% SDPNAL
    if useSDPNAL
        [blk,At,Bt,C,bb,dd] = MBQP_to_SDPNAL(Q,c,A,b,B,d,bidx,zid);
        LL = 0;
        OPTIONS.tol = tol;
        OPTIONS.maxtime = 3600;
        if size(At{1},2)+size(Bt{1},2) < 1000
            OPTIONS.AATsolve.method = 'direct';
        else
            OPTIONS.AATsolve.method = 'iterative';
        end
        OPTIONS.stopoption = 0;
        OPTIONS.maxiter = 1e5;
        tic;
        [obj_NAL,X_NAL,~,~,~,~,~,~,info_NAL,~] = ...
            sdpnalplus(blk,At,C,bb,LL,[],Bt,dd,[],OPTIONS);
        ttime_NAL = toc;
        fprintf('\n SDPNAL+ rank  = %3.0d\n',rank(X_NAL{1},1e-6));
        if record
            nlNAL = {fname{i},n,m,length(bidx),rank(X_NAL{1},1e-6),obj_NAL(1),obj_NAL(2),...
                max([info_NAL.etaRp,info_NAL.etaK1,info_NAL.etaK2]),...
                info_NAL.etaRd,max([info_NAL.etaC1,info_NAL.etaC1]),...
                info_NAL.relgap,info_NAL.iterSSN,info_NAL.iterSSNsub,info_NAL.iterADM,ttime_NAL,rhorate,MBQPtype};
            rrNAL = [rrNAL;nlNAL];
            cname = strcat(pwd,"/results/",probtype,'/',probtype,"-NAL-",string(datetime('now','Format','yyyy-MM-dd--HH-mm-ss')),".mat");
            save(cname,"rrNAL");
        end
    end
end


function [Q,c,A,b,B,d,bidx,zid,m,n,l,nc] = formulate_SparseStQP(Qtmp,ctmp,MBQPtype,rho)

n = size(Qtmp,1);
en = ones(n,1);
eyen = eye(n);
zn = zeros(n,n);

if MBQPtype == 1
    Q = zeros(4*n,4*n);
    Q(1:n,1:n) = Qtmp;
    c = [ctmp;zeros(3*n,1)];
    A = [en',zeros(1,3*n); ...
        zeros(1,n),en',zeros(1,2*n); ...
        eyen,-eyen,eyen,zn; ...
        zn,eyen,zn,eyen];
    b = [1;rho;zeros(n,1);en];
    bidx = (n+1):(2*n);
    zid = [];
    B = [];
    d = [];
elseif MBQPtype == 2
    Q = zeros(3*n,3*n);
    Q(1:n,1:n) = Qtmp;
    c = [ctmp;zeros(2*n,1)];
    A = [en',zeros(1,2*n); ...
        zeros(1,n),en',zeros(1,1*n); ...
        zn,eyen,eyen];
    b = [1;n-rho;en];
    bidx = (n+1):(2*n);
    G = zeros(3*n,3*n);
    G(1:n,n+1:2*n) = eyen;
    G(n+1:2*n,1:n) = eyen;
    M = [zeros(1,3*n+1);zeros(3*n,1),G];
    zid = find(M);
    B = [];
    d = [];
elseif MBQPtype == 3
    Q = zeros(2*n,2*n);
    Q(1:n,1:n) = Qtmp;
    c = [ctmp;zeros(1*n,1)];
    A = [en',zeros(1,1*n); ...
        zeros(1,n),en'];
    b = [1;n-rho];
    bidx = (n+1):(2*n);
    G = zeros(2*n,2*n);
    G(1:n,n+1:2*n) = eyen;
    G(n+1:2*n,1:n) = eyen;
    M = [zeros(1,2*n+1);zeros(2*n,1),G];
    zid = find(M);
    B = [];
    d = [];
elseif MBQPtype == 4
    Q = zeros(2*n,2*n);
    Q(1:n,1:n) = Qtmp;
    c = [ctmp;zeros(1*n,1)];
    A = [en',zeros(1,n); ...
        zeros(1,n),en'];
    b = [1;rho];
    bidx = (n+1):(2*n);
    zid = [];
    B = [zeros(n,n),eyen; ...
        eyen,-eyen];
    d = [ones(n,1);zeros(n,1)];
elseif MBQPtype == 5
    Q = zeros(2*n,2*n);
    Q(1:n,1:n) = Qtmp;
    c = [ctmp;zeros(1*n,1)];
    A = [en',zeros(1,n); ...
        zeros(1,n),en'];
    b = [1;n-rho];
    bidx = (n+1):(2*n);
    G = zeros(2*n,2*n);
    G(1:n,n+1:2*n) = eyen;
    G(n+1:2*n,1:n) = eyen;
    M = [zeros(1,2*n+1);zeros(2*n,1),G];
    zid = find(M);
    B = [zeros(n,n),eyen];
    d = [ones(n,1)];
elseif MBQPtype == 6
    Q = zeros(2*n,2*n);
    Q(1:n,1:n) = Qtmp;
    c = [ctmp;zeros(1*n,1)];
    A = [en',zeros(1,n)];
    b = [1];
    bidx = (n+1):(2*n);
    zid = [];
    B = [zeros(n,n),eyen; ...
        eyen,-eyen; ...
        zeros(1,n),en'];
    d = [ones(n,1);zeros(n,1);rho];
else
    fprintf('Type Error!\n');
end

A = sparse(A);
B = sparse(B);
Q = sparse(Q);
m = size(A,1);
n = size(Q,1);
l = length(bidx);
nc = m+m*n+l+1;

end