%% =======================================================================
%% Compare tightness of BIQ relaxations (requires Gurobi/SDPNAL+).
%% =======================================================================

%% preparation
close all;
clear
lastwarn('');
warning off;
rng('default');
setup_path; 
%%addpath(genpath('/gurobi1003/linux64'));
%%addpath(genpath('/Library/gurobi1001/macos_universal2/matlab'));
addpath(genpath('C:/gurobi1103/win64'));

%% solver
useSDPNAL = 1;
useGurobi = 1;

%% pars
probtype = 'BIQ';
MBQPtypes = 0:5;
% -------------------------
%          | BIQ  | BIQ-S |
% -------------------------
%  SHOR    |  0   |   5   |
%  SDP-RLT |  4   |   1   |
%  DNN     |  4   |   2   |
%  COMP    |  \   |   3   |
% -------------------------
tol = 1e-10;
record = 1;
fname = feval(strcat(['problems_' probtype]))';

%% record results
rrALM   = {"data","n","eq","binary","rank","pobj","dobj","pfeas","dfeas","comp","pdgap","iter","itersub","time","iq"};
rrCGAL  = {"data","n","eq","binary","rank","pobj","dobj","pfeas","dfeas","comp","pdgap","iter","itersub","time","iq"};
rrNAL   = {"data","n","eq","binary","rank","pobj","dobj","pfeas","dfeas","comp","pdgap","iter","itersub","iterA","time","iq"};

%% start test
for i = 11

    %% data
    load(strcat(fname{i},'.mat'),'data');
    Q = data.Q;
    c = data.c;
    A = data.A;
    b = data.b;
    B = [];
    d = [];
    bidx = data.bidx;
    zid  = [];
    % store lower bounds for each MBQPtype
    v_lb = NaN(1,6);
    v_star = NaN;
    [Q0,c0,A0,b0,B0,d0,bidx0,zid0,n,l,m1,m2,neq,niq] = formulate_MBQP(Q,c,A,b,B,d,bidx,zid,0);
    fprintf('\n ===> Start testing %s problems: %s \n',probtype,fname{i});

    %% Gurobi
    if useGurobi == 1
        vt = cell(n,1);
        vt(:) = {'C'};
        vt(bidx0)={'B'};
        vt = char(vt);
        model.Q = sparse(Q0);
        model.A = [sparse(A0);sparse(B0)];
        model.obj = 2*c0';
        model.rhs = [b0;d0];
        model.sense = [repmat('=',size(A0,1),1);repmat('<',size(B0,1),1)];
        model.vtype = vt;
        model.lb = zeros(n,1);
        params.NonConvex = 2;
        results = gurobi(model,params);
        v_star = results.objval;
    end

    %% SDPNAL for all MBQPtype choices
    if useSDPNAL
        for MBQPtype = MBQPtypes
            [Q,c,A,b,B,d,bidx,zid,n,l,m1,m2,neq,niq] = formulate_MBQP(Q0,c0,A0,b0,B0,d0,bidx0,zid0,MBQPtype);
            % convert data
            if MBQPtype == 0 || MBQPtype == 5
                [blk,At,Bt,C,bb,dd] = MBQP_to_SDPNAL_test_strength(Q,c,A,b,B,d,bidx,zid);
                LL = [];
            elseif MBQPtype == 1 || MBQPtype == 2 || MBQPtype == 3
                [blk,At,Bt,C,bb,dd] = MBQP_to_SDPNAL(Q,c,A,b,B,d,bidx,zid);
                LL = 0;
            elseif MBQPtype == 4
                [blk,At,Bt,C,bb,dd] = MBQP_to_SDPNAL_test_strength(Q,c,A,b,B,d,bidx,zid);
                LL = 0;
            end
            % solver
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
            fprintf('\n SDPNAL+ (type %d) rank  = %3.0d\n',MBQPtype,rank(X_NAL{1},1e-6));
            v_lb(MBQPtype+1) = obj_NAL(2);
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

    %% print summary table
    methods = {'SHOR','SDP-RLT','DNN','COMP'};
    biq_idx = [0 4 4 NaN];
    biqs_idx = [5 1 2 3];
    fprintf('\n%12s | %12s | %7s | %12s | %7s\n','method','BIQ v','%gap','BIQ-S v','%gap');
    for r = 1:numel(methods)
        if isnan(biq_idx(r))
            v_biq = NaN;
        else
            v_biq = v_lb(biq_idx(r)+1);
        end
        v_biqs = v_lb(biqs_idx(r)+1);
        if ~isnan(v_star)
            gap_biq = (v_star - v_biq) / max(1,abs(v_star)) * 100;
            gap_biqs = (v_star - v_biqs) / max(1,abs(v_star)) * 100;
        else
            gap_biq = NaN;
            gap_biqs = NaN;
        end
        fprintf('%12s | %12.4e | %7.3f | %12.4e | %7.3f\n',methods{r},v_biq,gap_biq,v_biqs,gap_biqs);
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

if type == 0 || type == 4

elseif type == 1 || type == 5
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
