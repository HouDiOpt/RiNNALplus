%% =======================================================================
%% Run full QMSTP experiments and record results (optional SDPNAL+ baseline).
%% =======================================================================

%% preparation
close all;
clear
lastwarn('');
warning off;
rng('default');
setup_path;

%% solver
useRiNNAL = 1;
useSDPNAL = 1; %%set to 1 if wants to run SDPNALplus

%% pars
probtype = 'QMSTP';
MBQPtype = 1; % 0. no strengen 1. (<=)
tol = 1e-6;
record = 1;
fname = feval(strcat(['problems_' probtype]))';

%% record results
rrALM   = {"data","n","eq","binary","rank","pobj","dobj","pfeas","dfeas","comp","pdgap","iter","itersub","time","iq"...
    "time_PG","PGiter","SSNiter","SSNCGiter","aveCG"};
rrNAL   = {"data","n","eq","binary","rank","pobj","dobj","pfeas","dfeas","comp","pdgap","iter","itersub","iterA","time","iq"};

%% start test
for i = [421:430 261:270 21:30 181:190 341:350]; % demo: [421]; 
    %   1:160 : esym (medium rank) 21:30   (m=1225) 101:110 (m=435)
    % 161:320 : sym  (high   rank) 181:190 (m=1225) 261:270 (m=435)
    % 321:480 : vsym (low    rank) 341:350 (m=1225) 421:430 (m=435)
    %%
    [Qtmp,ntmp,mtmp,edges] = load_data(fname{i},i);
    [Q,c,A,b,B,d,bidx,zid,n,l,m1,m2,neq,niq] = formulate_QMSTP(Qtmp,ntmp,mtmp,edges,MBQPtype);
    fprintf('\n ===> Start testing %s problems: %s \n',probtype,fname{i});

    %% RiNNAL+
    if useRiNNAL
        %%
        clear par
        par.tol = tol;
        r = [];
        %%
        par.projtype = 0;
        if contains(fname{i}, 'vsym')
            par.beta = 1e-5;
        else
            par.beta = 1e1;
        end
        %%
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
        [blk,At,Bt,C,bb,dd] = MBQP_to_SDPNAL(Q,c,A,b,B,d,bidx,zid);
        LL = 0;
        OPTIONS.tol = tol;
        OPTIONS.maxtime = 36000;
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



function [Q,c,A,b,B,d,bidx,zid,n,l,m1,m2,neq,niq] = formulate_QMSTP(Qtmp,ntmp,mtmp,edges,type)

n = size(Qtmp,1);
en = ones(n,1);
eyen = eye(n);
zn = zeros(n,n);

Q = Qtmp;
c = zeros(n,1);
A = [en'];
b = [ntmp-1];
bidx = 1:n;
zid = [];

if type == 1
    B = zeros(ntmp,mtmp);
    for i = 1:ntmp
        col_indices = find((edges(:,1) == i | edges(:,2) == i) & (edges(:,1) ~= edges(:,2)));
        B(i,col_indices) = -ones(1,length(col_indices));
    end
    d = -[ones(ntmp,1)];
elseif type == 0
    B = [];
    d = [];
else
    error('\n MBQP type error!');
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

m1 = size(A,1);
m2 = size(B,1);
n = size(Q,1);
l = length(bidx);
neq = m1+m1*n+l+1;
niq = m2*(2*n+m2+1)/2;

end

function [Q,n,m,edges] = load_data(rawname,i)

if i<= 600
    txt_name = strcat(rawname,'.txt');
    fileID = fopen(txt_name, 'r');

    % Check if the file opened successfully
    if fileID == -1
        error('Cannot open the file: %s', txt_name);
    end

    data = fscanf(fileID, '%d');
    n = data(2);
    m = data(1);
    Q = reshape(data(3:end),[m,m]);
    Q = (Q+Q')/2;
    edges = zeros(m,2);
    k = 1;
    for i = 1:n
        for j = i+1:n
            edges(k,:) = [i,j];
            k = k+1;
        end
    end
    fclose(fileID);
else
    txt_name = strcat(rawname,'.txt');
    [n, m, edges, Q] = readInput_aqmstp4qmstp(txt_name);
end

end

function [n, m, edges, Q] = readInput_aqmstp4qmstp(filepath)
% Open file
fid = fopen(filepath, 'r');

% Read n and m values
line = fgetl(fid);
n = sscanf(line, 'param n := %d ;');
line = fgetl(fid);
m = sscanf(line, 'param m := %d ;');

% Initialize variables
edges = [];
pos_edges = containers.Map('KeyType', 'char', 'ValueType', 'int32');

% Read edges
line = fgetl(fid);
edges_str = regexp(line, '\((\d+),(\d+)\)', 'tokens');
for i = 1:length(edges_str)
    edge = [str2double(edges_str{i}{1}), str2double(edges_str{i}{2})];
    edges = [edges; edge];
    pos_edges(mat2str(edge)) = i;
end

% Initialize cost matrix
Q = zeros(m, m);

% Read costs on edges
line = fgetl(fid);
cost_str = regexp(line, '\[(\d+),(\d+)\] (\d+)', 'tokens');
for i = 1:length(cost_str)
    edge = [str2double(cost_str{i}{1}), str2double(cost_str{i}{2})];
    cost = str2double(cost_str{i}{3});
    pos = pos_edges(mat2str(edge));
    Q(pos,pos) = cost;
end

% Read quadratic costs on pairs of adjacent edges
line = fgetl(fid);
line = fgetl(fid);
quad_cost_str = regexp(line, '\[(\d+),(\d+),(\d+),(\d+)\] (\d+)', 'tokens');
for i = 1:length(quad_cost_str)
    edge1 = [str2double(quad_cost_str{i}{1}), str2double(quad_cost_str{i}{2})];
    edge2 = [str2double(quad_cost_str{i}{3}), str2double(quad_cost_str{i}{4})];
    cost = str2double(quad_cost_str{i}{5});
    pos1 = pos_edges(mat2str(edge1));
    pos2 = pos_edges(mat2str(edge2));
    Q(pos1, pos2) = cost;
end

% Close file
fclose(fid);

% Ensure Q is symmetric
Q = (Q + Q') / 2;
end
