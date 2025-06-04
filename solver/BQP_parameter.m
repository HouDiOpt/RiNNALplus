function parPR = BQP_parameter(A,b,B,d,bidx,par)
%%
[m,n] = size(A);
[m2,~] = size(B);
AAt = A*A';
At = A';
Atid = At(bidx,:);
AAtinv = full(inv((A*At)));
diagAAtinv = 1./diag(AAt);
Ainv = A'/(A*A');%pinv(A);
Ainvid = Ainv(bidx,:);
Aid = A(:,bidx);
projA = eye(n)-Ainv*A; projA = (projA+projA')/2;
projAid = projA(bidx,bidx);
%% manifold type
if isempty(A)&&isempty(bidx)
    type = 3; % no constraints
elseif isempty(bidx)
    type = 2; % only affine
elseif isempty(A)
    type = 1; % only diag, but B=[n]
else
    type = 0; % both constraints
end
%% binary type
if length(bidx) == n
    id_full = 1;
else
    id_full = 0;
end
%% proj type
if isfield(par,'projtype')
    projtype = par.projtype;
    parPR.projtype = projtype;
end
%% Amap
if isfield(par,'Amap1')
    parPR.Amaptype = 1;
    parPR.Amap1 = par.Amap1;
else
    parPR.Amaptype = 0;
end
%% Atmap
if isfield(par,'Atmap1')
    parPR.Atmaptype = 1;
    parPR.Atmap1 = par.Atmap1;
else
    parPR.Atmaptype = 0;
end
%% Btype
if isempty(B)
    Btype = 0;
elseif size(B,1) == size(B,2)
    Btype = 2;
    if norm(B-eye(size(B,1)),'fro') == 0
        Btype = 1;
        par.Bmap1  = @(x) x;
        par.Btmap1 = @(x) x;
    end
else
    Btype = 2;
end
parPR.Btype = Btype;
%% Bmap
if isfield(par,'Bmap1')
    parPR.Bmaptype = 1;
    parPR.Bmap1 = par.Bmap1;
else
    parPR.Bmaptype = 0;
end
%% Btmap
if isfield(par,'Btmap1')
    parPR.Btmaptype = 1;
    parPR.Btmap1 = par.Btmap1;
else
    parPR.Btmaptype = 0;
end
%% dimension
parPR.m = m;
parPR.m2 = m2;
parPR.n = n;
parPR.l = length(bidx);
%% assign
parPR.b  = b;
parPR.A  = A;
parPR.At = A';
parPR.d  = d;
parPR.B  = B;
parPR.Bt = B';
parPR.Atid = Atid;
parPR.Ae = sum(A,2);
parPR.AAt = AAt;
parPR.Ainv = Ainv;%pinv(A)
parPR.Ainvid = Ainvid;
parPR.Aid = Aid;
parPR.Ainvt = Ainv';
parPR.AAtinv = AAtinv;
parPR.projA = projA;
parPR.projAid = projAid;
parPR.diagAAtinv = diagAAtinv;
parPR.bidx = bidx;
parPR.id_full = id_full;
parPR.bidxdiff = setdiff(1:n,bidx);
parPR.type = type;
parPR.tol = par.tol;
parPR.tune_rank = 0;
if isfield(par,'inc_step0'); parPR.inc_step0 = par.inc_step0; end
if isfield(par,'tolrrate'); parPR.tolrrate = par.tolrrate; end
if isfield(par,'incmax'); parPR.incmax = par.incmax; end
if isfield(par,'redratio'); parPR.redratio = par.redratio; end
if isfield(par,'rankmax'); parPR.rankmax = par.rankmax; end
if isfield(par,'tune_rank'); parPR.tune_rank = par.tune_rank; end

%% ============================== PG ==============================
%% pars
P = [b';-A'];
bidx1 = [0 bidx] + 1;
PtPinv = inv(full(P'*P)); PtPinv = (PtPinv+PtPinv')/2;
PtPinvPt = PtPinv*P';
%% compute I-P(P'P)^{-1}P'
[Q, ~] = qr(P, 0);
if any(P)
    J = eye(n+1) - Q * Q';
    J = (J+J')/2;
else
    J = speye(n+1);
end
%% assign
parPR.P_PG = P;
parPR.bidx1 = bidx1;
diagAAt = [1;3/2*ones(length(bidx),1)];
parPR.invdiagAAt = 1./diagAAt;
parPR.b_PG  = ones(length(bidx)+1,1);
%% operators
parPR.Jmap_PG  = @(X) Jmap_PG(X,P,PtPinvPt,J);
parPR.Amap_PG  = @(X) Amap_PG(X,Q);
parPR.JJmap_PG = @(X) JJmap_PG(X,Q);
parPR.ATmap_PG = @(X) ATmap_PG(X,n,Q);
K = sparse([1,ones(1,n)/2;zeros(n,1),speye(n)/2]);
parPR.KKtmap = @(X) AAtmap(X,K);
parPR.KtKmap = @(X) AAtmap(X,K');
parPR.KinvtKinvmap = @(X) AAtmap(X,inv(full(K))');
parPR.KinvKinvtmap = @(X) AAtmap(X,inv(full(K)));
N = K*P;
parPR.N = N;
parPR.NN = N*N';
%% estimate subproblem L
if any(B)
    B1 = full([d -B]);
    B1norm2 = norm(B1,2)^2;
    parPR.L = norm(full(K),2)^4*(1+B1norm2*(B1norm2+1));
else
    parPR.L = norm(full(K),2)^2;
end
%% perturbation
parPR.p = ones(length(bidx),1);
end


%% ===============================================================

%% Y = JX
function Y = Jmap_PG(X,P,Pinv,J)
[n,r] = size(P);
if r <= n/2
    Y = X-P*(Pinv*X);
else
    Y = J*X;
end
end

%% y = A(X)
function y = Amap_PG(X,Q)
JXJ = JJmap_PG(X,Q);
y = diag(JXJ);
end

%% X = A'(y)
function X = ATmap_PG(y,n,Q)
X = spdiags(y,0,n+1,n+1);
X = JJmap_PG(X,Q);
X = (X+X')/2;
end

%% Y = JXJ = (I-QQ')X(I-QQ')
function Y = JJmap_PG(X,Q)
%% v1
% Y = J*X*J;
% Y = (Y+Y')/2;
%% v2 (J = I-QQ')
Y = X - (X*Q)*Q';
Y = Y - Q*(Q'*Y);
Y = (Y+Y')/2;
end

%% Y = KXK'
function Y = AAtmap(X,K)
Y = K*X*K';
Y = (Y + Y')/2;
end
