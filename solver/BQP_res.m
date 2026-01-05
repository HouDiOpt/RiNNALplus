function [pfeas,dfeas,comp,pdgap,pfval,dfval,Veig,deig,alphaA,lam1,lam2] = BQP_res(X,Wslack,W,Uslack,U,Vslack,V,lam,Q,c,par)
incmax = 1;
if isfield(par,'incmax'); incmax = par.incmax; end
l = par.l;
m = par.m;
A = par.A;
B = par.B;
b = par.b;
d = par.d;
id = par.bidx;
Ainv = par.Ainv;
R0 = X.R; 
[n,r] = size(R0);
R = R0;
R(:,1) = R(:,1)+1;
R = R/2;
RR = R*R';
RR1 = [1,R(:,1)';R(:,1),RR];
RR1n = norm(RR1,"fro");
e1t = zeros(1,r); e1t(1) = 1;
R1 = [e1t;R];
dRR = diag(RR)-R(:,1);
ARR = (A*R)*R';
QR = Q*R;
Btype = par.Btype;
if Btype == 0
    de1tmBRRt = 0;
    de1tmBRsqr = 0;
elseif Btype == 1
    ee1tmR = -R; % (ee_1'-R)
    ee1tmR(:,1) = ee1tmR(:,1)+1;
    de1tmBRRt = R(:,1)'-RR; % (ee_1'-R)R'
    de1tmBRsqr = ee1tmR(:,1)-de1tmBRRt;
elseif Btype == 2
    de1tmBR = -B*R;
    de1tmBR(:,1) = de1tmBR(:,1)+d; % (be_1'-AR)
    de1tmBRRt = de1tmBR*R';        % (be_1'-AR)R'
    de1tmBRsqr = de1tmBR*de1tmBR'; % (be_1'-AR)*(be_1'-AR)'
end

%% recover subproblem dual variable: lambda_hat
muhat  = lam(1:l);
lamhat = lam(l+1:end);
lammat = reshape(lamhat,m,r); 
%% recover dual variables
JA = eye(n)-Ainv*A; JA = (JA+JA')/2;
muB = zeros(n,1); muB(id) = muhat;
L = Q-diag(muB)-W(2:end,2:end)-B'*V*B+(B'*U+U'*B)/2;
if any(B)
    h = 2*c+muB-2*W(2:end,1)+2*B'*V*d-U'*d;
else
    h = 2*c+muB-2*W(2:end,1);
end
S = [0,-zeros(1,n);-zeros(n,1),L];
S21 = h/2 + L*(Ainv*b);
alphaA = b'*lammat(:,1)/2 + h'*R(:,1)/2 - b'*Ainv'*h/2 - (h'/2 + b'*Ainv'*L)*(Ainv*b);
S = S + [-alphaA,S21';S21,zeros(n,n)];
S = [1,zeros(1,n);zeros(n,1),JA]*S*[1,zeros(1,n);zeros(n,1),JA];
S = (S+S')/2;
% S = [-R(:,1),eye(n)]'*JA*L*JA*[-R(:,1),eye(n)]; % not good when KKT approximately hold 
nS = norm(S,'fro');
%% escaping S direction (eigs)
if par.tune_rank == 1
    if n >= 3000
        opts.issym = true;
        [Veig,deig,flag] = eigs(S,incmax,'smallestreal',opts);
    else
        flag = 1;
    end
    if flag == 0
        deig = diag(deig);
        [deig,ind] = sort(deig);
        ind = ind(~isnan(deig));
        deig = deig(~isnan(deig));
        Veig = Veig(:,ind);
        Veig = JA*[-R(:,1),eye(n)]*Veig;
        deig = deig.*(deig<0);
    else
        [Veig,deig] = eig(S,'vector');
        [deig,ind] = sort(deig);
        Veig = Veig(:,ind);
        Veig = JA*[-R(:,1),speye(n)]*Veig;
        deig = deig.*(deig<0);
    end
else
    if n >= 3000
        opts.issym = true;
        [~,deig,flag] = eigs(S,1,'smallestreal',opts);
    else
        flag = 1;
    end
    if flag == 0
        deig = diag(deig);
        [deig,~] = sort(deig);
        deig = deig(~isnan(deig));
        Veig = [];
        deig = deig.*(deig<0);
    else
        deig = eig(S);
        [deig,~] = sort(deig);
        Veig = [];
        deig = deig.*(deig<0);
    end
end
%%
alpha = -S(1,1);
lam1 = Ainv'*(h+L*(Ainv*b));
% lam2 = Ainv'*L*(2*eye(n)-Ainv*A);
lam2 = [];
%% compute KKT
pfval = Produc(QR,R)+2*c'*R(:,1);
dfval = lam1'*b + alpha - Produc(V,d*d');
pfeas = sqrt( norm(A*R(:,1)-b,'fro')^2 ...
            + norm(ARR-b*R(:,1)','fro')^2 ...
            + norm(dRR(id))^2 ...
            + norm(RR1-Wslack,'fro')^2 ...
            + norm(de1tmBRRt-Uslack,'fro')^2 ...
            + norm(de1tmBRsqr-Vslack,'fro')^2 ...
            ) / ...
        (1 + sqrt(norm(b)^2 + norm(d*d','fro')^2 + 1));
%%
dfeas = norm(deig,'fro')/(1+nS);
comp = abs(Produc(RR1,S))/(1+RR1n+nS);
% no need for check (=0):
% abs(Produc(Wslack,W))/(1+nY+norm(W,'fro')))
pdgap = abs(pfval-dfval)/(1+abs(pfval)+abs(dfval));
end