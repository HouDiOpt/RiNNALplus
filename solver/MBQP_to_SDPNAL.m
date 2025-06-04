%%*********************************************************************
%% generate problem data for SDPNALplus
%%*********************************************************************
function [blk,At,Bt,C,bb,dd,flag]=MBQP_to_SDPNAL(Q,c,A,b,B,d,bidx,zid)
%% pars
flag = 0;
[m,~] = size(A);
[p,~] = size(B);
n = size(Q,1);
neq = m+m*n+length(bidx)+1+length(zid)/2;
niq = p*n+p*(p+1)/2;
%% vars
blk{1,1} = 's'; blk{1,2} = n+1;
Acell = cell(1,neq);
Bcell = cell(1,niq);
%% Ax=b
for i = 1:m
    Acell{i} = [0,A(i,:)/2;A(i,:)'/2,sparse(n,n)];
end
%% AX=bx'
for i = 1:m
    for j = 1:n
        bi = sparse(j,1,-b(i)/2,n,1);
        AX = sparse(j,1:n,A(i,:),n,n);
        AX = (AX+AX')/2;
        Acell{m+(i-1)*n+j} = [0,bi';bi,AX];
    end
end
%% diag_B(X)=x_B
for i = 1:length(bidx)
    ei = sparse(bidx(i),1,1,n,1); 
    eii = sparse(bidx(i),bidx(i),1,n,n);
    Acell{m+m*n+i} = sparse([0,ei'/2;ei/2,-eii]);
end
%% Y_{11} = 1
Acell{m+m*n+length(bidx)+1} = [1,sparse(1,n);sparse(n,1),sparse(n,n)];
%% X_{ij}=0, ij in zid
id0 = m+m*n+length(bidx)+1;
M = zeros(n+1,n+1);
M(zid) = 1;
M = triu(M);
[rowid,colid] = find(M);
for i = 1:length(zid)/2
    Acell{id0+i} = sparse([rowid(i),colid(i)],[colid(i),rowid(i)],[1/2,1/2],n+1,n+1);
end
%% [(d-Bx)x']_L>=0
for j = 1:n
    for i = 1:p
        rowid  = [j+1,1,2:(n+1),(j+1)*ones(1,n)];
        colid  = [1,j+1,(j+1)*ones(1,n),2:(n+1)];
        elevec = [d(i),d(i),-B(i,:),-B(i,:)];
        Bcell{(j-1)*p+i} = sparse(rowid,colid,elevec,n+1,n+1);
    end
end
%% [(d-Bx)(d-Bx)']_L>=0
idp = 1;
for j = 1:p
    for i = 1:j
        c1 = -(d(i)*B(j,:)'+d(j)*B(i,:)')/2;
        c2 = (B(i,:)'*B(j,:)+B(j,:)'*B(i,:))/2;
        Bcell{p*n+idp} = sparse([d(i)*d(j),c1';c1,c2]);
        idp = idp+1;
    end
end
if neq+niq<=1000*10000
    fprintf('\n start svec of At!');
    At = svec(blk,Acell,1);
    fprintf('\n start svec of Bt!');
    Bt = svec(blk,Bcell,1);
else
    At = [];
    error('svec too slow! stop!');
    flag = 1;
end
bb = zeros(neq,1);
bb(1:m,1) = b;
bb(m+m*n+length(bidx)+1) = 1;
dd = zeros(niq,1);
C = {[0,c';c,Q]};
end