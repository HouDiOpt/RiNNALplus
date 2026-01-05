function [blk,At,Bt,C,bb,dd,flag]=MBQP_to_SDPNAL_test_strength(Q,c,A,b,B,d,bidx,zid)
%% pars
flag = 0;
[m,~] = size(A);
[p,~] = size(B);
n = size(Q,1);
neq = m+length(bidx)+1+length(zid)/2;
niq = p+n;
%% vars
blk{1,1} = 's'; blk{1,2} = n+1;
Acell = cell(1,neq);
Bcell = cell(1,niq);
%% Ax=b
for i = 1:m
    Acell{i} = [0,A(i,:)/2;A(i,:)'/2,sparse(n,n)];
end
%% diag_B(X)=x_B
for i = 1:length(bidx)
    ei = sparse(bidx(i),1,1,n,1); 
    eii = sparse(bidx(i),bidx(i),1,n,n);
    Acell{m+i} = sparse([0,ei'/2;ei/2,-eii]);
end
%% Y_{11} = 1
Acell{m+length(bidx)+1} = [1,sparse(1,n);sparse(n,1),sparse(n,n)];
%% X_{ij}=0, ij in zid
id0 = m+length(bidx)+1;
M = zeros(n+1,n+1);
M(zid) = 1;
M = triu(M);
[rowid,colid] = find(M);
for i = 1:length(zid)/2
    Acell{id0+i} = sparse([rowid(i),colid(i)],[colid(i),rowid(i)],[1/2,1/2],n+1,n+1);
end
%% d-Bx >= 0
for i = 1:p
    Bmat = -[0,B(i,:)/2;B(i,:)'/2,zeros(n,n)];
    Bcell{i} = sparse(Bmat);
end
%% x >= 0
for i = 1:n
    Bcell{p+i} = sparse([1,i+1],[i+1,1],[1 1],n+1,n+1);
end
%%
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
bb(m+length(bidx)+1) = 1;
dd = [-d;zeros(n,1)];
C = {[0,c';c,Q]};
end