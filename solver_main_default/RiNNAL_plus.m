function [fval,X,info] = RiNNAL_plus(Q,c,A,b,B,d,p,r,E,par)

%% The RiNNAL+ solver aims to solve the SDP-RLT relaxation of MBQP problem.
% Authors: Di Hou, Tianyun Tang, Kim-Chuan Toh.
% Update : 2025.04.10

%% Mixed-bianry Quadratic Programming (MBQP):
% min  x'Qx + 2c'x
% s.t. Ax = b                                (equality)
%      Bx <= d                             (inequality)
%      x_{i} in {0,1}, i in p                  (binary)
%      x_i*x_j = 0, (i+1,j+1) in E    (complementarity)
%      x >= 0                             (nonnegative)

%% The SDP-RLT relaxation of (MBQP):
%   min  <Q,X> + 2<c,x>               -- dual variables --
%   s.t. 1. Ax = b                    (recovered from lam)
%        2. AX = bx'                  (recovered from lam)
%        3. diag(X)_B = x_B           (lam)
%        4. (Bx <= d)                 (usually redundant)
%        5. BX <= dx'                 ( U )
%        6. BXB'-Bxd'-dx'B'+dd' >= 0  ( V )
%        7. Y_{ij} = 0, (i,j) in E    ( W )
%        8. Y >= 0                    ( W )
%        9. Y := [1,x';x,X] PSD       ( S )
%
% 1. RiNNAL+ uses ALM to penalize constraints (5-8),
%    while preserving constraints (1-3,9) in the ALM subproblem.
% 2. The low-rank decomposition is applied to Y:
%    Y = [e1, (R'+e1*e')/2]'[e1, (R'+e1*e')/2].

%% Input:
% Q    - (n  x n) Symmetric matrix defining the quadratic objective term
% c    - (n  x 1) Vector defining the linear objective term
% A    - (m1 x n) Equality constraint coefficient matrix
% b    - (m1 x 1) Right-hand side vector for equality constraints
% B    - (m2 x n) Inequality constraint coefficient matrix
% d    - (m2 x 1) Right-hand side vector for inequality constraints
% p    - (1  x l) Indices for binary constraints
% r    - Scalar specifying the rank of the approximation (auto-set if empty)
% E    - Index vector containing index for complementarity constraints 
% par  - Structure containing parameters (optional fields listed below):
%        . tol       : Convergence tolerance (default = 1e-6)
%        . verbose   : Display output level (default = 1)
%        . beta      : Initial penalty parameter for ALM (default = 1)
%        . Amap1     : Operator to compute A(x)  (otherwise use A*x)
%        . Atmap1    : Operator to compute A'(y) (otherwise use A'*y)
%        . maxtime   : Maximum runtime in seconds (default = 3600)
%        . tune_rank : Enable adaptive rank adjustment (default = 0)
%        . projtype  : Projection subproblem solver
%          [ 0: schur_chol, 1: schur_pcg, 2: chol, 3: pcg (default = 3) ]
%        . SSN_frequence : Frequency of SSN updates in PG step (default = 5)

%% Output:
% fval - A 2-element vector containing:
%        fval(1): the final primal objective value (x'Qx + 2c'x)
%        fval(2): the final dual objective value from the SDP relaxation
%
% X    - A structure containing the low-rank factor R of the solution matrix Y,
%        and the reconstructed matrix RR1 = [1, x'; x, X] (i.e., X.RR1 = Y ≽ 0)
%
% info - A structure with diagnostic and performance information, including:
%        .pfeas        : final primal feasibility violation
%        .dfeas        : final dual feasibility violation
%        .comp         : final complementarity gap
%        .pdgap        : relative duality gap between fval(1) and fval(2)
%        .ttime        : total solver runtime (in seconds)
%        .ALMite       : total number of ALM outer iterations
%        .BBite        : total number of Riemannian BB iterations
%        .lambda       : final Lagrange multiplier W for constraint (9)
%        .lam          : dual multiplier for constraints (1)-(3)
%        .rankrecord   : history of rank values over ALM iterations
%        .ttime_PG     : total time spent on PG steps
%        .PGiterTotal  : total number of PG iterations
%        .SSNTotal     : total number of semismooth Newton iterations
%        .CGiterTotal  : total number of CG steps in SSN
%        .aveCGiter    : average CG steps per SSN iteration

%% Pars for ALM
tstart = clock;
rng('default');
verbose = 1;
vergap = 1;
tol = 1e-6;
maxiter = 1e5;
beta = 1;
ra = 1e-16;
scale = 1;
ramin = 1e-16;
ramax = 1e-16;
betamin = 1e-16;
betamax = 1e16;
subscale = 1;
uptol = 1e-16;
maxtime = 3600;
redfrequence = 1;
rankrecord = [];
goal_subratio = 0.9;
useSSN      = 1;
SSN_frequence = 5;
tune_rank = 0;
rankmin = 2;
adaptmaxiter = 0;
pdratio1 = 5;
pdratio2 = 0.5;
accuratePGsub = 0;
if isfield(par,'verbose'); verbose = par.verbose; end
if isfield(par,'vergap'); vergap = par.vergap; end
if isfield(par,'tol'); tol = par.tol; end
if isfield(par,'beta'); beta = par.beta; end
if isfield(par,'betamin'); betamin = par.betamin; end
if isfield(par,'betamax'); betamax = par.betamax; end
if isfield(par,'ra'); ra = par.ra; end
if isfield(par,'subscale'); subscale = par.subscale; end
if isfield(par,'ramin'); ramin = par.ramin; end
if isfield(par,'ramax'); ramax = par.ramax; end
if isfield(par,'scale'); scale = par.scale; end
if isfield(par,'uptol'); uptol = par.uptol; end
if isfield(par,'maxtime'); maxtime = par.maxtime; end
if isfield(par,'redfrequence'); redfrequence = par.redfrequence; end
if isfield(par,'goal_subratio'); goal_subratio = par.goal_subratio; end
if isfield(par,'SSN_frequence'); SSN_frequence = par.SSN_frequence; end
if isfield(par,'tune_rank'); tune_rank = par.tune_rank; end
if isfield(par,'useSSN'); useSSN = par.useSSN; end
if isfield(par,'rankmin'); rankmin = par.rankmin; end
if isfield(par,'adaptmaxiter'); adaptmaxiter = par.adaptmaxiter; end
if isfield(par,'pdratio1'); pdratio1 = par.pdratio1; end
if isfield(par,'pdratio2'); pdratio2 = par.pdratio2; end
if isfield(par,'accuratePGsub'); accuratePGsub = par.accuratePGsub; end

%% check ra range
if ramin>ramax
    if isfield(par,'rankmin')
        ramax = ramin;
    else
        ramin = ramax;
    end
end

%% Scale the problem (Q,c)
if (scale==1)
    nCorg = sqrt(norm(Q,'fro')^2+norm(c,'fro')^2);
    nC = max(1,nCorg);
elseif (scale==0)
    nC = 1;
end
Corg = [0,c';c,Q];
Qs = Q/nC;
cs = c/nC;

%% Check format of input data
n = size(Q,1);
if size(B,1) == 0; B = zeros(0,n); d = zeros(0,1); end
if size(A,1) == 0; A = zeros(0,n); b = zeros(0,1); end
if size(c,1) == 0; c = zeros(n,1); end
m1 = size(A,1);
m2 = size(B,1);
l = length(p);
neq = m1+m1*n+l+1;
niq = m2*(2*n+m2+1)/2;

%% Print header
if verbose
    fprintf('\n ---------------------------------------------------------------------------------------------------');
    fprintf('\n RiNNAL+: a Riemannian ALM for SDP-RPT / DNN relaxations')
    fprintf('\n ---------------------------------------------------------------------------------------------------');
    fprintf('\n MBQP (dim,binary,equality,inequality) = (%d,%d,%d,%d)',n,l,m1,m2);
    fprintf('\n SDP  (dim,       equality,inequality) = (%d,%d,%d)',n,neq,niq);
    fprintf('\n ---------------------------------------------------------------------------------------------------');
    fprintf('\n ite |     fval     |  beta  |  pfeas    dfeas    comp     pdgap |   ra   |gradnorm|ites| time |rank|');
end

%% Construct operator
parPR = BQP_parameter(A,b,B,d,p,par);
Rie_proj = @(C,X,parPR)BQP_proj(C,X,parPR);
Rie_retrac = @(X,H,s,tol,parPR)BQP_retrac_perturb(X,H,s,tol,parPR);
Bmap0  = @(X) Bmap(X,parPR,0);
Btmap0 = @(X) Btmap(X,parPR,0);

%% Initialize R
if isempty(r)
    r = min(200,max(5,round(n/5)));
end
parPR.checkres = 1;
if isfield(par,'x0')
    R0 = 1e-5*randn(n,r);
    R0(:,1) = 2*par.x0-1;
    X.R = R0;
    X = Rie_retrac(X,0,0,1e-5,parPR);
else
    X.R = randn(n,r);
    X = Rie_retrac(X,0,0,1e-5,parPR);
end

%% Pars for BB
parBB.tolg = 1e2;
parBB.verbose = 0;
parBB.safegard = 1;
parBB.eigratio = 1e2;
parBB.maxiter = 50;
parBB.alp = 1;
parBB.adstop = 1;
parBB.checkrank = 0;
parBB.connec = 1;
parBB.plotyes = 0;
parBB.scaler = nC;
parBB.uptol  = uptol;
parBB.comp_pfeas_fre = 1e4;
if isfield(par,'safegard'); parBB.safegard = par.safegard; end
if isfield(par,'adstop');parBB.adstop = par.adstop; end
if isfield(par,'BBmaxiter');parBB.maxiter = par.BBmaxiter; end
if isfield(par,'comp_pfeas_fre');parBB.comp_pfeas_fre = par.comp_pfeas_fre; end
if isfield(par,'BBverbose');parBB.verbose = par.BBverbose; end

%% Initialization
W = zeros(n+1,n+1);
U = zeros(m2,n);
V = zeros(m2,m2);
pfeas   = Inf;
subfeas = Inf;
res     = Inf;
BBiterTotal = 0;
PGiterTotal = 0;
SSNTotal    = 0;
CGiterTotal = 0;
retrac_flag = 0;
useRie      = 1;
ttime_PG    = 0;
retrac_singular_count_sum = 0;
parPR.checkres = 0;
singular_flag = 0;
maxiter0 = parBB.maxiter;
singular_all = [];
pfeas_all = [];
dfeas_all = [];
comp_all  = [];
subfeas_all = [];
res_all = [];

%% Main loop
for i = 1:maxiter

    %% pre
    W0 = W; U0 = U; V0 = V;
    pfeas0 = pfeas;
    parBB.ALMite = i;
    rankrecord = [rankrecord, size(X.R,2)];

    %% =========================== Rie GD =============================
    if useRie

        %% construct low-rank augmented Lagrangian function
        f_g = @(X,flag)f_and_g_Rie(X,Qs,cs,W0,U0,V0,d,beta,E,flag,Bmap0,Btmap0,parPR.Btype);

        %% remove perturbation
        parPR.p = ones(length(p),1);
        Rie_retrac = @(X,H,s,tol,parPR)BQP_retrac_perturb(X,H,s,tol,parPR);

        %% Rie GD
        [X,lam,iter,~,~,gradnorm,retrac_flag,retrac_singular_count] = ...
           Rie_BB_Sub_perturb(f_g,Rie_proj,Rie_retrac,X,ra,parBB,parPR);
        BBiterTotal = BBiterTotal + iter;
        retrac_singular_count_sum = retrac_singular_count_sum + retrac_singular_count;

        %% add perturbation
        perturb_scale = 1;
        while retrac_flag
            pert = rand(length(p),1);
            pconst = min(1e-3,max(tol));
            pconst = pconst*perturb_scale*sqrt(l)*pert/norm(pert);
            fprintf('\n perturbation_const: %3.1e',norm(pconst));
            parPR.p = (1+pconst).*ones(length(p),1);
            %%
            Rie_retrac = @(X,H,s,tol,parPR)BQP_retrac_perturb(X,H,s,tol,parPR);
            [X,lam,iter,~,~,gradnorm,retrac_flag,retrac_singular_count] = ...
               Rie_BB_Sub_perturb(f_g,Rie_proj,Rie_retrac,X,ra,parBB,parPR);
            %%
            if retrac_flag == 0
                fprintf('\n perturbation work!');
            elseif perturb_scale >= 1e5
                fprintf('\n perturbation too large, cannot work!');
                break
            end
            perturb_scale = perturb_scale*10;
        end

        %% auxilary variables
        R = X.R;
        R(:,1) = R(:,1)+1;
        R = R/2;
        RR = R*R';
        RR1 = [1,R(:,1)';R(:,1),RR];
        de1tmBR = -B*R;
        de1tmBR(:,1) = de1tmBR(:,1)+d; % (de_1'-BR)
        de1tmBRRt = de1tmBR*R';        % (de_1'-BR)R'
        de1tmBRsqr = de1tmBR*de1tmBR'; % (de_1'-BR)*(de_1'-BR)'

        %% update  dual  variable
        Wtmp = RR1-W0/beta;
        Wslack = Projpoly(Wtmp,E);
        W = (Wslack-Wtmp)*beta;

        Utmp = U0-beta*(de1tmBRRt);
        U = max(0,Utmp);
        Uslack = (U-Utmp)/beta;

        Vtmp = V0-beta*(de1tmBRsqr);
        V = max(0,Vtmp);
        Vslack = (V-Vtmp)/beta;

        %% compute residue
        [pfeas,dfeas,comp,pdgap,pfval,dfval,Veig,deig,lam11,lam1,lam2] = ...
            BQP_res(X,Wslack,W*nC,Uslack,U*nC,Vslack,V*nC,2*lam*nC,Q,c,parPR);
        subfeas  = subscale*max(dfeas,comp);
        pdratio  = subfeas/pfeas;
        res = max([pfeas,dfeas,comp]);
        ttime = etime(clock,tstart);
        final_rank = size(X.R,2);
        if verbose && mod(i,vergap) == 0
            fprintf('\n %3d |%+6.7e|%3.2e|%3.2e %3.2e %3.2e %3.2e|%3.2e|%3.2e|%4d|%6.1f|',...
                i,pfval,beta,pfeas,dfeas,comp,pdgap,ra,gradnorm,iter,ttime);
            fprintf('%4d|', final_rank);
        end
    end
    
    %% check stopping criterion
    if res < tol
        fprintf('\n ---------------------------------------------------------------------------------------------------');
        fprintf('\n Residue = %3.2e < %3.2e\n SDP is solved to the required accuracy',res,tol);
        fprintf('\n ---------------------------------------------------------------------------------------------------');
        info.retrac_flag = 0;
        break
    end
    if ttime > maxtime
        fprintf('\n\n Reach time limit, Stop! \n');
        break
    end

    %% tune rank
    if tune_rank == 1
        %% increase rank
        inc_flag = 0;
        if pdratio > goal_subratio
            [X.R,inc_flag,~,retrac_singular_count] = ...
            rankinc(X.R,Veig,deig/nC,f_g,Rie_proj,Rie_retrac,parPR);
        end
        if inc_flag == 0
            fprintf(' ');
        end
        retrac_singular_count_sum = retrac_singular_count_sum + retrac_singular_count;

        %% reduce rank
        if size(X.R,2)>1 && mod(i,redfrequence) == 0
            [X,~] = rankred(X,Rie_retrac,parPR);
        end

    end

    %% ============================= PG ===============================
    if useSSN && (mod(i-1,SSN_frequence) == 0 || singular_flag)
        tstart_PG = clock;

        %% record
        W1 = W;
        U1 = U;
        V1 = V;
        X1 = X;
        pfeas1 = pfeas;
        dfeas1 = dfeas;
        comp1  = comp;
        subfeas1 = subfeas;
        res1 = res;

        %% Parameter for PG
        parPG.verbose = 0;
        parPG.maxiter = 1;
        parPG.beta    = beta;
        parPG.pfeas   = pfeas;
        parPG.subfeas = subfeas;
        parPG.tol     = tol;
        parPG.lam     = lam;
        parPG.lam11   = lam11;
        parPG.lam1    = lam1;
        parPG.lam2    = lam2;
        parPG.nC      = nC;
        parPG.printPG = 0;
        parPG.accuratePGsub = accuratePGsub;
        if isfield(par,'PGmaxiter');   parPG.maxiter     = par.PGmaxiter; end
        if isfield(par,'PGtol');       parPG.tol         = par.PGtol; end
        if isfield(par,'printPG');     parPG.printPG     = par.printPG; end
        if isfield(par,'printSSN');    parPG.printSSN    = par.printSSN; end
        if isfield(par,'printSSNsub'); parPG.printSSNsub = par.printSSNsub; end

        %% retrac accurately (optional)
        X = Rie_retrac(X,0,0,1e-8,parPR);

        %% recover RR'
        R = X.R;
        R(:,1) = R(:,1)+1;
        R = R/2;
        X.RR1 = [1,R(:,1)';R(:,1),R*R'];
        Xp = X;

        %% construct f and g
        f_g = @(X,cp_flag,lr_flag)f_and_g_PG(X,Qs,cs,W0,U0,V0,B,d,beta,E,...
            cp_flag,Bmap0,parPR.Btype,parPR.KKtmap,parPR.KtKmap);

        %% record hist
        rp = size(X.R,2);

        %% PG step
        [X,~,~,y,~,info_PG] = PG(f_g,X,parPG,parPR);
        parPG.y = y;
        SSNTotal    = SSNTotal    + info_PG.iterSSN;
        PGiterTotal = PGiterTotal + info_PG.iter;
        CGiterTotal = CGiterTotal + info_PG.iterCG;

        %% compute residue of new X
        RR1 = X.RR1;
        R1  = RR1(2:end,1);
        RR  = RR1(2:end,2:end);
        de1tmBRRt = d*R1'-B*RR;
        de1tmBRsqr = d*d'-d*(R1'*B')-B*R1*d'+B*RR*B'; % (de_1'-BR)*(de_1'-BR)'

        Wtmp = RR1-W0/beta;
        Wslack = Projpoly(Wtmp,E);
        W = (Wslack-Wtmp)*beta;

        Utmp = U0-beta*(de1tmBRRt);
        U = max(0,Utmp);
        Uslack = (U-Utmp)/beta;

        Vtmp = V0-beta*(de1tmBRsqr);
        V = max(0,Vtmp);
        Vslack = (V-Vtmp)/beta;

        dRR = diag(RR)-R1;

        %% compute sub res
        pfeas = sqrt(norm(A*R1-b,'fro')^2 ...
            + norm(A*RR-b*R1','fro')^2 ...
            + norm(dRR(p))^2 ...
            + norm(RR1-Wslack,'fro')^2 ...
            + norm(de1tmBRRt-Uslack,'fro')^2 ...
            + norm(de1tmBRsqr-Vslack,'fro')^2 ...
            ) / ...
            (1 + sqrt(norm(b)^2 + norm(d*d','fro')^2 + 1));
        subpfeas = info_PG.subpfeas;
        subdfeas = info_PG.subdfeas;
        subcomp  = info_PG.subcomp;
        subfval  = full(Produc(Corg,X.RR1));
        pfeas    = max(pfeas,subpfeas);
        dfeas    = subdfeas;
        comp     = subcomp;
        subfeas  = subscale*max(dfeas,comp);
        res = max([pfeas,subpfeas subdfeas subcomp]);

        %% check res
        if res < tol
            final_rank = info_PG.rank;
            %% print
            if verbose && mod(i,vergap) == 0
                ttime_PG = ttime_PG + etime(clock,tstart_PG);
                fprintf('\n   * |%+6.7e|%3.2e|%3.2e %3.2e %3.2e %3.2e|%3.2e|%3.2e|%4d|%6.1f|',...
                    subfval,beta,pfeas,subdfeas,subcomp,pdgap,ra,gradnorm,...
                    info_PG.iterSSN,ttime+etime(clock,tstart_PG));
                fprintf('%4d|', final_rank);
            end
            %%
            fprintf('\n ---------------------------------------------------------------------------------------------------');
            fprintf('\n (PG) Residue = %3.2e < %3.2e\n SDP is solved to the required accuracy',res,tol);
            fprintf('\n ---------------------------------------------------------------------------------------------------');
            info.retrac_flag = 0;
            break
        end

        %% check singular
        if subfeas < 5e-2*subfeas1
            singular_all = [singular_all 1];
        else
            singular_all = [singular_all 0];
        end

        if i>=6 && subfeas1 >= 0.1*tol && subfeas < 1e-2*subfeas1
            singular_flag = 1;
        elseif length(singular_all)>=3
            if sum(singular_all(end-2:end))==0
                singular_flag = 0;
            end
        end

        %% check update or not
        if abs(info_PG.rank-size(Xp.R,2))>=1 || info_PG.decrease %(pfeas<1.1*pfeas1 && subfeas<0.9*subfeas1)

            %% update rank
            [VV,DD] = eig(X.RR1);
            dd = diag(DD);
            [dd,did] = sort(dd,'descend');
            VV = VV(:,did);

            %% plot
            % figure1 = figure();
            % semilogy(sort(eig(Xp.RR1),'descend'),'*');
            % hold on
            % semilogy(dd,'o');
            % legend('eig-old','eig-new')
            % delete(figure1)

            %% find rank
            l2 = info_PG.rank;
            if isempty(l2)
                rn = rp;
            else
                rn = l2(1);
            end
            rn = max(rankmin,rn);
            rn = min([rn,rp+round(n/5),n]);
            tmp1 = VV(:,1:rn)*spdiags(sqrt(dd(1:rn)),0,rn,rn);

            % Householder
            a1 = tmp1(1,:)';
            b1 = [1;zeros(rn-1,1)];
            a_norm = a1 / norm(a1);
            b_norm = b1 / norm(b1);
            u = a_norm - b_norm;
            u = u / norm(u);
            P = eye(length(a1)) - 2 * (u * u');
            %
            tmp2 = tmp1*P;
            tmp2 = tmp2*2;
            tmp2(:,1) = tmp2(:,1)-1;
            X.R  = tmp2(2:end,:);

            % keep original KKT
            if res >= res1
                pfeas = pfeas1;
                dfeas = dfeas1;
                comp  = comp1;
                subfeas = subfeas1;
            end

        else
            W = W1;
            U = U1;
            V = V1;
            X = X1;
            pfeas = pfeas1;
            dfeas = dfeas1;
            comp  = comp1;
            subfeas = subfeas1;
            res = res1;
        end

        ttime_PG = ttime_PG + etime(clock,tstart_PG);

        %% print
        if verbose && mod(i,vergap) == 0
            fprintf('\n   * |%+6.7e|%3.2e|%3.2e %3.2e %3.2e %3.2e|%3.2e|%3.2e|%4d|%6.1f|',...
                subfval,beta,pfeas,dfeas,comp,pdgap,ra,gradnorm,info_PG.iterSSN,ttime+etime(clock,tstart_PG));
            fprintf('%4d|', size(X.R,2));
        end
    end

    %% update beta and ra
    if (pdratio1*pfeas < subfeas || pfeas < 0.5*tol)
        beta = beta/1.5;
    elseif (pdratio2*pfeas > subfeas && pfeas > 0.7*pfeas0) || (0.01*pfeas > subfeas)
        beta = beta*1.5;
    end

    %% record KKT info
    pfeas_all = [pfeas_all pfeas];
    dfeas_all = [dfeas_all dfeas];
    comp_all  = [comp_all  comp ];
    subfeas_all = [subfeas_all subfeas];
    res_all = [res_all res];

    %% tune maxiter
    if adaptmaxiter && i>=50
        if pfeas_all(end)  > 0.9*min(pfeas_all(end-40:end-20)) && ...
                subfeas_all(end)>0.9*min(subfeas_all(end-40:end-20))
            parBB.maxiter = min(parBB.maxiter*2,maxiter0*16);
        end
    end
    parBB.maxiter = max(parBB.maxiter,1);
    parBB.maxiter = min(parBB.maxiter,10000);
    beta = max(betamin,min(betamax,beta));
    ra = max(ramin,min(ramax,ra));
end

%% Record results
fval(1) = pfval;
fval(2) = dfval;
X.RR1   = RR1;
info.pfeas = pfeas;
info.dfeas = dfeas;
info.comp  = comp;
info.pdgap = pdgap;
info.ttime = ttime;
info.ALMite = i;
info.BBite  = BBiterTotal;
info.lambda = W*nC;
info.lam = 2*lam*nC;
info.rankrecord = rankrecord;
info.ttime_PG = ttime_PG;
info.PGiterTotal = PGiterTotal;
info.SSNTotal = SSNTotal;
info.CGiterTotal = CGiterTotal;
info.aveCGiter = CGiterTotal/max(1,SSNTotal);

%% Print
fprintf('\n max dim. of SDP variable  = %d',n+1);
fprintf('\n num. of equality   constraints  = %d',neq);
fprintf('\n num. of inequality constraints  = %d',niq);
fprintf('\n Computing time (total)          = %.1f',ttime);
fprintf('\n Computing time (RieBB)          = %.1f',ttime-ttime_PG);
fprintf('\n Computing time (PG)             = %.1f',ttime_PG);
fprintf('\n number of RieBB iter            = %d',BBiterTotal);
fprintf('\n number of PG    iter            = %d',PGiterTotal);
fprintf('\n number of SSN   iter            = %d',SSNTotal);
fprintf('\n total number of cg step         = %d',CGiterTotal);
fprintf('\n average number of cg step per SSN subproblem = %.1f',CGiterTotal/max(1,SSNTotal));
fprintf('\n -----------------------------------------------------');
fprintf('\n primal objval  = %6.7e',fval(1));
fprintf('\n dual   objval  = %6.7e',fval(2));
fprintf('\n relative gap   = %3.2e',pdgap);
fprintf('\n primfeasorg    = %3.2e',pfeas);
fprintf('\n dualfeasorg    = %3.2e',dfeas);
fprintf('\n complement     = %3.2e',comp);
fprintf('\n final rank     = %d',final_rank);
fprintf('\n ---------------------------------------------------------------------------------------------------\n\n');
end

%%**********************************************************************
%% Functions
function [f,eugf,pfeas,RR1] = f_and_g_Rie(X,Q,c,W,U,V,d,beta,zid,comp_pfeas_flag,Bmap0,Btmap0,Btype)

R = X.R;
R(:,1) = R(:,1)+1;
R = R/2;
r = size(R,2);
e1t = zeros(1,r); e1t(1) = 1;
R1 = [e1t;R];
RR1 = R1*R1';
RRt = RR1(2:end,2:end);

if Btype == 0
    de1tmARRt  = 0;
    de1tmARsqr = 0;
elseif Btype == 1
    de1tmARRt  = R(:,1)'-RRt; % (be_1'-AR)R'
    de1tmARsqr = (d-R(:,1))-de1tmARRt; % (be_1'-AR)*(be_1'-AR)'
elseif Btype == 2
    de1tmARRt  = d*R(:,1)'-Bmap0(RRt); % (be_1'-AR)R'
    de1tmARsqr = (d-(Bmap0(R(:,1))))*d'-Bmap0(de1tmARRt')'; % (be_1'-AR)*(be_1'-AR)'
end
%%
Wtmp = RR1-W/beta;
Wslack = Projpoly(Wtmp,zid);
Wp = (Wslack-Wtmp)*beta;

Utmp = U-beta*(de1tmARRt);
Up = max(0,Utmp);

Vtmp = V-beta*(de1tmARsqr);
Vp = max(0,Vtmp);

%%
f = Produc(Q,RRt)+2*c'*R(:,1)+norm(Wp,'fro')^2/(2*beta) ...
    +norm(Up,'fro')^2/(2*beta)+norm(Vp,'fro')^2/(2*beta);
%%
if Btype == 0
    eugf = (Q-Wp(2:end,2:end))*R;
    eugf(:,1) = eugf(:,1)+c-Wp(2:end,1);
elseif Btype == 1
    eugf = (Q-Wp(2:end,2:end)-Vp+(Up+Up')/2)*R;
    eugf(:,1) = eugf(:,1)+c-Wp(2:end,1)+sum(Vp,2)-sum(Up,1)'/2;
elseif Btype == 2
    BtUp = Btmap0(Up);
    eugf = (Q-Wp(2:end,2:end)-Btmap0(Btmap0(Vp)')'+(BtUp+BtUp')/2)*R;
    eugf(:,1) = eugf(:,1)+c-Wp(2:end,1)+Btmap0((Vp*d))-Up'*d/2;
end

%%
if comp_pfeas_flag == 1
    Uslack = (Up-Utmp)/beta;
    Vslack = (Vp-Vtmp)/beta;
    pfeas = max([ ...
        norm(RR1-Wslack,'fro')/(1+norm(RR1,'fro')), ...
        norm(de1tmARRt-Uslack,'fro')/(1+norm(Uslack,'fro')), ...
        norm(de1tmARsqr-Vslack,'fro')/(1+norm(Vslack,'fro')), ...
        ]);
else
    pfeas = inf;
end

end

%%**********************************************************************
function [f,eugf,pfeas] = f_and_g_PG(X,Q,c,W,U,V,B,d,beta,zid,comp_pfeas_flag,Bmap0,Btype,KKtmap,KtKmap)

C   = [0,c';c,Q];
dB  = [d -B];
n   = size(B,2);
zI  = [sparse(1,n);speye(n)];
RR1 = X.RR1;
RR1 = KtKmap(RR1);
R1  = RR1(2:end,1);
RR  = RR1(2:end,2:end);

if Btype == 0
    de1tmBRRt  = 0;
    de1tmBRsqr = 0;
elseif Btype == 1
    de1tmBRRt  = R1'-RR;
    de1tmBRsqr = RR -R1-R1'+RR1(1,1);
    de1tmBRsqr = (de1tmBRsqr+de1tmBRsqr')/2;
elseif Btype == 2
    de1tmBRRt = d*R1'-Bmap0(RR); % (de_1'-BR)R'
    de1tmBRsqr = dB*RR1*dB'; % (de_1'-BR)*(de_1'-BR)'
    de1tmBRsqr = (de1tmBRsqr+de1tmBRsqr')/2;
end

Wtmp = RR1-W/beta;
Wslack = Projpoly(Wtmp,zid);
Wp = (Wslack-Wtmp)*beta;

Utmp = U-beta*(de1tmBRRt);
Up = max(0,Utmp);

Vtmp = V-beta*(de1tmBRsqr);
Vp = max(0,Vtmp);

f = Produc(C,RR1)+norm(Wp,'fro')^2/(2*beta) ...
    +norm(Up,'fro')^2/(2*beta)+norm(Vp,'fro')^2/(2*beta);

if Btype == 0
    eugf = C-Wp;
    eugf(1,1) = 0;% in fact no need
    eugf = KKtmap(eugf);
else
    tmp1 = dB'*Up*zI';
    tmp1 = (tmp1+tmp1')/2;
    %
    eugf = C-Wp-dB'*Vp*dB-tmp1;
    eugf(1,1) = 0;
    eugf = KKtmap(eugf);
end

%%
if comp_pfeas_flag == 1
    Uslack = (Up-Utmp)/beta;
    Vslack = (Vp-Vtmp)/beta;
    pfeas = max([ ...
        norm(RR1-Wslack,'fro')/(1+norm(RR1,'fro')), ...
        norm(de1tmBRRt-Uslack,'fro')/(1+norm(Uslack,'fro')), ...
        norm(de1tmBRsqr-Vslack,'fro')/(1+norm(Vslack,'fro')), ...
        ]);
else
    pfeas = inf;
end

end
