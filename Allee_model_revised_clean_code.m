%%  delayed Allee paper
clear; close all
base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
ddebiftool_path(base,'sym');
format compact; format short

%% ========================================================================
% Parameters and model
% =========================================================================
parnames = {'a','b','lambda','tau'};
cind = [parnames; num2cell(1:length(parnames))];
in = struct(cind{:});
par([in.a,in.b,in.lambda,in.tau]) = [1, 0.5, 0, 3.5];
a = par(in.a);  b = par(in.b);
beta  = a-b;    Uf = sqrt(beta/3);
lam_f = (2/(3*sqrt(3)))*beta^(3/2);
tau_BT = 1/b;
fprintf('analytic folds: lambda = +/- %.8f,  BT delay tau = %.8f\n',lam_f,tau_BT);

funcs = set_symfuncs(@sym_allee_model,'sys_tau',@()in.tau);
getp  = @(br,name) arrayfun(@(p)p.parameter(in.(name)),br.point);
X = struct();                                   % collected results

%% ========================================================================
% (A) Equilibria and PO branches in lambda at tau = 3.5
% =========================================================================
par(in.lambda) = -1;  par(in.tau) = 3.5;
rr = roots([1,0,-beta,-par(in.lambda)]);  
u0 = real(rr(abs(imag(rr))<1e-8));
[eqlam,~] = SetupStst(funcs,'parameter',par,'x',u0(1),'contpar',in.lambda,...
    'step',0.05,'min_bound',[in.lambda,-1],'max_bound',[in.lambda,1],...
    'max_step',[in.lambda,0.05;0,0.05]);
eqlam = br_contn(funcs,eqlam,3000); 
eqlam = br_rvers(eqlam);
eqlam = br_contn(funcs,eqlam,3000);
[eqlam_ref,nunst_lam,~,p_bif_lam] = MonitorChange(funcs,eqlam,'min_iterations',5);
lam_eq = getp(eqlam_ref,'lambda');
u_eq   = arrayfun(@(p)p.x(1),eqlam_ref.point);

is_fold     = abs(abs(lam_eq(p_bif_lam))-lam_f) < 5e-3;
foldidx_lam = p_bif_lam(is_fold);   [~,o] = sort(lam_eq(foldidx_lam)); foldidx_lam = foldidx_lam(o);
hopfidx_lam = p_bif_lam(~is_fold);  [~,o] = sort(lam_eq(hopfidx_lam)); hopfidx_lam = hopfidx_lam(o);
fprintf('folds at lambda = %s;  Hopf at lambda = %s\n',...
    mat2str(lam_eq(foldidx_lam),6),mat2str(lam_eq(hopfidx_lam),6));

psol_lam = cell(1,numel(hopfidx_lam));
for j = 1:numel(hopfidx_lam)
    [funcspo1,pb,suc] = SetupPsol(funcs,eqlam_ref,hopfidx_lam(j),'degree',6,'intervals',500,...
        'max_step',[in.lambda,0.05],'min_bound',[in.lambda,-0.2],...
        'max_bound',[in.lambda,0.2],'step',0.02,'outputfuncs',true);
    if ~suc, warning('SetupPsol failed at lambda-Hopf %d',j); continue; end
    [pb,succpo1] = br_contn(funcspo1,pb,34);
    [pb_ref,nun,~,pbf] = MonitorChange(funcspo1,pb,'range',2:length(pb.point),...
        'exclude_trival',true,'min_iterations',5);
    psol_lam{j} = struct('branch',pb_ref,'nunst',nun,'p_bif',pbf,...
        'lam',getp(pb_ref,'lambda'),...
        'max',arrayfun(@(p)max(p.profile(1,:)),pb_ref.point),...
        'min',arrayfun(@(p)min(p.profile(1,:)),pb_ref.point),...
        'per',arrayfun(@(p)p.period,pb_ref.point));
end

%% ========================================================================
% (B) Equilibria and Hopf-born PO branch in tau at lambda = -0.05
% =========================================================================
par(in.lambda) = -0.05;  par(in.tau) = 0;
ur = sort(real(roots([1,0,-beta,-par(in.lambda)])));
u_low = ur(1);  u_mid = ur(2);  u_up = ur(3);
Aup = a-3*u_up^2;
tau_H_analytic = acos(Aup/b)/sqrt(b^2-Aup^2);
fprintf('lambda=-0.05: u_low=%.6f u_mid=%.6f u_up=%.6f, tau_H=%.6f\n',...
    u_low,u_mid,u_up,tau_H_analytic);

[eqtau,~] = SetupStst(funcs,'parameter',par,'x',u_up,'contpar',in.tau,'step',0.05,...
    'min_bound',[in.tau,0],'max_bound',[in.tau,10],'max_step',[in.tau,0.05;0,0.05]);
eqtau = br_contn(funcs,eqtau,3000);
[eqtau_ref,nunst_tau,~,p_bif_tau] = MonitorChange(funcs,eqtau,'min_iterations',5);
tau_eq = getp(eqtau_ref,'tau');
u_eqB  = arrayfun(@(p)p.x(1),eqtau_ref.point);
[~,jj] = min(tau_eq(p_bif_tau));  hopf_tau_idx = p_bif_tau(jj);
fprintf('numerical tau_H = %.8f (analytic %.8f)\n',tau_eq(hopf_tau_idx),tau_H_analytic);

[funcspo2,pbt,suc] = SetupPsol(funcs,eqtau_ref,hopf_tau_idx,'degree',6,'intervals',120,...
    'max_step',[in.tau,0.05],'max_bound',[in.tau,10],'min_bound',[in.tau,0],'step',0.1,...
    'outputfuncs',true);
if ~suc, error('SetupPsol failed at the tau-Hopf point'); end
[pbt,sucpo2] = br_contn(funcspo2,pbt,50);
[pbt_ref,nunst_po_tau,~,~] = MonitorChange(funcspo2,pbt,'range',2:length(pbt.point),...
    'exclude_trival',true,'min_iterations',5);
tau_po     = getp(pbt_ref,'tau');
max_po_tau = arrayfun(@(p)max(p.profile(1,:)),pbt_ref.point);
min_po_tau = arrayfun(@(p)min(p.profile(1,:)),pbt_ref.point);
per_po_tau = arrayfun(@(p)p.period,pbt_ref.point);
fprintf('Hopf-born PO branch: tau in [%.4f, %.4f], T in [%.2f, %.2f]\n',...
    min(tau_po),max(tau_po),min(per_po_tau),max(per_po_tau));

%% ========================================================================
% (C) Two-parameter Hopf and fold curves, BT point
% =========================================================================
parbd2 = {'min_bound',[in.lambda,-0.2;in.tau,0],'max_bound',[in.lambda,0.2;in.tau,12],...
          'max_step',[in.lambda,0.01;in.tau,0.1;0,0.1]};
[fhopf,hopf_br,~] = SetupHopf(funcs,eqtau_ref,hopf_tau_idx,'contpar',[in.lambda,in.tau],...
    'dir',in.tau,'step',0.05,'outputfuncs',true,parbd2{:});
hopf_br = br_contn(fhopf,hopf_br,300);  hopf_br = br_rvers(hopf_br);
hopf_br = br_contn(fhopf,hopf_br,300);
[hopf_br,hopf_tests,hopf_indbif,hopf_biftype] = LocateSpecialPoints(fhopf,hopf_br);
lam_hopf = getp(hopf_br,'lambda');  tau_hopf = getp(hopf_br,'tau');
om_hopf  = arrayfun(@(p)p.omega,hopf_br.point);
L1_hopf  = hopf_tests.genh(1,:);

[~,jf] = min(lam_eq(foldidx_lam));
[ffold,fold_br,~] = SetupFold(funcs,eqlam_ref,foldidx_lam(jf),'contpar',[in.lambda,in.tau],...
    'dir',in.tau,'step',0.05,'outputfuncs',true,parbd2{:});
fold_br = br_contn(ffold,fold_br,300); 
fold_br = br_rvers(fold_br);
fold_br = br_contn(ffold,fold_br,300);
lam_fold = getp(fold_br,'lambda'); 
tau_fold = getp(fold_br,'tau');

% BT point and normal-form coefficients (paper: a2=-6bU_f, b2=-4U_f)
indBT = hopf_indbif(strcmp(hopf_biftype,'BT'));  btbr = hopf_br;
if isempty(indBT)
    [fold_br,~,fold_indbif,fold_biftype] = LocateSpecialPoints(ffold,fold_br);
    indBT = fold_indbif(strcmp(fold_biftype,'BT'));  btbr = fold_br;
end
if ~isempty(indBT)
    X.bt = btbr.point(indBT(1));
    fprintf('BT at (lambda,tau) = (%.8f, %.8f)\n',X.bt.parameter(in.lambda),X.bt.parameter(in.tau));
    if isfield(X.bt,'nmfm'), disp(X.bt.nmfm); end
end
fprintf('analytic a2 = %.4f, b2 = %.4f\n',-6*b*Uf,-4*Uf);

%% ========================================================================
% (D) Large-amplitude cycle: seed, branch, fold of POs, two-parameter fold
% =========================================================================
% dde23 is used only to obtain a starting orbit (Poincare section u=u_mid).
p0 = pbt_ref.point(1).parameter;  p0(in.tau) = 4.5;
opts = ddeset('RelTol',1e-8,'AbsTol',1e-9,'Events',@(t,y,Z)section_event(y,u_mid));
seed = dde23(@(t,y,Z)funcs.wrap_rhs(cat(2,y,Z),p0),p0(in.tau),0.15,[0,1500],opts);

[funcspo,brL,suc] = branch_from_simulations(funcs,seed,in.tau,p0,...
    'degree',6,'intervals',120,'step',-0.01,'extra_condition',0,'outputfuncs',true,...
    'min_bound',[in.tau,3],'max_bound',[in.tau,6.5],'max_step',[in.tau,0.05]);
fprintf('large cycle seed: success = %d, T = %.4f\n',suc,brL.point(1).period);
brL = br_contn(funcspo,brL,250);  brL = br_rvers(brL);  brL = br_contn(funcspo,brL,120);
[nun,~,~,brL.point] = GetStability(brL,'funcs',funcs,'exclude_trivial',true);
TT = arrayfun(@(p)p.period,brL.point);
ok = TT < 200;                                  % drop near-homoclinic points
brL.point = brL.point(ok);
X.large_br    = brL;
X.large_nunst = nun(ok);
X.large_tau   = getp(brL,'tau');
X.large_T     = TT(ok);
X.large_max   = arrayfun(@(p)max(p.profile(1,:)),brL.point);
X.large_min   = arrayfun(@(p)min(p.profile(1,:)),brL.point);
[X.tau_PO,indfold] = min(X.large_tau);          % turning point = fold of POs
fprintf('fold of periodic orbits: tau_PO = %.5f\n',X.tau_PO);

[pffuncs,pfbr] = SetupPOfold(funcs,brL,indfold,'contpar',[in.tau,in.lambda],...
    'dir',in.lambda,'step',2e-3,'max_step',[in.lambda,0.005;in.tau,0.05],...
    'min_bound',[in.lambda,-0.2;in.tau,1/b],'max_bound',[in.lambda,0.2;in.tau,12],'outputfuncs',true);
pfbr.method.continuation.plot = 1;
pfbr = br_contn(pffuncs,pfbr,200); 
pfbr = br_rvers(pfbr);
pfbr = br_contn(pffuncs,pfbr,200);
X.pofold_lam = getp(pfbr,'lambda');  X.pofold_tau = getp(pfbr,'tau');
fprintf('PO fold curve: min tau = %.4f at lambda = %.4f\n',...
    min(X.pofold_tau),X.pofold_lam(X.pofold_tau==min(X.pofold_tau)));

%% ========================================================================
% (E) Homoclinic curve, second homoclinic, saddle eigenvalues
% =========================================================================
% tau_hom(lambda): the upper branch of the unstable manifold of U_mid
% converges to U_up for tau < tau_hom and leaves its basin for tau > tau_hom.
lam_list = [-0.1358 -0.1355 -0.135 -0.1325 -0.13 -0.125 -0.12 -0.11 -0.10 ...
            -0.09 -0.08 -0.07 -0.06 -0.05 -0.04 -0.03 -0.02 -0.015 -0.01];
tau_hom = nan(size(lam_list));
parfor k = 1:numel(lam_list)
    tau_hom(k) = tau_hom_bisect(lam_list(k),a,b,24,3000);
end
X.hom = struct('lam',num2cell(lam_list),'tau_fit',num2cell(tau_hom));
disp([lam_list(:),tau_hom(:)])
nu = linspace(0,1.5,200);                       % BT asymptotics (paper eq. homasym)
X.hom_asym_tau = (1+nu)/b;
X.hom_asym_lam = -lam_f + 147/100*b^2/Uf*nu.^2;

k05 = find(abs(lam_list+0.05)<1e-9,1);
X.tau_hom = tau_hom(k05);
X.tau_hom2 = tau_hom2_bisect(-0.05,a,b,4.25,4.40,24,2000);
fprintf('lambda=-0.05: tau_hom = %.5f, tau_hom2 = %.5f\n',X.tau_hom,X.tau_hom2);

pmid = dde_stst_create('parameter',pbt_ref.point(1).parameter,'x',u_mid);
pmid.parameter(in.tau) = X.tau_hom;
pmid.stability = p_stabil(funcs,pmid,eqtau_ref.method.stability);
ev = pmid.stability.l1;
X.sig_u = max(real(ev));  X.sig_s = max(real(ev(real(ev)<0)));
d  = tau_po - X.tau_hom;  m = d > 0 & d < 0.1;
cf = polyfit(log(d(m)),per_po_tau(m),1);
X.loglaw = struct('tau_hom',X.tau_hom,'slope',cf(1),'c',cf(2),'mask',m);
fprintf('saddle: sigma_u=%.4f sigma_s=%.4f; log-law slope %.3f (theory %.3f, %d pts)\n',...
    X.sig_u,X.sig_s,cf(1),-1/abs(X.sig_s),nnz(m));

%% ========================================================================
% (F) Floquet multipliers and threshold experiments
% =========================================================================
[nunst_po,dom_po] = GetStability(pbt_ref,'funcs',funcs,'exclude_trivial',true);
simopts = ddeset('RelTol',1e-8,'AbsTol',1e-9);
sims = struct();

% F1: perturbed stable upper equilibrium, point picked from eqtau_ref
[~,ie] = min(abs(tau_eq-1.0));  peq = eqtau_ref.point(ie);
sims.stable = dde23(@(t,y,Z)funcs.wrap_rhs(cat(2,y,Z),peq.parameter),...
    peq.parameter(in.tau),peq.x+0.03,[0,150],simopts);
sims.stable_point = peq;

% F2: periodic orbit picked from pbt_ref; perturbations +/- eps
X.eps_list = [0.005 0.01 0.02 0.03 0.05 0.1];
X.eps_tau  = [4.315, 4.0];
X.eps_out  = zeros(2*numel(X.eps_list),numel(X.eps_tau));
for j = 1:numel(X.eps_tau)
    [~,i] = min(abs(tau_po-X.eps_tau(j)));  po = pbt_ref.point(i);
    [~,ie] = min(abs(tau_eq-tau_po(i)));    ue = eqtau_ref.point(ie).x(1);
    X.floquet(j) = dom_po(i);
    fprintf('tau=%.4f: nunst=%d, dominant multiplier=%s\n',tau_po(i),nunst_po(i),num2str(dom_po(i)));
    hp = @(th) po_history(po,th);
    f  = @(t,y,Z) funcs.wrap_rhs(cat(2,y,Z),po.parameter);
    r = 0;
    for e = X.eps_list
        for s = [1,-1]
            r = r+1;
            sol = dde23(f,tau_po(i),@(th)ue+(1+s*e)*(hp(th)-ue),[0,600],simopts);
            X.eps_out(r,j) = fate(sol,u_mid);
            if j==1 && e==0.03               % runs shown in the threshold figure
                if s==1, sims.po_plus = sol; else, sims.po_minus = sol; end
            end
        end
    end
    if j==1, sims.po_point = po; sims.po_ueq = ue; end
end
disp('eps scan (rows +eps,-eps; cols tau): +1 U_up, -1 U_low, 0 large cycle');
disp([kron(X.eps_list(:),[1;1]).*repmat([1;-1],numel(X.eps_list),1), X.eps_out])

%% ========================================================================
% (G) Closed-form first Lyapunov coefficient (paper eq. c1)
% =========================================================================
Us = linspace(Uf+1e-5,sqrt((a+b)/3)-1e-5,3000);
A  = a-3*Us.^2;  w = sqrt(b^2-A.^2);  tH = acos(A/b)./w;
c1 = (-3+36*Us.^2./(b-A)+18*Us.^2./((b-A)/b.*(2i*w-2*A-b)))./(1-tH.*A+1i*tH.*w);
X.l1_tau = tH;  X.l1_val = real(c1)./w;
X.l1_asym = 2*beta./((b-A).*sqrt(2*b*(b-A)));
fprintf('min l1 on the Hopf curve = %.4f\n',min(X.l1_val));

%% ========================================================================
% (H) Basin map of constant histories and critical level (lambda = -0.05)
% =========================================================================
X.basin_umid = u_mid;  X.basin_uup = u_up;
X.basin_tau  = 3.0:0.02:5.4;
X.basin_c    = linspace(u_mid+0.002,u_up-0.002,200);
X.basin      = zeros(numel(X.basin_c),numel(X.basin_tau));
for it = 1:numel(X.basin_tau)
    X.basin(:,it) = fate_rk4(X.basin_c,-0.05,X.basin_tau(it),a,b,u_mid,150,900).';
end
X.ccrit_tau = 3.20:0.02:4.30;  X.ccrit = nan(size(X.ccrit_tau));
cfine = linspace(u_mid+1e-4,u_up-1e-4,2000);
for k = 1:numel(X.ccrit_tau)
    fk = fate_rk4(cfine,-0.05,X.ccrit_tau(k),a,b,u_mid,200,900);
    j = find(fk==1,1,'first');
    if ~isempty(j) && ~any(fk==0), X.ccrit(k) = cfine(j); end
end
disp([X.ccrit_tau(:),X.ccrit(:)])

%% ========================================================================
save('allee_results_final_23October.mat')
fprintf('saved allee_results.mat\n');

%% ======================== local functions ===============================
function u = po_history(po,th)
% periodic history from a DDE-Biftool psol point (time in physical units)
[t,iu] = unique(po.mesh*po.period);
u = interp1(t,po.profile(1,iu),mod(th,po.period),'pchip');
end

function [val,isterm,dir] = section_event(y,lev)  % dde23 event: u = lev, upward
val = y(1)-lev;  isterm = 0;  dir = 1;
end

function c = fate(sol,um)            % +1 upper, -1 lower, 0 large cycle
t = linspace(sol.x(end)-250,sol.x(end),2000); u = deval(sol,t);
if max(u)-min(u) > 0.5, c = 0; elseif u(end) > um, c = 1; else, c = -1; end
end

function c = wu_fate(lam,tau,a,b,A,um,opts,Tend)
% fate of the upper branch of the unstable manifold of U_mid
su  = fzero(@(s) s-A+b*exp(-s*tau),1);
sol = dde23(@(t,y,Z) a*y-y^3+lam-b*Z,tau,@(t) um+1e-6*exp(su*t),[0,Tend],opts);
c   = fate(sol,um);
end

function th = tau_hom_bisect(lam,a,b,nit,Tend)
% homoclinic delay: upper W^u(U_mid) switches from U_up (+1) to not U_up
r = sort(real(roots([1,0,-(a-b),-lam])));  um = r(2);  A = a-3*um^2;
Au = a-3*r(3)^2;  tH = acos(Au/b)/sqrt(b^2-Au^2);
lo = (1/b)*(1+1e-4);  hi = tH-1e-3;
opts = ddeset('RelTol',1e-9,'AbsTol',1e-11);
side = @(tau) wu_fate(lam,tau,a,b,A,um,opts,Tend);
if side(lo)~=1 || side(hi)==1, th = NaN; return; end
for n = 1:nit, m = (lo+hi)/2; if side(m)==1, lo = m; else, hi = m; end, end
th = (lo+hi)/2;
end

function th = tau_hom2_bisect(lam,a,b,lo,hi,nit,Tend)
% second homoclinic: upper W^u(U_mid) switches from U_low (-1) to the large cycle (0)
r = sort(real(roots([1,0,-(a-b),-lam])));  um = r(2);  A = a-3*um^2;
opts = ddeset('RelTol',1e-9,'AbsTol',1e-11);
side = @(tau) wu_fate(lam,tau,a,b,A,um,opts,Tend);
if side(lo)~=-1 || side(hi)~=0, th = NaN; return; end
for n = 1:nit, m = (lo+hi)/2; if side(m)==-1, lo = m; else, hi = m; end, end
th = (lo+hi)/2;
end

function cls = fate_rk4(cvec,lam,tau,a,b,um,N,Tend)
% fate of all constant histories u=c at once (vectorised RK4, step tau/N)
c = cvec(:).';  M = numel(c);  dt = tau/N;  ns = round(Tend/dt);
B = repmat(c,N+1,1);  cur = N+1;
mx = -inf(1,M);  mn = inf(1,M);  nwin = ns-round(250/dt);
f = @(x,xd) a*x - x.^3 + lam - b*xd;
for n = 1:ns
    old = mod(cur,N+1)+1;  nxt = mod(old,N+1)+1;
    x = B(cur,:);  d0 = B(old,:);  d1 = B(nxt,:);  dm = 0.5*(d0+d1);
    k1 = f(x,d0); k2 = f(x+0.5*dt*k1,dm); k3 = f(x+0.5*dt*k2,dm); k4 = f(x+dt*k3,d1);
    xn = x + dt/6*(k1+2*k2+2*k3+k4);
    B(old,:) = xn;  cur = old;
    if n > nwin, mx = max(mx,xn); mn = min(mn,xn); end
end
cls = sign(xn-um);  cls(mx-mn > 0.5) = 0;
end