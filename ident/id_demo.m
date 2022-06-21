% Robust recovery example for system 
%   u(t+1)=X*[u(t);r(t)]
% via noisy observations of states u(t) (dimension nu) and inputs r(t) (dimension nr) for t=0,1,...,N-1
% observations of entries in u(t) and r(t) are subject to additive noise
% of absolute or semi-relative magnitude epsilon.
%
% Section 4.3.2 of "Tight computationally efficient approximation of matrix
% norms with applications" by A. Juditsky, G. Kotsalis and A. Nemirovski
%
% Author: A. Nemirovski (2022)

clear all
cvx_solver mosek
cvx_quiet(true)

verbose=0;
eps_s=0.01;

data=BoeingIdData(eps_s);
disp(['Boeing system identification, eps=',num2str(eps_s)])
disp(' ')
disp('Minimizing Gamma over E...')
disp(' ')

tstart=clock;
UpsBnd=inf;
GammaBest=inf;
RiskBest=inf;
xRec=zeros(data.nx,1);
itr=0;
flagOk=1;
UpsMax=1;
Coeffs=inf*ones(20,1);
fprintf(' CPU  |  Gamma  | Ups min | Ups max | RiskUB  |\n')
while(1==1)
    itr=itr+1;
    resE=Upsilons(data,UpsBnd,verbose);
    if resE.riskUB<RiskBest
        RiskBest=resE.riskUB;
%        fprintf('RiskUB set to %5.4f\n',RiskBest);
        xRec=resE.x;
    end
    tused=etime(clock,tstart);
    fprintf('%5.1f | %7.4f |',tused,resE.Gamma)
    if resE.Gamma==inf
        disp('Failure')
        flagOk=0;
        break;
    end
    if itr==1
        UpsilonStar=resE.Upsi;
        UpsMin=resE.Upsi;
    end
    Coeff(itr)=resE.Gamma;
    if resE.Gamma<GammaBest
        GammaBest=resE.Gamma;
        E=resE.E;
        UpsilonE=resE.Upsi;
        UpsMin=UpsilonE;
    else
        UpsMax=resE.Upsi;
    end
    fprintf('  %5.4f |  %5.4f |  %5.4f |\n',UpsMin,UpsMax, RiskBest)
    if UpsMax-UpsMin<1.e-2
        break    
    end
    UpsBnd=0.5*(UpsMin+UpsMax);
end
if flagOk
    coeff=GammaBest;
    disp(' ')
    disp('Specifying H...');
    tstart=clock;
    resH=GetH(data,coeff);
    tused=etime(clock,tstart);
    if verbose
        fprintf('CPU=%5.1f; CVX status %s; Risk=%5.4f\n',tused,resH.status,resH.obj)
    end
    % recovery of Bx in BxRecovered;
    BxRecovered=data.bary+resH.H'*data.barq;
    Err=data.B*data.x-BxRecovered;
    RecErr.two=norm(Err);
    RecErr.one=norm(Err,1);
    RecErr.inf=norm(Err,inf);
    disp(' ')
    fprintf('Actual norms of recovery error: l1=%5.4f l2=%5.4f linf=%5.4f\n',...
        RecErr.one,RecErr.two,RecErr.inf)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw figures
nu=data.nu;
nr=data.nr;
nx=data.nx;
Xrec=zeros(nu,nu+nr);
for Xell=1:nx,
    i=data.Xell2ij(Xell,1);
    j=data.Xell2ij(Xell,2);
    Xrec(i,j)=xRec(Xell);
end
T=256;
trajT=randn(nu+nr,T+1);
trajR=trajT;
for t=1:T
    trajT(1:nu,t+1)=data.X*trajT(:,t);
    trajR(1:nu,t+1)=Xrec*trajR(:,t);
end
close all
r=groot; 
scrs=r.MonitorPositions(3:4);
ff=figure('name','  True and recovered system trajectories',...
    'NumberTitle','off');
ff.Position=[(scrs(1)-560)/2, (scrs(2)-820)/2, 560,820];
tiledlayout(nu,1)
for i=1:nu
    ax = nexttile;
    plot(ax,0:1:T,trajT(i,:),':b',0:1:T,trajR(i,:),'--r')
    grid on
    xlim(ax, [1,256])
    xs=['$x_',num2str(i),'$'];
    title(ax,xs,'interpreter','latex')
    if i==1
        legend('true','recovered', 'Location','northwest')
    end
end
    
    
    