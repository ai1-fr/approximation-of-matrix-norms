% Peak-to-peak control simulation example 
% Section 3.3.3 of "Tight computationally efficient approximation of matrix
% norms with applications" by A. Juditsky, G. Kotsalis and A. Nemirovski
%
% Author: A. Nemirovski (2022)

clear all

cvx_solver mosek
cvx_quiet(true)

% generating Boeing747 data
[data,cntr]=GetBoeing747;
T=cntr.T;
ImpR=FullImpulseResponse(data,cntr);

% building p2p controller on horizon T=256
disp('Building POB controller...')
disp(' ')
tstart=clock;
POBT=GetOptimalH_W(data,cntr,ImpR,[],1,1);
tused=etime(clock,tstart);
fprintf('CPU=%5.1f, POB controller built on  horizon T=%d\n',tused,T)
fprintf('Peak2peak gains: states=%.3e; outputs=%.3e; controls=%.3e\n',POBT.obj(1),POBT.obj(2),POBT.obj(3));

% checking p2p controller on horizon T=512
cntr2=cntr; T2=2*T; cntr2.T=T2;

ImpRR=FullImpulseResponse(data, cntr2);
HH=POBT.Hijs;
disp(' ')
disp('Checking POB controller...')
disp(' ')
tstart=clock;
POBTT=GetOptimalH_W(data,cntr2,ImpRR,HH,0,1);
tused=etime(clock,tstart);
fprintf('CPU=%5.1f, POB controler checked on horizon T=%d\n',tused,T2)
fprintf('Peak2peak gains: states=%.3e; outputs=%.3e; controls=%.3e\n',POBTT.obj(1),POBTT.obj(2),POBTT.obj(3));

% checking trivial controller on horizon T=512
disp(' ')
disp('Checking trivial controller...')
disp(' ')
if max(abs(eig(data.x2x)))<1
    HH=zeros(data.nu,data.ny,cntr.d);
    tstart=clock;
    ZROTT=GetOptimalH_W(data,cntr2,ImpRR,HH,0,0);
    tused=etime(clock,tstart);
    fprintf('CPU=%5.1f, Trivial controler checked on horizon T=%d\n',  tused,T2)
    fprintf('Peak2peak gains: states=%.3e; outputs=%.3e; controls=%.3e\n',ZROTT.obj(1),ZROTT.obj(2),ZROTT.obj(3));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Drawing figures
close all
%Flagzero=1; Flag='y'; FFlagShow='y';

ZRH=cell(cntr.d,1);
for i=1:cntr.d
    ZRH{i}=zeros(data.nu,data.ny);
end
ntrajH.x=zeros(T,1);
ntrajH.y=zeros(T,1);
ntrajH.u=zeros(T,1);
% plots of random disturbances
% FlagB=0
z=zeros(data.nx,1);
ds=randn(data.nd,T);
for t=1:T
    zz=data.d2x'*data.x2x*z;
    nzz=norm(zz,cntr.pd);
    if nzz==0
        e=randn(data.nd,1);
        e=e/norm(e,cntr.pd);
    else
        e=zz/nzz;
    end
    ds(:,t)=e;
    z=data.x2x*z+data.d2x*ds(:,t)+0.1*randn(data.nx,1);
end
for t=1:T
    ds(:,t)=ds(:,t)/norm(ds(:,t),cntr.pd);
end
trajH=GetTraj(data,cntr,POBT.H,ds,T);
for t=1:T
    ntrajH.x(t)=norm(trajH.x(:,t),cntr.px);
    ntrajH.y(t)=norm(trajH.y(:,t),cntr.py);
    ntrajH.u(t)=norm(trajH.u(:,t),cntr.pu);
end
trajZ=GetTraj(data,cntr,ZRH,ds,T);
for t=1:T
    ntrajZ.x(t)=norm(trajZ.x(:,t),cntr.px);
    ntrajZ.y(t)=norm(trajZ.y(:,t),cntr.py);
    ntrajZ.u(t)=norm(trajZ.u(:,t),cntr.pu);
end

r=groot; 
scrs=r.MonitorPositions(3:4);

ff1=figure('name','  1. System with random disturbances',...
    'NumberTitle','off');
ff1.Position=[(scrs(1)-760)/2+10, (scrs(2)-720)/2+10, 860, 720];
tiledlayout(3,2)
% tl plot
ax1 = nexttile;
p0=plot(ax1,1:1:T,ntrajH.x,'-g','LineWidth',1);
hold on
for i=1:data.nx
    p(i)=plot(ax1,1:1:T,trajH.x(i,:),'-b');
end
grid on
legend([p(1), p0], {'states', 'state norm'}, 'Location','southeast')
title(ax1,'POB states')
xlim(ax1, [1,256])
hold off

% tr plot
ax2 = nexttile;
plot(ax2,1:1:T,ntrajZ.x,'-g','LineWidth',1)
hold on
for i=1:data.nx
    plot(ax2,1:1:T,trajZ.x(i,:),'-b');
end
grid on
title(ax2,'Trivial states')
xlim(ax2, [1,256])
hold off

% ml plot
ax3 = nexttile;
plot(ax3, 1:1:T,ntrajH.u,'-g','LineWidth',1);
hold on;
for i=1:data.nu
    plot(ax3, 1:1:T,trajH.u(i,:),'-r');
end
grid on;
title(ax3,'POB controls')
xlim(ax3, [1,256])
hold off

% mr plot
ax4 = nexttile;
plot(ax4, 1:1:T,ntrajZ.u,'-g','LineWidth',1);
hold on
for i=1:data.nu
    plot(ax4, 1:1:T,trajZ.u(i,:),'-r')
end
grid on;
title(ax4,'Trivial controls')
xlim(ax4, [1,256])
hold off

% bl plot
ax5 = nexttile;
plot(ax5,1:1:T,ntrajH.y,'-g','LineWidth',1);
hold on;
for i=1:data.ny
    plot(ax5,1:1:T,trajH.y(i,:),'-c');
end
grid on
title(ax5,'POB outputs')
xlim(ax5, [1,256])
hold off

% br plot
ax6 = nexttile;
plot(ax6,1:1:T,ntrajZ.y,'-g','LineWidth',1);
hold on
for i=1:data.ny
    plot(ax6,1:1:T,trajZ.y(i,:),'-c');
end
grid on;
title(ax6,'Trivial outputs')
xlim(ax6, [1,256])
hold off
disp(' ')
disp('Actual gains, random disturbances')
fprintf('states: %.3e/%.3e; outputs: %.3e/%.3e; controls: %.3e/%.3e\n',...
            max(ntrajH.x),max(ntrajZ.x),max(ntrajH.y),max(ntrajZ.y),max(ntrajH.u),max(ntrajZ.u))

%plots of worst-case disturbances
ds=GetBadDs(data,cntr,T);
trajH=GetTraj(data,cntr,POBT.H,ds,T);
for t=1:T
    ntrajH.x(t)=norm(trajH.x(:,t),cntr.px);
    ntrajH.y(t)=norm(trajH.y(:,t),cntr.py);
    ntrajH.u(t)=norm(trajH.u(:,t),cntr.pu);
end
trajZ=GetTraj(data,cntr,ZRH,ds,T);
for t=1:T
    ntrajZ.x(t)=norm(trajZ.x(:,t),cntr.px);
    ntrajZ.y(t)=norm(trajZ.y(:,t),cntr.py);
    ntrajZ.u(t)=norm(trajZ.u(:,t),cntr.pu);
end
ff2=figure('name','  2. System with worst-case disturbances',...
    'NumberTitle','off');
ff2.Position=[(scrs(1)-760)/2, (scrs(2)-720)/2, 860, 720];
tiledlayout(3,2)
% tl plot
ax1 = nexttile;
p0=plot(ax1,1:1:T,ntrajH.x,'-g','LineWidth',1);
hold on
for i=1:data.nx
    p(i)=plot(ax1,1:1:T,trajH.x(i,:),'-b');
end
grid on
legend([p(1), p0], {'states', 'state norm'}, 'Location','southeast')
title(ax1,'POB states')
xlim(ax1, [1,256])
hold off

% tr plot
ax2 = nexttile;
plot(ax2,1:1:T,ntrajZ.x,'-g','LineWidth',1)
hold on
for i=1:data.nx
    plot(ax2,1:1:T,trajZ.x(i,:),'-b');
end
grid on
title(ax2,'Trivial states')
xlim(ax2, [1,256])
hold off

% ml plot
ax3 = nexttile;
plot(ax3, 1:1:T,ntrajH.u,'-g','LineWidth',1);
hold on;
for i=1:data.nu
    plot(ax3,1:1:T,trajH.u(i,:),'-r');
end
grid on;
title(ax3,'POB controls')
xlim(ax3, [1,256])
hold off

% mr plot
ax4 = nexttile;
plot(ax4, 1:1:T,ntrajZ.u,'-g','LineWidth',1);
hold on
for i=1:data.nu
    plot(ax4, 1:1:T,trajZ.u(i,:),'-r')
end
grid on;
title(ax4,'Trivial controls')
xlim(ax4, [1,256])
hold off

% bl plot
ax5 = nexttile;
plot(ax5,1:1:T,ntrajH.y,'-g','LineWidth',1);
hold on;
for i=1:data.ny
    plot(ax5,1:1:T,trajH.y(i,:),'-c');
end
grid on
title(ax5,'POB outputs')
xlim(ax5, [1,256])
hold off

% br plot
ax6 = nexttile;
plot(ax6,1:1:T,ntrajZ.y,'-g','LineWidth',1);
hold on
for i=1:data.ny
    plot(ax6,1:1:T,trajZ.y(i,:),'-c');
end
grid on;
title(ax6,'Trivial outputs')
xlim(ax6, [1,256])
hold off

disp(' ')
disp('Actual gains, worst-case disturbances')
fprintf('states: %.3e/%.3e; outputs: %.3e/%.3e; controls: %.3e/%.3e\n',...
            max(ntrajH.x),max(ntrajZ.x),max(ntrajH.y),max(ntrajZ.y),max(ntrajH.u),max(ntrajZ.u))
