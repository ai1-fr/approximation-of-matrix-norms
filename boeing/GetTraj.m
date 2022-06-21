function traj=GetTraj(S,cntr,H,ds,T)
% trace sample trajectory
%
% Author: A. Nemirovski (2022)

ys=zeros(S.ny,T);
vs=zeros(S.ny,T);
xs=zeros(S.nx,T+1);
bxs=zeros(S.nx,T+1);
bys=zeros(S.ny,T);
us=zeros(S.nu,T);
for t=0:T
    it=t+1;
    if t==0
        xs(:,it)=0;
        bxs(:,it)=0;
    else
        xs(:,it)=S.x2x*xs(:,t)+S.u2x*us(:,t)+S.d2x*ds(:,t);
        bxs(:,it)=S.x2x*bxs(:,t)+S.u2x*us(:,t);     
    end
    if it<=T
        ys(:,it)=S.x2y*xs(:,it)+S.d2y*ds(:,it);
        bys(:,it)=S.x2y*bxs(:,it);
        vs(:,it)=ys(:,it)-bys(:,it);
        
        %%%
        u=zeros(S.nu,1);
        for tau=0:cntr.d-1
            itau=tau+1;
            if tau<=t
                u=u+H{itau}*vs(:,t-tau+1);
            end
        end
        us(:,it)=u;
    end
end
traj.x=xs(:,2:T+1);
traj.u=us;
traj.y=ys;
traj.by=bys;
traj.bx=bxs;