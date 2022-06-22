function res=GetOptimalH_W(S,cntr,ImpR,HH,flag,flagW)
% Synthesis of peak-to-peak controller
% cf. Section 3.3.3 of "Tight computationally efficient approximation of matrix
% norms with applications" by A. Juditsky, G. Kotsalis and A. Nemirovski
% See boeing_demo.m example
% Syntax:
%   res=GetOptimalH_W(Sys,Cntr,ImpR,H,flag,flagW)
% Inputs: 
%   Sys:    structure defining dynamical system with fields
%       nx: state vector dimension
%       ny: output vector dimension
%       nu: control vector dimension
%       nd: disturbance vector dimension
%       x2x: nxXnx state to state matrix
%       u2x: nxXnu control to state matrix
%       d2x: nxXnd disturbance to state matrix
%       d2y: nyXnd disturbance to output matrix
%   Cntr:   control parameters structure with fields
%       pd: in [1,2], norm exponent on the space of disturbances
%       px: in [1,2], norm exponent on the state space
%       py: in [1,2], norm exponent on the space of outputs
%       pu: in [1,2], norm exponent on the space of controls
%       xyu:[cx,cy,cu] nonegative vector, defines the optimization criterion
%           cx*Gain(d->x) + cy*Gain(d->y) + cu*Gain(d->u)')
%       d:  controller memory depth (whatever it is) 
%       T:  controller training horizon
%       W.A: diagonal matrix with nonnegative elements, matrix of constaints on individual gains
%       W.b: nonegative vector, r.h.s. of the system of constraints on individual gains
%           W.A*[Gain(d->x); Gain(d->y); Gain(d->u)]<=W.b
%   ImpR:   Impulse response structure built by FullImpulseResponse.m
%   H:      when nonempty, contains the controller celle array for which the peak-to-peak gain is computed    
%   flag:  {0,1}:   if 1, controller is synthesized, 
%                   if 0, the peak-to-peak gain is computed  for conroller in H  
%   flagW: ???? flag to check the gain constraints ????
%
% Output: 
%   res:  structure with fields
%       H: cell array containing synthesized controller
%       Hijs: 3d array containg synthesized controller
%       obj:optimal gain bound
%       outx:
%       outy:
%       outu:
%
% Author: A. Nemirovski (2022)

SMALL=0.01;
p=cntr.pd; d=cntr.d; T=cntr.T;
nx=S.nx; ny=S.ny; nd=S.nd; nu=S.nu;
% nz=nx*(d+1);
nh=nu*ny*d;
if isempty(HH)
    HH=zeros(nu,ny,d);
end
%
flagd=0;
if p==2
    ps=inf; flagd=1;
elseif p==inf
    ps=1;
else
    ps=p/(p-2);
end
if cntr.px==2
    rg.x=inf; dg.x=1;
else
    rg.x=cntr.px/(2-cntr.px);
    dg.x=0;
end
if cntr.py==2
    rg.y=inf;dg.y=1;
else
    rg.y=cntr.py/(2-cntr.py);
    dg.y=0;
end
if cntr.pu==2
    rg.u=inf;dg.u=1;
else
    rg.u=cntr.pu/(2-cntr.pu);
    dg.u=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% building CVX problem
cvx_begin;
variable H(nu,ny,d);
if flagd==0
    variable mlmx(T,nd);
    variable mlmy(T,nd);
    variable mlmu(T,nd);
    variable nmlmx(T,1);
    variable nmlmy(T,1);
    variable nmlmu(T,1);
else
    variable mlmx(T,1);
    variable mlmy(T,1);
    variable mlmu(T,1);
end 
if dg.x==0
    variable dgx(nx,1);
end
variable pupx(nx,nx,T);
if dg.y==0
    variable dgy(ny,1);
end
variable pupy(ny,ny,T);
if dg.u==0
    variable dgu(nu,1);
end
variable pupu(nu,nu,T);
variable outx(nx,nd,T+1);
variable outy(ny,nd,T);
variable outu(nu,nd,T);
variable obj(3,1);
variable tx;
variable ty;
variable tu;
variable objw;
expression totx;    
expression toty;
expression totu;
variable opt;
if flag==0
    H == HH;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outx == reshape(reshape(ImpR.x,nx*nd*(T+1),nh+1)*[reshape(H,nh,1,1);1],nx,nd,T+1);
for t=1:T,
    pupx(:,:,t) == pupx(:,:,t)';
end;
for t=1:T,
    if flagd,
        [pupx(:,:,t),outx(:,:,t+1);outx(:,:,t+1)',mlmx(t)*eye(nd)] == semidefinite(nx+nd);
    else
        [pupx(:,:,t),outx(:,:,t+1);outx(:,:,t+1)',diag(mlmx(t,:))] == semidefinite(nx+nd);
    end;
end;
totx=pupx(:,:,1);
for j=2:T,
    totx=totx+pupx(:,:,j);
end;
if rg.x==inf,
    tx >= norm(totx);
else
    diag(dgx)-totx == semidefinite(nx);
    tx >= norm(dgx,rg.x);
end;
if flagd,
    obj(1) >= 0.5*(tx+sum(mlmx));
else
    for t=1:T,
        nmlmx(t) >= norm(mlmx(t,:),ps);
    end;
    obj(1) >= 0.5*(tx+sum(nmlmx));
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outy == reshape(reshape(ImpR.y,ny*nd*T,nh+1)*[reshape(H,nh,1,1);1],ny,nd,T);
for t=1:T,
    pupy(:,:,t) == pupy(:,:,t)';
end;
for t=1:T,
    if flagd,
        [pupy(:,:,t),outy(:,:,t);outy(:,:,t)',mlmy(t)*eye(nd)] == semidefinite(ny+nd);
    else
        [pupy(:,:,t),outy(:,:,t);outy(:,:,t)',diag(mlmy(t,:))] == semidefinite(ny+nd);
    end;
end;
toty=pupy(:,:,1);
for j=2:T,
    toty=toty+pupy(:,:,j);
end;
if rg.y==inf,
    ty >= norm(toty);
else
    diag(dgy)-toty == semidefinite(ny);
    ty >= norm(dgy,rg.y);
end;
if flagd,
    obj(2) >= 0.5*(ty+sum(mlmy));
else
    for t=1:T,
        nmlmy(t) >= norm(mlmy(t,:),ps);
    end;
    obj(2) >= 0.5*(ty+sum(nmlmy));
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outu == reshape(reshape(ImpR.u,nu*nd*T,nh+1)*[reshape(H,nh,1,1);1],nu,nd,T);
for t=1:T,
    pupu(:,:,t) == pupu(:,:,t)';
end;
for t=1:T,
    if flagd,
        [pupu(:,:,t),outu(:,:,t);outu(:,:,t)',mlmu(t)*eye(nd)] == semidefinite(nu+nd);
    else
        [pupu(:,:,t),outu(:,:,t);outu(:,:,t)',diag(mlmu(t,:))] == semidefinite(nu+nd);
    end;
end;
totu=pupu(:,:,1);
for j=2:T,
    totu=totu+pupu(:,:,j);
end;
if rg.u==inf,
    tu >= norm(totu);
else
    diag(dgu)-totu == semidefinite(nu);
    tu >= norm(dgu,rg.u);
end;
if flagd,
    obj(3) >= 0.5*(tu+sum(mlmu));
else
    for t=1:T,
        nmlmu(t) >= norm(mlmu(t,:),ps);
    end;
    obj(3) >= 0.5*(tu+sum(nmlmu));
end;
if (~isempty(cntr.W.A))&flagW,
    cntr.W.A*obj <= cntr.W.b;
end;
opt >= cntr.xyu*obj;
minimize opt+SMALL*sum(obj);
cvx_end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res.H=cell(d,1);
for i=1:d,
    res.H{i}=full(H(:,:,i));
end;
res.status=cvx_status;
res.obj=obj;
res.outx=outx;
res.outy=outy;
res.outu=outu;
res.Hijs=H;