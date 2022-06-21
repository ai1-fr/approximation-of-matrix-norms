function ImpR=FullImpulseResponse(S,cntr)
% Syntax:
%   ImpR=FullImpulseResponse(Sys,Cntr)
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
% Output: 
%   ImpR:   Impulse response structure 
%
% Author: A. Nemirovski (2022)
nx=S.nx;
ny=S.ny;
nd=S.nd;
nu=S.nu;
d=cntr.d;
T=cntr.T;
x2x=S.x2x;
u2x=S.u2x;
d2x=S.d2x;
x2y=S.x2y;
d2y=S.d2y;
nz=(d+1)*nx;
nh=nu*ny*d;
%%%
ImpR.F=zeros(nz,nd,T+1,nh+1);
z2z=zeros(nz,nz);
for ir=1:d+1
    z2z((ir-1)*nx+1:ir*nx,(ir-1)*nx+1:ir*nx)=x2x;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=zeros(nu,ny,d);
basec=nx;
for ig=1:d
    z2z(1:nx,basec+1:basec+nx)=u2x*H(:,:,ig)*x2y;
    basec=basec+nx;
end
%%%
R=zeros(nx,nd,d);
R(:,:,1)=u2x*H(:,:,1)*d2y+d2x;
for ig=2:d,
    R(:,:,ig)=u2x*H(:,:,ig)*d2y;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for idl=1:nd,
    Delta=zeros(nd,1);
    Delta(idl)=1;
    for t=0:T,
        it=t+1;
        if t==0,
            zeta=zeros(nz,1);
        else
            zeta=z2z*zeta;
            if t>d,
                ImpR.F(:,idl,it,nh+1)=zeta;
                continue;
            end;
            zeta(1:nx)=zeta(1:nx)+R(:,:,t)*Delta;
            zeta((it-1)*nx+1:it*nx)=zeta((it-1)*nx+1:it*nx)+d2x*Delta;
            ImpR.F(:,idl,it,nh+1)=zeta;
        end;
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ih=0;
for id=1:d,
    for iy=1:ny,
        for iu=1:nu,
            ih=ih+1;
            H=zeros(nu,ny,d);
            if id<=d,
                H(iu,iy,id)=1;
            end;
            basec=nx;
            for ig=1:d,
                z2z(1:nx,basec+1:basec+nx)=u2x*H(:,:,ig)*x2y;
                basec=basec+nx;
            end;
            %%%
            R=zeros(nx,nd,d);
            R(:,:,1)=u2x*H(:,:,1)*d2y+d2x;
            for ig=2:d,
                R(:,:,ig)=u2x*H(:,:,ig)*d2y;
            end;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for idl=1:nd,
                Delta=zeros(nd,1);
                Delta(idl)=1;
                for t=0:T,
                    it=t+1;
                    if t==0,
                        zeta=zeros(nz,1);
                    else      
                        zeta=z2z*zeta;
                        if t>d,
                            ImpR.F(:,idl,it,ih)=zeta;
                            continue;
                        end;
                        zeta(1:nx)=zeta(1:nx)+R(:,:,t)*Delta;
                        zeta((it-1)*nx+1:it*nx)=zeta((it-1)*nx+1:it*nx)+d2x*Delta;
                        ImpR.F(:,idl,it,ih)=zeta;
                    end;
                end;
                ImpR.F(:,idl,:,ih)= ImpR.F(:,idl,:,ih)-ImpR.F(:,idl,:,nh+1);
            end;        
        end;
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ImpR.x=ImpR.F(1:nx,:,:,:);
ImpR.y=zeros(ny,nd,T,nh+1);
ImpR.by=zeros(ny,nd,T,nh+1);
for t=0:T-1,
    it=t+1;
    if t==0,
        ImpR.y(:,:,it,nh+1)=d2y;
        ImpR.by(:,:,it,nh+1)=0;
    end;
    ImpR.y(:,:,it,nh+1)=ImpR.y(:,:,it,nh+1)+x2y*ImpR.x(:,:,it,nh+1);
    for ih=1:nh,
        ImpR.y(:,:,it,ih)=x2y*ImpR.x(:,:,it,ih);
        ImpR.by(:,:,it,ih)= ImpR.y(:,:,it,ih)-x2y*ImpR.F(nx+1:2*nx,:,it,ih);
    end;
end;
ImpR.u=zeros(nu,nd,T,nh+1);
ih=0;
for tau=0:d-1,
    for iy=1:ny,
        for iu=1:nu,
            ih=ih+1;
            for t=0:T-1,
                it=t+1;
                if t<tau,
                    continue;
                end;
                ImpR.u(iu,:,it,ih)=ImpR.y(iy,:,it-tau,nh+1)+ImpR.y(iy,:,it-tau,ih)-ImpR.by(iy,:,it-tau,ih);
            end;
        end;
    end;
end;

        
        
        