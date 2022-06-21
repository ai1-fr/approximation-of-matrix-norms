function [data,cntr]=GetBoeing747
% generating Boeing747 demo data
%
% Author: A. Nemirovski (2022)

data.nx=4;
data.ny=2;
data.nu=2;
data.nd=2;
x2xc=[-0.003,0.039,0,-0.322;...
    -0.065,-0.319,7.74,0;...
    0.020,-0.101,-0.429,0;...
    0,0,1,0];
d2xc=[0.003,-0.039;0.065,0.319;-0.020,0.101;0,0];
u2xc=[0.01,1;-0.18,-0.04;-1.16,0.598;0,0];
x2y=[1,0,0,0;0,-1,0,7.74];
dt=1;
Nt=200;
ddt=dt/Nt;
x2x=expm(x2xc);
u2x=zeros(data.nx,data.nu);
d2x=zeros(data.nx,data.nd);
d2y=zeros(data.ny,data.nd);
for it=1:Nt-1
    F=expm((Nt-it+0.5)*ddt*x2xc);
    u2x=u2x+F*u2xc*ddt;
    d2x=d2x+F*d2xc*ddt;
end
data.x2x=x2x;
data.u2x=u2x;
data.d2x=d2x;
data.x2y=x2y;
data.d2y=d2y;


cntr.pd=2;
cntr.px=2;
cntr.py=2;
cntr.pu=2;

cntr.xyu=[1,0,0];
cntr.W.A=[];
cntr.W.b=[];

cntr.d=16;
cntr.T=256;
end % endof GetBoeing747

