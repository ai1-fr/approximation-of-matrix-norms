function ds=GetBadDs(data,cntr,Targ)
% Generating worst-case disturbances 
%
% Author: A. Nemirovski (2022)

p=cntr.pd;
nx=data.nx;
ny=data.ny;
if cntr.xyu(1)>cntr.xyu(2),
    ind='x';
    r=cntr.px;
    nt=nx;
else
    ind='y';
    r=cntr.py;
    nt=ny;
end;
if r==1,
    rs=inf;
else
    rs=r/(r-1);
end;
nd=data.nd;
nx=data.nx;
ny=data.ny;
d2targ=zeros(nt,Targ*nd);
for i=1:nd*Targ,
    e=zeros(nd*Targ,1);
    e(i)=1;
    ds=reshape(e,nd,Targ);
    z=zeros(nx,1);
    for t=1:Targ,
        z=data.x2x*z+data.d2x*ds(:,t);
        y=data.x2y*z;
    end;
    if ind=='x',
        d2targ(:,i)=z;
    else
        d2targ(:,i)=y;
    end;
end;
ds=randn(nd,cntr.T);
for t=1:t,
    ds(:,t)=ds(:,t)/norm(ds(:,t),r);
end;
d2t=randn(nt,1);
d2t=d2t/norm(d2t,rs);
old=0;
while(1==1)
    g=reshape(d2targ'*d2t,nd,Targ);
    for t=1:Targ,
        e=g(:,t);
        if p==inf,
            ds(:,t)=sign(e);
        else
            ds(:,t)=sign(e).*(abs(e).^(1/(p-1)));
        end;
        ds(:,t)=ds(:,t)/norm(ds(:,t),p);
    end;
    img=d2targ*reshape(ds(:,1:Targ),nd*Targ,1);
    new=norm(img,r);
    %disp(sprintf('%5.4f -> %5.4f',old,new));
    if new-old<1.e-4,
        break;
    end;
    if r==1,
        d2t=sign(img);
    else
        d2t=sign(img).*(abs(img).^(r-1));
        d2t=d2t/norm(d2t,rs);
    end;
    old=new;
end;
        
    
        
        