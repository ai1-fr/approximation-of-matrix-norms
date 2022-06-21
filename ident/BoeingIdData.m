function data=BoeingIdData(r_eps)
% Generating data for boeing system identification example 
% u(t+1)=X*[u(t);r(t)], 
%
% Author A. Nemirovski (2022)

N=12; 
data.N=N;
data.X=[0.9957,0.0339,-0.0211,-0.3214,0.0140,0.9886,0.0043,-0.0337;...
    0.0076,0.4699,4.6604,0.0022,-3.4373,1.6648,-0.0079,0.5285;...
    0.0168,-0.0605,0.4038,-0.0029,-0.8219,0.4378,-0.0167,0.0600;...
    0.0091,-0.0370,0.7194,0.9990,-0.4735,0.2491,-0.0091,0.0370];
[nu,dum]=size(data.X);
data.nu=nu;
data.nr=dum-nu;
nr=data.nr;

%
data.ij2ellX=zeros(nu,nu+nr);
data.ellX2ij=zeros(nu*(nu+nr),2);
ell=0;
for j=1:nu+nr
    for i=1:nu
        ell=ell+1;
        data.ij2Xell(i,j)=ell;
        data.Xell2ij(ell,:)=[i,j];
    end
end

%
nx=nu*(nu+nr);
data.nx=nx;
data.P=[];
data.p=[];
data.B=speye(nx);

%
data.md='r';
data.vl=r_eps;
S=(N+1)*nu+N*nr;
data.S=S;
data.ell2ti=zeros(S,2);
data.ti2ell=zeros(N+1,nu+nr);
ell=0;
for t=1:N
    for i=1:nu+nr
        ell=ell+1;
        data.ell2ti(ell,:)=[t,i];
        data.ti2ell(t,i)=ell;
    end
end
for i=1:nu
    ell=ell+1;
    data.ti2ell(N+1,i)=ell;
    data.ell2ti(ell,:)=[N+1,i];
end
data.r2ti=zeros(N*nu,2);
data.ti2r=zeros(N,nu);
r=0;
for t=1:N
    for i=1:nu
        r=r+1;
        data.r2ti(r,:)=[t,i];
        data.ti2r(t,i)=r;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating observations
data.urs=randn(nu+nr,N+1);
data.eps=zeros(S,1);
for t=2:N+1
    data.urs(:,t)=[data.X*data.urs(:,t-1);randn(nr,1)]; 
end
data.obs=zeros(S,1);
data.xis=zeros(S,1);
for ell=1:S
    t=data.ell2ti(ell,1);
    i=data.ell2ti(ell,2);
    tru=data.urs(i,t);
    if rand(1,1)>0.3
        err=data.vl;
    else
        err=-data.vl;
    end
       data.obs(ell)=tru-err*max(1,abs(tru));
        data.xis(ell)=err*max(1,abs(tru));
        data.eps(ell)=data.vl*max(1,abs(data.obs(ell)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.nx=nu*(nu+nr);
data.nq=N*nu;
nq=data.nq;
ii=zeros(N*nu*(nu+nr),1);
jj=zeros(N*nu*(nu+nr),1);
vv=zeros(N*nu*(nu+nr),1);
iis=zeros(nu,nu+nr,N);
jjs=zeros(nu,nu+nr,N);
vvs=zeros(nu,nu+nr,N);
inds=zeros(nu+nr,N);
data.q=zeros(nq,1);
data.qs=zeros(nq,S);
data.Qs=cell(S,1);
ind=0;
for t=1:N
    for i=1:nu
        ir=data.ti2r(t,i);
        ell=data.ti2ell(t+1,i);
        data.q(ir)=data.obs(ell);
        data.qs(ir,ell)=-1;
        for j=1:nu+nr
            ic=data.ij2Xell(i,j);
            ell=data.ti2ell(t,j);
            vl=data.obs(ell);
            ind=ind+1;
            ii(ind)=ir;
            jj(ind)=ic;
            vv(ind)=vl;
            inds(j,t)=inds(j,t)+1;
            is=inds(j,t);
            iis(is,j,t)=ir;
            jjs(is,j,t)=ic;
            vvs(is,j,t)=-1;
        end
    end
end
data.Q=sparse(ii,jj,vv,nq,nx);
for j=1:nu+nr
    for t=1:N
        k=inds(j,t);
        s=data.ti2ell(t,j);
        data.Qs{s}=sparse(iis(1:k,j,t),jjs(1:k,j,t),vvs(1:k,j,t),nq,nx);
    end
end
for j=1:nu
    data.Qs{S-j+1}=sparse([],[],[],nq,nx);
end
    
%
% if ~isempty(data.P)
%     Tmp=(data.P*data.P')^(-1)*data.P;
%     data.Pi=eye(nx)-data.P'*Tmp;
%     data.barx=Tmp'*data.p;
% else
% Defining matrix and r.h.s. of the system of linear constraints
    data.Pi=eye(nx);
    data.barx=zeros(nx,1);
% end
data.bary=data.B*data.barx;
data.barq=data.q-data.Q*data.barx;
data.barqs=data.qs;
data.barQs=cell(S,1);
for ell=1:S
    data.barQs{ell}=data.eps(ell)*data.Qs{ell};
    data.barqs(:,ell)=data.barQs{ell}*data.barx-data.eps(ell)*data.qs(:,ell);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.x=zeros(nx,1);
for i=1:nu
    for j=1:nu+nr
        Xell=data.ij2Xell(i,j);
        data.x(Xell)=data.X(i,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    


            
        
            
        
    

