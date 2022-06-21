function res=Upsilons(data,UpsBnd,AUTHOR)
% Upsilons computes linear estimate
% cf. Section 4.3.2 of "Tight computationally efficient approximation of matrix
% norms with applications" by A. Juditsky, G. Kotsalis and A. Nemirovski
%
% Syntax:
%   res=Upsilons(data,UpsBnd,verbose)
% Inputs:
%   data: data structure containing training data
%   see BoeingIdData.m for example of created data structure
%   UpsBnd: bound on the Ypsilon parameter of matrix
%
% Output: 
%   res: structure with fields
%       E: E matrix of the estimate
%       obj: optimal value of the optimization problem
%       Upsi: computed value of Ypsilon
%
% Author: A. Nemirovski (2022)

nx=data.nx;
Pi=data.Pi;
n=nx;
m=data.nq;
S=data.S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Put=cell(S,1);
SuppC=cell(S,1);
nC=zeros(S,1);
for s=1:S,
    tmp=ones(1,m)*abs(data.barQs{s}*Pi);
    [ii,jj,vv]=find(tmp);
    nc=nnz(jj);
    nC(s)=nc;
    Put{s}=sparse([1:1:nc],jj,ones(nc,1),nc,n);
    SuppC{s}=zeros(m,nc);
    Tmp=data.barQs{s}*Pi;
    it=0;
    for i=1:n,
        if norm(Tmp(:,i))>0,
            it=it+1;
            SuppC{s}(:,it)=Tmp(:,i);
        end;
    end;
end;
mxnC=max(nC);

%
nvar=1;
G2V=zeros(n,n);
Hs=zeros(mxnC,mxnC,S);
nM=sum(nC)+n;
M2V=ones(nM,nM);
for i=1:n,
    for j=1:i,
        nvar=nvar+1;
        G2V(i,j)=nvar;
        G2V(j,i)=nvar;
        M2V(i,j)=nvar;
        M2V(j,i)=nvar;
    end
end
base=n;
for s=1:S
    for i=1:nC(s)
        for j=1:i
            nvar=nvar+1;
            Hs(i,j,s)=nvar;
            Hs(j,i,s)=nvar;
            M2V(base+i,base+j)=nvar;
            M2V(base+j,base+i)=nvar;
        end
    end
    base=base+nC(s);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cvx_begin
variable Vars(nvar,1);
Vars(1) == 0;
variable x(n,1);
variable E(m,n);
variable upsilon;
variable wupsilon;
wupsilon >= 0;
variable mus(S,1);
variable lm;
variable Upsi;
variable bG(n,n) symmetric;
variable bH(n,n) symmetric;
expression Z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x == Pi*E'*data.barq;
%%%%%%%%%%%%%%%%%%%%%%%%
[upsilon*eye(n),Pi*E'*data.barqs;
    data.barqs'*E*Pi,diag(mus)] == semidefinite(n+S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wupsilon*eye(n)-bG-Vars(G2V) == semidefinite(n);
Z=bH;
for s=1:S
    p=nC(s);
    Z=Z+Put{s}'*Vars(Hs(1:p,1:p,s))*Put{s};
end
lm*eye(n) - Z == semidefinite(n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[bG,Pi*(E'*data.Q-speye(n))*Pi;Pi*(data.Q'*E-speye(n))*Pi,bH] == semidefinite(2*n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z=Vars(M2V);
base=n;
for s=1:S
    nc=nC(s);
    Z(1:n,base+1:base+nc) = Z(1:n,base+1:base+nc)+Pi*E'*SuppC{s};
    Z(base+1:base+nc,1:n) = Z(base+1:base+nc,1:n)+SuppC{s}'*E*Pi;
    base=base+nc;
end
Z == semidefinite(nM);
0.5*(wupsilon+lm) <= Upsi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if UpsBnd==inf
    minimize Upsi;
else
    Upsi <= UpsBnd;
    minimize norm(x)+0.5*(upsilon+sum(mus));
end;
cvx_end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res.status=cvx_status;
res.UpsBnd=UpsBnd;
res.E=E;
res.obj=cvx_optval;
res.Ups0=0.5*(upsilon+sum(mus));
res.Ups=0.5*(wupsilon+lm);
res.Upsi=Upsi;
res.x=Pi*E'*data.barq+data.barx;
res.PETbarq=norm(x);
if res.Upsi<1,
    res.riskUB=(Upsi/(1-Upsi))*(res.PETbarq+res.Ups0)+res.Ups0;
else
    res.riskUB=inf;
end
if res.Upsi<1
    res.Gamma=(norm(x)+0.5*(upsilon+sum(mus)))/(1-res.Upsi);
else
    res.Gamma=inf;
end
if AUTHOR
    res.x=x+data.barx;
    res.bG=bG;
    res.bH=bH;
    res.G=Vars(G2V);
    res.Hs=cell(S,1);
    res.nC=nC;
    for s=1:S,
        p=nC(s);
        res.Hs{s}=Vars(Hs(1:p,1:p,s));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    zero.x=norm(x-Pi*E'*data.barq);
    Tmp=full([upsilon*eye(n),Pi*E'*data.barqs;data.barqs'*E,diag(mus)]);
    pos=min(eig(Tmp));
    zero.ups=max(-pos,0);
    Tmp=full([bG,Pi*(E'*data.Q-eye(n))*Pi;Pi*(data.Q'*E-eye(n))*Pi,bH]);
    pos=min(eig(Tmp));
    zero.bG=max(-pos,0);
    totG=bG;
    totH=bH;
    zero.GH=0;
    for s=1:S
        p=nC(s);
        Hc=Put{s}'*Vars(Hs(1:p,1:p,s))*Put{s};
        Up=Pi*E'*data.barQs{s}*Pi;
        Gc=Up*(1.e-12*eye(n)+Hc)^(-1)*Up';
        Tmp=full([Gc,Up;Up',Hc]);
        pos=min(eig(Tmp));
        zero.GH=max(zero.GH,-pos);
        totG=totG+Gc;
        totH=totH+Hc;res.PETbarq=norm(x);
    end
    pos=min(eig(wupsilon*eye(n)-totG));
    zero.totG=max(0,-pos);
    pos=min(eig(lm*eye(n)-totH));
    zero.totH=max(0,-pos);
    pos=res.Upsi-0.5*(wupsilon+lm);
    zero.Ups=max(-pos,0);
    res.zero=zero;
    risk=norm(res.x-data.x);
    fprintf('***** UpsBnd=%5.1f',UpsBnd)
    fprintf('      Ups0=%5.4f Ups=%5.4f Add=%5.4f RiskUB=%5.4f risk=%5.4f',...
        res.Ups0,res.Ups,res.PETbarq,res.riskUB,risk)
end




