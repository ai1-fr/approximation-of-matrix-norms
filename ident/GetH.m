function res=GetH(data,coeff)
% GetH computes linear estimate by solving optimization problem in (39b)
% cf. Section 4.3.2 of "Tight computationally efficient approximation of matrix
% norms with applications" by A. Juditsky, G. Kotsalis and A. Nemirovski
%
% Syntax:
%   res=GetH(data,coeff)
% Inputs:
%   data: data structure containing training data
%       see BoeingIdData.m for example of created data structure
%   coeff: value of the selected Gamma parameter
%
% Output: 
%   res: structure with fields
%       E: E matrix of the estimate
%       obj: optimal value of the optimization problem
%       Upsi: computed value of Ypsilon
%       H:  computed H estimation matrix 
%
% Author: A. Nemirovski (2022)

nx=data.nx;
Pi=data.Pi;
n=nx;
m=data.nq;
S=data.S;    
nu=size(data.B,1);
B=data.B;
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
        if norm(Tmp(i,:))>0,
            it=it+1;
            SuppC{s}(:,it)=Tmp(:,i);
        end;
    end;
end;
mxnC=max(nC);
%%%
nvar=1;
G2V=zeros(nu,nu);
Hs=zeros(mxnC,mxnC,S);
nM=sum(nC)+nu;
M2V=ones(nM,nM);
for i=1:nu,
    for j=1:i,
        nvar=nvar+1;
        G2V(i,j)=nvar;
        G2V(j,i)=nvar;
        M2V(i,j)=nvar;
        M2V(j,i)=nvar;
    end;
end;
base=nu;
for s=1:S,
    for i=1:nC(s),
        for j=1:i,
            nvar=nvar+1;
            Hs(i,j,s)=nvar;
            Hs(j,i,s)=nvar;
            M2V(base+i,base+j)=nvar;
            M2V(base+j,base+i)=nvar;
        end;
    end;
    base=base+nC(s);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cvx_begin;
variable Vars(nvar,1);
Vars(1) == 0;
variable H(m,nu);
variable upsilon;
variable wupsilon;
wupsilon >= 0;
variable mus(S,1);
variable lm;
variable Upsi;
variable bG(nu,nu) symmetric;
variable bH(n,n) symmetric;
expression Z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[upsilon*eye(nu),H'*data.barqs;
    data.barqs'*H,diag(mus)] == semidefinite(nu+S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wupsilon*eye(nu)-bG-Vars(G2V) == semidefinite(nu);
Z=bH;
for s=1:S,
    p=nC(s);
    Z=Z+Put{s}'*Vars(Hs(1:p,1:p,s))*Put{s};
end;
lm*eye(n) - Z == semidefinite(n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[bG,(H'*data.Q-B)*Pi;Pi*(data.Q'*H-B'),bH] == semidefinite(nu+n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z=Vars(M2V);
base=nu;
for s=1:S,
    nc=nC(s);
    Z(1:nu,base+1:base+nc) = Z(1:nu,base+1:base+nc)+H'*SuppC{s};
    Z(base+1:base+nc,1:nu) = Z(base+1:base+nc,1:nu)+SuppC{s}'*H;
    base=base+nc;
end;
Z == semidefinite(nM);
0.5*(wupsilon+lm) <= Upsi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minimize coeff*Upsi+0.5*(upsilon+sum(mus));;
cvx_end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res.status=cvx_status;
res.Upsi=Upsi;
res.H=H;
res.bG=bG;
res.bH=bH;
res.G=Vars(G2V);
res.Hs=cell(S,1);
res.nC=nC;
for s=1:S,
    p=nC(s);
    res.Hs{s}=Vars(Hs(1:p,1:p,s));
end;
res.obj=cvx_optval;
res.BXrec=data.bary+H'*data.barq;

    
    
    
    