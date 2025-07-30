function [x_min,f_min,iteration] = my_BFGS(fu,gf,x0,stop_tol,varargin)
iteration=0;
xk=x0;
xk2=zeros(1,length(x0));
gk=gf(xk);
n=length(gk);
ck=eye(n,n);
if isempty(varargin)
%X=zeros(100,6);
 pk=-ck*gk';
    alpha =bracketing(fu,gf,pk',xk);
    xk2=xk;
    xk=xk2+alpha*pk';
    gk2=gk;
    gk=gf(xk);
    dk=xk-xk2;
    qk=gk-gk2;
    ck=((dk*qk')/(qk*qk'))*eye(n,n);
    %gk=gf(x0);
    %xk=x0;

while(norm(gk)>stop_tol)
    
pk=-ck*gk';
    alpha =bracketing(fu,gf,pk',xk);
xk2=xk;
    xk=xk2+alpha*pk';
    gk2=gk;
    gk=gf(xk);
    dk=xk-xk2;
    qk=gk-gk2;
a=1/(qk*dk');
    ck=(eye(n,n)-a*dk'*qk)*ck*(eye(n,n)-a*qk'*dk)+a*(dk')*dk;

%X(iteration,1:4)=x_k;
%X(iteration,5)=alpha;
%X(iteration,6)=fu(x_k);
iteration=iteration+1;
if iteration >100
    x_min=xk;
f_min=fu(xk);
break
end
end
x_min=xk;
f_min=fu(xk);
else 
    var=varargin;
    var=cell2mat(var);
    while(norm(xk2-xk)>stop_tol)
    mu=1;
    hk=hf(xk,var);
    [R,flag] = chol(hk);
while(flag~=0)
    hk=hk+mu*eye(length(hk));
    [R,flag] = chol(hk);
   mu=mu+10;
end
pk=-R\(R'\gf(xk,var)');
alpha = Sectioning(fu,gf,a0,b0,pk,xk);
line=@(t) (xk+t*pk');
objfun=@(t,var) (fu(line(t),var));
[a_min,~] = GSS_96105275(objfun,a0,b0,GSS_tol,var);
x_k2=xk;
xk=xk+a_min*pk';
iteration=iteration+1;
end
x_min=xk;
f_min=fu(xk,var);
end
end
