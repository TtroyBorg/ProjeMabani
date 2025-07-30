function [x_min,f_min,iteration] = newton_GSS(fu,gf,hf,x0,stop_tol,varargin)
iteration=0;
x_k=x0;
x_k2=zeros(1,length(x0));
if isempty(varargin)
%X=zeros(100,6);

while(norm(x_k2-x_k)>stop_tol)
    mu=1;
    hk=hf(x_k);
    [R,flag] = chol(hk);
while(flag~=0)
    hk=hk+mu*eye(length(hk));
    [R,flag] = chol(hk);
    mu=mu+1;
end
pk=-R\(R'\gf(x_k)');
%pk=-hk\gf(x_k)';
alpha =bracketing(fu,gf,pk',x_k);
%objfun=@(a) (fu(x_k+a*pk'));
%lpha =fminunc(objfun,a0);
%alpha = Sectioning(fu,gf,a0,b0,pk',x_k);
x_k2=x_k;
x_k=x_k+alpha*pk';
%X(iteration,1:4)=x_k;
%X(iteration,5)=alpha;
%X(iteration,6)=fu(x_k);
iteration=iteration+1;
if iteration >1000 
    x_min=x_k;
f_min=fu(x_k);
break
end
end
x_min=x_k;
f_min=fu(x_k);
else 
    var=varargin;
    var=cell2mat(var);
    while(norm(x_k2-x_k)>stop_tol)
    mu=1;
    hk=hf(x_k,var);
    [R,flag] = chol(hk);
while(flag~=0)
    hk=hk+mu*eye(length(hk));
    [R,flag] = chol(hk);
   mu=mu+10;
end
pk=-R\(R'\gf(x_k,var)');
alpha = Sectioning(fu,gf,a0,b0,pk,x_k);
line=@(t) (x_k+t*pk');
objfun=@(t,var) (fu(line(t),var));
[a_min,~] = GSS_96105275(objfun,a0,b0,GSS_tol,var);
x_k2=x_k;
x_k=x_k+a_min*pk';
iteration=iteration+1;
end
x_min=x_k;
f_min=fu(x_k,var);
end
end

