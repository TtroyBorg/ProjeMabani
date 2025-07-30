function [x_min,f_min,iteration] = SD_GSS(f,gradiantf,x0,stop_tol,varargin)
iteration=0;
%a0=0;
%b0=10;%for rosenbrock please set b0<2.3
xk=x0;
xk2=zeros(1,length(xk));
if isempty(varargin)
while(norm(xk-xk2)>stop_tol)
gk=-gradiantf(xk);
%objfun=@(a) (f(xk+a*gk));
%[a_min,~] = GSS_96105275(objfun,a0,b0,GSS_tol);
alpha =bracketing(f,gradiantf,gk,xk);
xk2=xk;
xk=xk+alpha*gk;
iteration=iteration+1;
if iteration >1000 
    x_min=xk;
f_min=f(xk);
break
end
end
x_min=xk;
f_min=f(xk);

else
  var=varargin;
  var=cell2mat(var);
 while(norm(xk-xk2)>stop_tol)

gk=-gradiantf(xk,var);
%objfun=@(a) (f(xk+a*gk));
%[a_min,~] = GSS_96105275(objfun,a0,b0,GSS_tol);
alpha =bracketing(f,gradiantf,gk,xk,var);
xk2=xk;
xk=xk+alpha*gk;
iteration=iteration+1;
if iteration >1000 
    x_min=xk;
f_min=f(xk,var);
break
end
end
x_min=xk;
f_min=f(xk,var);
end
end


