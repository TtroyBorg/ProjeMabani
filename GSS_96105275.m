function [x_min,N] = GSS_96105275(f,a0,b0,epsilon,varargin)

%h=varargin{1};
rho=(3-sqrt(5))/2;
N=ceil((log10(epsilon/(b0-a0)))/log10(1-rho));
a1=a0+rho*(b0-a0);
b1=a0+(1-rho)*(b0-a0);
y_help1=f(b1);
y_help2=f(a1);

if isempty(varargin)
    for i=1:N
         if y_help2<y_help1 %y_help2<y_help1;f(a1)<f(b1)y_help2<y_help1
         b0=b1;
         b1=a1;
        a1=a0+rho*(b0-a0);
        y_help1=y_help2;
        y_help2=f(a1);
    else
         a0=a1;
         a1=b1;
        b1=a0+(1-rho)*(b0-a0);
        y_help2=y_help1;
        y_help1=f(b1);
         end
    end
    x_min=(a0+b0)/2;
else 
    var=varargin;
    var=cell2mat(var);
    for i=1:N
  if  f(a1,var)<f(b1,var)
      b0=b1;
      b1=a1;
      a1=a0+rho*(b0-a0);
  else
      a0=a1;
      a1=b1;
      b1=a0+(1-rho)*(b0-a0);
  end
    end
x_min=(a0+b0)/2;
end
end
