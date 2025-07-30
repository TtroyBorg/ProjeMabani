function alpha = bracketing(fu,gf,pk,x_k,varargin)
f_bar=10^(-9);
c1=10^(-3);
c2=0.9;
line=@(t) (x_k+t*pk);
fphi=@(t) (fu(line(t)));
Dphi=@(t)(gf(line(t))*pk');
amax=(f_bar-fphi(0))/(c1*Dphi(0));
a=0;
b=1;
while(1)
    fphi_b=fphi(b);
    fphi_a=fphi(a);
    fphi_0=fphi(0);
    Dphi_0=Dphi(0);
    Dphi_b=Dphi(b);
   if fphi_b<=f_bar
       alpha=b;
       break
   end
   if fphi_b>fphi_0+c1*b*Dphi_0 || fphi_b>=fphi_a
       alpha = Sectioning(fu,gf,a,b,pk,x_k);
    break
   end
   if abs(Dphi_b)<= -c2*Dphi_0
       alpha=b;
       break
   end
   if (Dphi_b)>=0
       alpha = Sectioning(fu,gf,b,a,pk,x_k);
    break
   end
   c=b;
   b=fchoose(b,a,amax);
   a=c;
end
end

