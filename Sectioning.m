function alphaj = Sectioning(fu,gf,a,b,pk,x_k)
c1=10^(-3);
c2=0.9;
%tau1=0.1;
%tau2=0.5;
epsilon=10^(-10);
a_istrue=1;
line=@(t) (x_k+t*pk);
fphi=@(t) (fu(line(t)));
Dphi=@(t)(gf(line(t))*pk');
candid_a=0;
candid_b=0;
fphi_b=0;
fphi_a=0;

while(1)
    Dphi_a=Dphi(a);
    fphi_0=fphi(0);
    Dphi_0=Dphi(0);
    a_istrue=1;
    candid_a=a;
    candid_b=b;
    phi_1=fphi(a);
    phi_2=fphi(b);
    Dphi_1=(b-a)*Dphi_a;
    x1=a;
    x2=b;
    %x1=a+tau1*(b-a);
    %x2=b-tau2*(b-a);
    x3=(-Dphi_1)/(2*(phi_2-phi_1-Dphi_1));
    
   if x3>x2 || x3<x1
   x3=x2; 
    end
    xt=[x1;x2;x3];
    q=@(z) (phi_1+Dphi_1*z+(phi_2-phi_1-Dphi_1)*(z.^2));
    [~,I] = min([q(x1),q(x2),q(x3)]); 
    alpha=xt(I);
    if alpha==x3
        alpha=a+x3*(b-a);
    end
   %alpha=(a+b)/2;
    if (a-b)*Dphi_a<epsilon
     alphaj=a; 
      break
    end
   fphi_alpha= fphi(alpha);
   Dphi_alpha=Dphi(alpha);
   
if(fphi_alpha>fphi_0+c1*alpha*Dphi_a|| fphi_alpha>=phi_1)
    candid_b=alpha;
    if abs(candid_b-candid_a)<0.5*abs(b-a)
        a_istrue=0;
        b=candid_b;
        a=candid_a;
    end
else
    if abs(Dphi_alpha)<-c2*Dphi_a
        alphaj=alpha;
        break
    end
 
    if (candid_b-candid_a)*Dphi_alpha>=0
        candid_b=candid_a;
    end
    candid_a=alpha;
    if abs(candid_b-candid_a)<0.5*abs(b-a)
        a_istrue=0;
        b=candid_b;
        a=candid_a;
    end
end
    if (a_istrue)
        alpha=(a+b)/2;
        fphi_alpha= fphi(alpha);
        Dphi_alpha=Dphi(alpha);
        if(fphi_alpha>fphi_0+c1*alpha*Dphi_0|| fphi_alpha>=phi_1)
    b=alpha;
        else
    if abs(Dphi_alpha)<-c2*Dphi_0
        alphaj=alpha;
        break
    end
 
    if (b-a)*Dphi_alpha>=0
        b=a;
    end
    a=alpha;
        
        end
    
    
    
    end
end
end

