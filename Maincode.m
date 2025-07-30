stop_tol=10^(-3);
%x0=[1,2];%rosenberg
x0=[1,2,2,2];%powell
GSS_tol=10^(-5);
global i_f;
global i_g;
%global i_h;
i_f=0;
i_g=0;
%i_h=0;
%[x_min,f_min,iteration] = newton_GSS(@f,@gradiantf,@hessianf,x0,stop_tol);
%[x_min,f_min,iteration] = SD_GSS(@f,@gradiantf,x0,stop_tol);
[x_min,f_min,iteration] = my_BFGS(@f,@gradiantf,x0,stop_tol);