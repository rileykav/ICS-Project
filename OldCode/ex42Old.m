1;









T=2;
Nt=1001;
delta_t=T/(Nt-1);
delta_x=0.01;
d=0.001;
gamma_const=5;
Nx=int32(1+(1/delta_x));
t=linspace(0,T,Nt);
x=linspace(0,1,Nx);




a=-2;
b=1;
c=1;
D = diag(a*ones(1,Nx)) + diag(b*ones(1,Nx-1),1) + diag(c*ones(1,Nx-1),-1);
D(1,:)=0;
D(Nx,:)=0;
D=D/(delta_x^2);
der_a=-2;
der_b=1;
der_c=1;
der_d=0;
der_f=0;
D2 = diag(der_a*ones(1,Nx)) + diag(der_b*ones(1,Nx-1),1) +diag(der_f*ones(1,Nx-2),2)+diag(der_d*ones(1,Nx-2),-2)+ diag(der_c*ones(1,Nx-1),-1);
D2(1,Nx)=der_c;
D2(Nx,1)=der_b;
D2=D2 / delta_x^2;

u_eulfor=zeros(Nx,Nt);
for i=1:Nx
    u_eulfor(i,1)=exp(-50*(x(i)-.5)^2);
end
for i=2:Nt
    g = u_eulfor(:,i-1) + delta_t*gamma_const*(u_eulfor(:,i-1).*(1 - u_eulfor(:,i-1)));
    u_eulfor(:,i) = (eye(Nx) - delta_t*d*D2)\g;
end

xi=3;
tau=4;




listoftimes=1:(Nt-1)/20:Nt
clf
hold on
for i=1:20
    tau=listoftimes(i)
    plot(x,u_eulfor(:,tau));
end
title("\\Delta x=0.01")
h=get(gcf, "currentaxes");
set(h, "fontsize", 12, "linewidth", 1);
