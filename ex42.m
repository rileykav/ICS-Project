1;









T=1;
Nt=10;
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
u_eulfor=zeros(Nx,Nt);
for i=1:Nx
    u_eulfor(i,1)=exp(-50*(x(i)-.5)^2);
end
for i=2:Nt
    g = u_eulfor(:,i-1) + delta_t*gamma_const*(u_eulfor(:,i-1).*(1 - u_eulfor(:,i-1)));
    u_eulfor(:,i) = (eye(Nx) - delta_t*d*D)\g;
end

xi
tau=10





clf
plot(x,u_eulfor(:,tau));
title("\\Delta x=0.01")
h=get(gcf, "currentaxes");
set(h, "fontsize", 12, "linewidth", 1);
