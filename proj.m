1;


delta_x=0.01;
T=0.03;
Nt=100;
d=0.001;
gamma_const=5;
delta_t=T/(Nt-1);
Nx=int32(1+(1/delta_x));
t=linspace(0,T,Nt);
x=linspace(0,1,Nx);




der_a=0;
der_b=-1;
der_c=1;
der_d=-1/2;
der_f=1/2;
D3 = diag(der_a*ones(1,Nx)) + diag(der_b*ones(1,Nx-1),1) +diag(der_f*ones(1,Nx-2),2)+diag(der_d*ones(1,Nx-2),-2)+ diag(der_c*ones(1,Nx-1),-1);
D3(1,Nx)=der_c;
D3(1,Nx-1)=der_d;
D3(2,Nx)=der_d;
D3(Nx,1)=der_b;
D3(Nx,2)=der_f;
D3(Nx-1,1)=der_f;

der_a=-2;
der_b=1;
der_c=1;
der_d=0;
der_f=0;
D2 = diag(der_a*ones(1,Nx)) + diag(der_b*ones(1,Nx-1),1) +diag(der_f*ones(1,Nx-2),2)+diag(der_d*ones(1,Nx-2),-2)+ diag(der_c*ones(1,Nx-1),-1);
D2(1,Nx)=der_c;
D2(Nx,1)=der_b;
D2=D2 / delta_x^2;

D2minuscubed=(-D2)^3;
D2minuscubedover4minus=sqrtm(sqrtm(D2minuscubed));





% D3=D3/(delta_x^2);
D3;
D3over2=sqrtm(sqrtm(D3^2 / delta_x^3));
D3minusover2=sqrtm(-D3 / delta_x^3);
D3over2=(sqrtm(D3 / delta_x^3));
% D3over2=(D3over2 + D3over2minus)/(sqrt(2))






u_initial=zeros(Nx,1);
for i=1:Nx
    u_initial(i)=exp(-100*(x(i)-0.5)^2);
end

u_exmp61=zeros(Nx,Nt);
for i=1:Nt
    u_exmp61(:,i)=(expm(t(i).*D3over2)*u_initial);
end
u_exmpproj=zeros(Nx,Nt);
for i=1:Nt
    u_exmpproj(:,i)=(expm(-t(i).*D2minuscubedover4minus)*u_initial);
end
% 
% u_eulfor=zeros(Nx,Nt);
% u_eulfor(:,1)=u_initial;
% for i=1:Nt
%     u_eulfor(:,i+1)=(eye(Nx)+delta_t.*D3over2)*u_eulfor(:,i);
% end
% 
% u_eulback=zeros(Nx,Nt);
% u_eulback(:,1)=u_initial;
% for i=1:Nt
%     A=eye(Nx)-delta_t.*D3over2;
%     u_eulback(:,i+1)=A\u_eulfor(:,i);
% end

% for i=2:Nt
%     g = u_eulfor(:,i-1) + delta_t*gamma_const*(u_eulfor(:,i-1).*(1 - u_eulfor(:,i-1)));
%     A = -(d/sqrt(2))*(D3over2 + D3minusover2);
%     u_eulfor(:,i) = (eye(Nx) - delta_t*d*D3over2)\g;
% end













clf 
hold on
labels=[''];
listoftimes=[30];
% listoftimes=[10,20,30,40,50,60];
% listoftimes=[3]
for ti=1:size(listoftimes)(2)
    tau=listoftimes(ti);
%     plot(x, u_exmp(:,tau))
    plot(x, u_exmp61(:,tau))
    plot(x, u_exmpproj(:,tau))
%     plot(x, u_eulback(:,tau))
%     labels=[labels;sprintf("time(%d)=%d",tau,t(tau))];
end
% legend(labels=[sprintf()])
title(sprintf('\\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt))
h=get(gcf, "currentaxes");






Nx=10
der_a=-2;
der_b=1;
der_c=1;
der_d=0;
der_f=0;
D2 = diag(der_a*ones(1,Nx)) + diag(der_b*ones(1,Nx-1),1) +diag(der_f*ones(1,Nx-2),2)+diag(der_d*ones(1,Nx-2),-2)+ diag(der_c*ones(1,Nx-1),-1);
D2(1,Nx)=der_c;
D2(Nx,1)=der_b;

D2

der_a=0;
der_b=1;
der_c=-1;
der_d=0;
der_f=0;
D1C = diag(der_a*ones(1,Nx)) + diag(der_b*ones(1,Nx-1),1) +diag(der_f*ones(1,Nx-2),2)+diag(der_d*ones(1,Nx-2),-2)+ diag(der_c*ones(1,Nx-1),-1);
D1C(1,Nx)=der_c;
D1C(Nx,1)=der_b;

D1C=D1C/2

der_a=-2;
der_b=2;
der_c=0;
der_d=0;
der_f=0;
D1plus = diag(der_a*ones(1,Nx)) + diag(der_b*ones(1,Nx-1),1) +diag(der_f*ones(1,Nx-2),2)+diag(der_d*ones(1,Nx-2),-2)+ diag(der_c*ones(1,Nx-1),-1);
D1plus(1,Nx)=der_c;
D1plus(Nx,1)=der_b;

D1plus=D1plus/2

der_a=2;
der_b=0;
der_c=-2;
der_d=0;
der_f=0;
D1minus = diag(der_a*ones(1,Nx)) + diag(der_b*ones(1,Nx-1),1) +diag(der_f*ones(1,Nx-2),2)+diag(der_d*ones(1,Nx-2),-2)+ diag(der_c*ones(1,Nx-1),-1);
D1minus(1,Nx)=der_c;
D1minus(Nx,1)=der_b;

D1minus=D1minus/2


(D1plus+D1minus)/2
D1plus*D1minus
