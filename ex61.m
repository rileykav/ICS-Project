1;


delta_x=0.01;
T=0.03;
Nt=3001;
delta_t=T/(Nt-1);
stability=delta_t/(delta_x)^(3/2);
if stability>0.5
    printf("EF Unstable, stability ratio: %.1f\n", stability)
end

Nx=int32(1+(1/delta_x));
t=linspace(0,T,Nt);
x=linspace(0,1,Nx);




a=0;
b=-1;
c=1;
d=-1/2;
f=1/2;
D3 = diag(a*ones(1,Nx)) + diag(b*ones(1,Nx-1),1) +diag(f*ones(1,Nx-2),2)+diag(d*ones(1,Nx-2),-2)+ diag(c*ones(1,Nx-1),-1);
D3(1,Nx)=c;
D3(1,Nx-1)=d;
D3(2,Nx)=d;
D3(Nx,1)=b;
D3(Nx,2)=f;
D3(Nx-1,1)=f;
% D3;
D3=D3/(delta_x^3);
D3over2=-sqrtm(D3);
% D3over2=-sqrtm(-D3);
% D3over2=-(sqrtm(D3) + sqrtm(-D3))/sqrt(2);



u_initial=zeros(Nx,1);
for i=1:Nx
    u_initial(i)=exp(-100*(x(i)-0.5)^2);
end


u_expm=zeros(Nx,Nt);
for i=1:Nt
    u_expm(:,i)=(expm(t(i).*D3over2)*u_initial);
end

u_eulfor=zeros(Nx,Nt);
u_eulfor(:,1)=u_initial;
for i=1:Nt-1
    u_eulfor(:,i+1)=(eye(Nx)+delta_t.*D3over2)*u_eulfor(:,i);
end

u_eulback=zeros(Nx,Nt);
u_eulback(:,1)=u_initial;
for i=1:Nt-1
    A=eye(Nx)-delta_t.*D3over2;
    u_eulback(:,i+1)=A\u_eulfor(:,i);
end

clf 
hold on
labels=[''];
% listoftimes=[3];
listoftimes=1:round((Nt-1)/5):Nt
% listoftimes=[3]
%
%

function Plot3D(plot_x, plot_y, plot_z)
    % Assume plot_x and plot_y are Nx, Ny lenght vectors and
    % plot_z is a Nx*Ny matrix
    [XX,YY]=meshgrid(plot_x,plot_y);
    surf(XX,YY,plot_z)
endfunction

N=50;x_temp=linspace(-1,1,N);y_temp=linspace(-1,1,N);z_temp=zeros(N,N);for i=1:N;for j=1:N;z_temp(i,j)=exp(-(x_temp(i)^2+ y_temp(j)^2));end;end;

clf
hold on 
% Plot3D(x_temp,y_temp,z_temp)
cmap=colormap();
for ti=1:size(listoftimes)(2)
    tau=listoftimes(ti);
%     Plot3D(x,t,real(u_eulback'))
%     xlabel('x')
%     ylabel('t')
    colour_i=round(1+(tau*63/Nt));
    plot(x, u_expm(:,tau),'Color', color=cmap(colour_i,:),'linewidth',2)
%     plot(x, u_eulfor(:,tau))
%     plot(x, u_eulback(:,tau))
    axis([0 1 0 1])
%     labels=[labels;sprintf("time(%d)=%d",tau,t(tau))];
end
xlabel('x')
ylabel('u')
% legend(labels=labels)
h=get(gcf, "currentaxes");
set(h, "fontsize", 12, "linewidth", 1);
cbar=colorbar()
caxis([0,T])
ylabel(cbar,'t','Rotation',0)
title(sprintf('EXPM, D^{3/2}=-(D^3)^{1/2}) \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
filename=sprintf('Ex61/Ex61-EXPM-posDer-multi-t-delx=%d.png',delta_x);
% title(sprintf('EXPM, D^{3/2}=-(-D^3)^{1/2}) \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
% filename=sprintf('Ex61/Ex61-EXPM-minusDer-multi-t-delx=%d.png',delta_x);
% title(sprintf('EXPM, D^{3/2}=-((D^3)^{1/2}+(-D^3)^{1/2})/sqrt(2) \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
% filename=sprintf('Ex61/Ex61-EXPM-mixedDer-multi-t-delx=%d.png',delta_x);
print(1,filename, '-dpng','-S600,400')
% 

clf
hold on 
% Plot3D(x_temp,y_temp,z_temp)
cmap=colormap();
for ti=1:size(listoftimes)(2)
    tau=listoftimes(ti);
%     Plot3D(x,t,real(u_eulback'))
%     xlabel('x')
%     ylabel('t')
    colour_i=round(1+(tau*63/Nt));
    plot(x, u_eulfor(:,tau),'Color', color=cmap(colour_i,:),'linewidth',2)
%     plot(x, u_eulfor(:,tau))
%     plot(x, u_eulback(:,tau))
    axis([0 1 0 1])
%     labels=[labels;sprintf("time(%d)=%d",tau,t(tau))];
end
xlabel('x')
ylabel('u')
% legend(labels=labels)
h=get(gcf, "currentaxes");
set(h, "fontsize", 12, "linewidth", 1);
cbar=colorbar()
caxis([0,T])
ylabel(cbar,'t','Rotation',0)
title(sprintf('EF, D^{3/2}=-(D^3)^{1/2}) \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
filename=sprintf('Ex61/Ex61-EF-posDer-multi-t-delx=%d.png',delta_x);
% title(sprintf('EF, D^{3/2}=-(-D^3)^{1/2}) \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
% filename=sprintf('Ex61/Ex61-EF-minusDer-multi-t-delx=%d.png',delta_x);
% title(sprintf('EF, D^{3/2}=-((D^3)^{1/2}+(-D^3)^{1/2})/sqrt(2) \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
% filename=sprintf('Ex61/Ex61-EF-mixedDer-multi-t-delx=%d.png',delta_x);
print(1,filename, '-dpng','-S600,400')
clf
hold on 
% Plot3D(x_temp,y_temp,z_temp)
cmap=colormap();
for ti=1:size(listoftimes)(2)
    tau=listoftimes(ti);
%     Plot3D(x,t,real(u_eulback'))
%     xlabel('x')
%     ylabel('t')
    colour_i=round(1+(tau*63/Nt));
    plot(x, u_expm(:,tau),'Color', color=cmap(colour_i,:),'linewidth',2)
%     plot(x, u_eulfor(:,tau))
%     plot(x, u_eulback(:,tau))
    axis([0 1 0 1])
%     labels=[labels;sprintf("time(%d)=%d",tau,t(tau))];
end
xlabel('x')
ylabel('u')
% legend(labels=labels)
h=get(gcf, "currentaxes");
set(h, "fontsize", 12, "linewidth", 1);
cbar=colorbar()
caxis([0,T])
ylabel(cbar,'t','Rotation',0)
title(sprintf('EB, D^{3/2}=-(D^3)^{1/2}) \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
filename=sprintf('Ex61/Ex61-EB-posDer-multi-t-delx=%d.png',delta_x);
% title(sprintf('EB, D^{3/2}=-(-D^3)^{1/2}) \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
% filename=sprintf('Ex61/Ex61-EB-minusDer-multi-t-delx=%d.png',delta_x);
% title(sprintf('EB, D^{3/2}=-((D^3)^{1/2}+(-D^3)^{1/2})/sqrt(2) \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
% filename=sprintf('Ex61/Ex61-EB-mixedDer-multi-t-delx=%d.png',delta_x);
print(1,filename, '-dpng','-S600,400')

