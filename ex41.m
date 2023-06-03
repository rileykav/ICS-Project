1;



function u=HeatEquationExact(x,t)
%%% Solving ax^2+bx+c
    u=exp(-(pi^2)*t).*sin(pi*x);
%     u=x+t;
endfunction

t=[0,1,2,3,4,5,6,7,8,9,10]./100;
x_continuous=linspace(0,1);


delta_x=0.01;
Nt=2001;
% delta_x=0.1;
% Nt=21;
T=0.1;
delta_t=T/(Nt-1);

stability=delta_t/(delta_x)^2
if stability>0.5
    printf("EF Unstable, stability ratio: %.1f\n", stability)
end

Nx=int32(1+(1/delta_x));
t=linspace(0,T,Nt);
x=linspace(0,1,Nx);



der_a=-2;
der_b=1;
der_c=1;
der_d=0;
der_f=0;
D2 = diag(der_a*ones(1,Nx)) + diag(der_b*ones(1,Nx-1),1) +diag(der_f*ones(1,Nx-2),2)+diag(der_d*ones(1,Nx-2),-2)+ diag(der_c*ones(1,Nx-1),-1);
% D2(1,Nx)=der_c;
% D2(Nx,1)=der_b;
D2(1,:)=0;
D2(Nx,:)=0;
D2=D2 / delta_x^2;

u_initial=sin(pi*x);

u_expm=zeros(Nx,Nt);
u_eulfor=zeros(Nx,Nt);
u_eulback=zeros(Nx,Nt);

u_eulfor(:,1)=u_initial;
u_eulback(:,1)=u_initial;

% EXPM 
for i=1:Nt
    u_expm(:,i)=(expm(t(i).*D2)*u_initial')';
end

% Euler Forward
for i=1:Nt-1
    u_eulfor(:,i+1) = ((eye(Nx)+delta_t.*D2)*u_eulfor(:,i));
endfor 

% Euler Backward
for i=1:Nt-1
    u_eulback(:,i+1) = ((eye(Nx)-delta_t.*D2)\u_eulback(:,i));
endfor 


% Plotting
% clf 
% hold on
% listoftimes=([0,1,2,3,4,5,6,7,8,9,10])*((Nt-1)/10)+1;
% h=get(gcf, "currentaxes");
% set(h, "fontsize", 12, "linewidth", 1);
% for ti=1:size(listoftimes)(2)
%     clf
%     hold on
%     labels=[''];
%     tau=listoftimes(ti)
%     filename=sprintf("Ex41/ex41-all-delx=%d-t=%.2f.png", delta_x, t(tau));
%     plot(x, u_eulfor(:,tau))
%     labels=[labels;"EF"];
% 
%     plot(x, u_eulback(:,tau))
%     labels=[labels;"EB"];
% 
%     plot(x, u_expm(:,tau))
%     labels=[labels;"EXPM"];
% 
%     plot(x_continuous,HeatEquationExact(x_continuous,t(tau)))
%     labels=[labels;"EXACT"];
% 
%     title(sprintf('\\Delta x=%d, T=%d, Nt=%d, t=%.02f', delta_x, T, Nt, t(tau)))
%     legend(labels=labels)
%     axis([0 1 0 1])
% 
%     saveas(1, filename)
% end
% 
% 
% clf
% hold on
% labels=[''];
% for ti=1:size(listoftimes)(2)
%     tau=listoftimes(ti);
%     filename=sprintf("Ex41/ex41eulfor-delx=%d.png", delta_x );
% 
%     plot(x, u_eulfor(:,tau))
%     labels=[labels;sprintf("t=%d",t(tau))];
% 
% 
%     title(sprintf('Euler Forward: \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt))
%     legend(labels=labels)
%     axis([0 1 0 1])
% 
%     saveas(1, filename)
% end
% clf
% hold on
% labels=[''];
% for ti=1:size(listoftimes)(2)
%     tau=listoftimes(ti);
%     filename=sprintf("Ex41/ex41eulback-delx=%d.png", delta_x );
% 
%     plot(x, u_eulback(:,tau))
%     labels=[labels;sprintf("t=%d",t(tau))];
% 
% 
%     title(sprintf('Euler Backward: \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt))
%     legend(labels=labels)
%     axis([0 1 0 1])
% 
%     saveas(1, filename)
% end
% 
% clf
% hold on
% labels=[''];
% for ti=1:size(listoftimes)(2)
%     tau=listoftimes(ti);
%     filename=sprintf("Ex41/ex41expm-delx=%d.png", delta_x);
% 
%     plot(x, u_expm(:,tau))
%     labels=[labels;sprintf("t=%d",t(tau))];
% 
% 
%     title(sprintf('Expm: \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt))
%     legend(labels=labels)
%     axis([0 1 0 1])
% 
%     saveas(1, filename)
% end
listoftimes=1:round((Nt-1)/5):Nt
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
title(sprintf('EXPM, \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
filename=sprintf('Ex41/Ex41-EXPM-multi-t-delx=%d.png',delta_x);
print(1,filename, '-dpng','-S1200,800')
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
% legend(labels=labels)
xlabel('x')
ylabel('u')
h=get(gcf, "currentaxes");
set(h, "fontsize", 12, "linewidth", 1);
cbar=colorbar()
caxis([0,T])
ylabel(cbar,'t','Rotation',0)
title(sprintf('EF, \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
filename=sprintf('Ex41/Ex41-EF-multi-t-delx=%d.png',delta_x);
print(1,filename, '-dpng','-S1200,800')
clf
hold on 
cmap=colormap();
for ti=1:size(listoftimes)(2)
    tau=listoftimes(ti);
    colour_i=round(1+(tau*63/Nt));

    plot(x, u_eulback(:,tau),'Color', color=cmap(colour_i,:),'linewidth',2)

end
xlabel('x')
ylabel('u')
axis([0 1 0 1])
h=get(gcf, "currentaxes");
set(h, "fontsize", 12, "linewidth", 1);
cbar=colorbar()
caxis([0,T])
ylabel(cbar,'t','Rotation',0)
title(sprintf('EB, \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
filename=sprintf('Ex41/Ex41-EB-multi-t-delx=%d.png',delta_x);
print(1,filename, '-dpng','-S1200,800')
