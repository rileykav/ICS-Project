1;


%%%%%%%%%%%%%% Initial Values
% Here we define our global values, the total time, the number of time steps
% and the space interval size
T=0.03;
Nt=101;
delta_x=0.01;
delta_t=T/(Nt-1);


% We use these to generate the number of space steps, then define a 
% discrete vector for space and time
Nx=int32(1+(1/delta_x));
t=linspace(0,T,Nt);
x=linspace(0,1,Nx);

%%%%%%%%%%%%%% Stability check
% Here we check the stability for the EF method; not this only warns when the 
% code is unstable, it will still run
stability=delta_t/(delta_x)^(3/2);
if stability>0.5
    printf("EF Unstable, stability ratio: %.1f\n", stability)
end




%%%%%%%%%%%%%% Defining the finite-differnce Matrices
% Here we define the finite-differnce matrix we want to use in this problem. 
% For later convience the Periodic condisiton is
% also included
% D3
% This consists of 0 along the diagonal, ±2's on the1-off diagonal and  ±1 2-off diagonals.
a=0;    % Diagonal
b=-1;   % Right-1-diagonal
c=1;    % Left-1-diagonal
d=-1/2; % Left-2-diagonal
f=1/2;  % Right-2-diagonal
D3 = diag(a*ones(1,Nx)) + diag(b*ones(1,Nx-1),1) +diag(f*ones(1,Nx-2),2)+diag(d*ones(1,Nx-2),-2)+ diag(c*ones(1,Nx-1),-1);
% B.C.
D3(1,Nx)=c;
D3(1,Nx-1)=d;
D3(2,Nx)=d;
D3(Nx,1)=b;
D3(Nx,2)=f;
D3(Nx-1,1)=f;
D3;
D3=D3/(delta_x^3);

% (a)
D3over2=-sqrtm(D3);
%
% (b)
D3over2=-sqrtm(-D3);

% (c)
D3over2=-(sqrtm(D3) + sqrtm(-D3))/sqrt(2);


%% D2
% This consists of -2 along the diagonal and 1's on the two 1-off diagonals.
der_a=-2;   % Diagonal
der_b=1;    % Right-1-off Diagonal
der_c=1;    % Left-1-off Diagonal
der_d=0;    % Left-2-off Diagonal
der_f=0;    % left-2-off Diagonal
D2 = diag(der_a*ones(1,Nx)) + diag(der_b*ones(1,Nx-1),1) +diag(der_f*ones(1,Nx-2),2)+diag(der_d*ones(1,Nx-2),-2)+ diag(der_c*ones(1,Nx-1),-1);
D2(1,Nx)=der_c;
D2(Nx,1)=der_b;
D2=D2 / delta_x^2;
D2Fixed=D2;
D2Fixed(1,:)=0;
D2Fixed(Nx,:)=0;

% (Comparison)
% D3over2=D2


%%%%%%%%%%%%%% Initial Conditions
% We initialise our three integration method solution matrices. 
u_initial=zeros(Nx,1);
for i=1:Nx
    u_initial(i)=exp(-100*(x(i)-0.5)^2);
end
u_expm=zeros(Nx,Nt);
u_eulfor=zeros(Nx,Nt);
u_eulback=zeros(Nx,Nt);
u_eulfor(:,1)=u_initial;
u_eulback(:,1)=u_initial;


%%%%%%%%%%%%%% Integration Steps
% EXPM - Solve the ODE system analytically via the EXPM
for i=1:Nt
    u_expm(:,i)=(expm(t(i).*D3over2)*u_initial)
end

% Euler Forward - Solve using the EF integration techniques
for i=1:Nt-1
    u_eulfor(:,i+1)=(eye(Nx)+delta_t.*D3over2)*u_eulfor(:,i);
end

% Euler Backward - Solve using the EB integration techniques
for i=1:Nt-1
    A=eye(Nx)-delta_t.*D3over2;
    u_eulback(:,i+1)=A\u_eulfor(:,i);
end



%%%%%%%%%%%%%% Plotting
% Some of this is not used, and is keep as archival for future usage if required
% Also note the plots require a folder Ex61 to be present.
%%% Colorplot EXPM
clf
hold on 
listoftimes=1:(Nt-1)/5:Nt;
cmap=colormap();
for ti=1:size(listoftimes)(2)
    tau=listoftimes(ti);
    colour_i=round(1+(tau*63/Nt));
    plot(x, u_expm(:,tau),'Color', color=cmap(colour_i,:),'linewidth',2)
    axis([0 1 0 1])
end
xlabel('x')
ylabel('u')
h=get(gcf, "currentaxes");
set(h, "fontsize", 12, "linewidth", 1);
cbar=colorbar();
caxis([0,T])
ylabel(cbar,'t','Rotation',0)

% title(sprintf('EXPM, D^{3/2}=-(D^3)^{1/2} \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
% filename=sprintf('Ex61/Ex61-EXPM-posDer-multi-t-delx=%d.png',delta_x);
% title(sprintf('EXPM, D^{3/2}=-(-D^3)^{1/2} \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));;
% filename=sprintf('Ex61/Ex61-EXPM-minusDer-multi-t-delx=%d.png',delta_x);
% title(sprintf('EXPM, D^{3/2}=-((D^3)^{1/2}+(-D^3)^{1/2})/sqrt(2) \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
% filename=sprintf('Ex61/Ex61-EXPM-mixedDer-multi-t-delx=%d.png',delta_x);
% 

title(sprintf('EXPM, D^2, Periodic BCs, \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
filename=sprintf('Ex61/Ex61-EXPM-2Der-PBCs-multi-t-delx=%d.png',delta_x);


print(1,filename, '-dpng','-S600,400')


%%% Colorplot EF
% clf
% hold on 
% cmap=colormap();
% for ti=1:size(listoftimes)(2)
%     tau=listoftimes(ti);
%     colour_i=round(1+(tau*63/Nt));
%     plot(x, u_eulfor(:,tau),'Color', color=cmap(colour_i,:),'linewidth',2)
%     axis([0 1 0 1])
% end
% xlabel('x')
% ylabel('u')
% h=get(gcf, "currentaxes");
% set(h, "fontsize", 12, "linewidth", 1);
% cbar=colorbar();
% caxis([0,T])
% ylabel(cbar,'t','Rotation',0)

% title(sprintf('EF, D^{3/2}=-(D^3)^{1/2} \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
% filename=sprintf('Ex61/Ex61-EF-posDer-multi-t-delx=%d.png',delta_x);
% title(sprintf('EF, D^{3/2}=-(-D^3)^{1/2} \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
% filename=sprintf('Ex61/Ex61-EF-minusDer-multi-t-delx=%d.png',delta_x);
% title(sprintf('EF, D^{3/2}=-((D^3)^{1/2}+(-D^3)^{1/2})/sqrt(2) \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
% filename=sprintf('Ex61/Ex61-EF-mixedDer-multi-t-delx=%d.png',delta_x);
% 
% print(1,filename, '-dpng','-S600,400')



%%% Colorplot EB
% clf
% hold on 
% cmap=colormap();
% for ti=1:size(listoftimes)(2)
%     tau=listoftimes(ti);
%     colour_i=round(1+(tau*63/Nt));
%     plot(x, u_expm(:,tau),'Color', color=cmap(colour_i,:),'linewidth',2)
%     axis([0 1 0 1])
% end
% xlabel('x')
% ylabel('u')
% h=get(gcf, "currentaxes");
% set(h, "fontsize", 12, "linewidth", 1);
% cbar=colorbar();
% caxis([0,T])
% ylabel(cbar,'t','Rotation',0)
% 
% title(sprintf('EB, D^{3/2}=-(D^3)^{1/2} \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
% filename=sprintf('Ex61/Ex61-EB-posDer-multi-t-delx=%d.png',delta_x);
% title(sprintf('EB, D^{3/2}=-(-D^3)^{1/2} \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
% filename=sprintf('Ex61/Ex61-EB-minusDer-multi-t-delx=%d.png',delta_x);
% title(sprintf('EB, D^{3/2}=-((D^3)^{1/2}+(-D^3)^{1/2})/sqrt(2) \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
% filename=sprintf('Ex61/Ex61-EB-mixedDer-multi-t-delx=%d.png',delta_x);
% 
% print(1,filename, '-dpng','-S600,400')

% listoftimes=1:(Nt-1)/10:Nt;
% clf
% hold on
% for i=1:10
%     tau=listoftimes(i)
%     plot(x,u_eulfor(:,tau));
% end
% title("\\Delta x=0.01")
% h=get(gcf, "currentaxes");
% set(h, "fontsize", 12, "linewidth", 1);
% print(1,'diff-openbound-D3over2.png', '-dpng','-S1200,800')
% 
