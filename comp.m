1;



%%%%%%%%%%%%%% Initial Values
% Here we define our global values, the total time, the number of time steps
% and the space interval size
% Here is various methods to up the list of delta x's 
% to loop through

%%% Set number of steps and uper and lower bounds
% Nsteps=100;
% delta_xList=linspace(0.01,0.1,Nsteps);

% Similar but with two different linspaces, alowing us to probe across a wider range of delta x's
delta_xList=[linspace(0.001,0.01,10),linspace(0.01,0.1,10)];
Nsteps=size(delta_xList)(2);
differenceEXPM=zeros(1,Nsteps);
for j=1:Nsteps
    delta_x=delta_xList(j)

% delta_x=0.01;
    T=0.03;
    Nt=81;
    delta_t=T/(Nt-1);
    stability=delta_t/(delta_x)^(3/2);
    if stability>0.5
        printf("EF Unstable, stability ratio: %.1f\n", stability)
    end

    Nx=int32(1+(1/delta_x));
    t=linspace(0,T,Nt);
    x=linspace(0,1,Nx);



%%%%%%%%%%%%%% Defining the finite-differnce Matrices
%     Here we define the finite-differnce matrix we want to use in this problem. 
%     For later convience the Periodic condisiton is
%     also included
%     D3
%     This consists of 0 along the diagonal, ±2's on the1-off diagonal and  ±1 2-off diagonals.
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

    D23over4=-sqrtm(sqrtm((-D2)^3));






%%%%%%%%%%%%%% Initial Conditions
%   We initialise our three integration method solution 
%   matrices. 
    u_initial=zeros(Nx,1);
    for i=1:Nx
        u_initial(i)=exp(-100*(x(i)-0.5)^2);
    end;


%%%%%%%%%%%%%% Final Conditions
%   We find the vectors at the final T
    uofT=expm(t(end).*D3over2)*u_initial;
    vofT=expm(t(end).*D23over4)*u_initial;
%   Compute their differeces
    differ=uofT - vofT;
    differenceEXPM(j)=sqrt(differ'*differ);
    
end

%%%%%%%%%%%%%% Plotting
% Some of this is not used, and is keep as archival for future usage if required
% Also note the plots require a folder Comp to be present.

% Difference Plotting
clf
hold on
cmap=colormap
loglog(delta_xList, differenceEXPM, 'Color', color=cmap(1,:), ':o', 'linewidth', 2)
legend(labels=['EXPM'])
xlabel('\Delta x')
ylabel('||u-v||')
title('Comparing -(sqrt(D3)+sqrt(-D3))/sqrt(2) vs. -sqrt(sqrt((-D2)3))')
print(1,'Comp/Comp-diff-EXPM.png', '-dpng','-S600,400')



%%% Colorplot EXPM
% Plotting of u(x,t)
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
% 
title(sprintf('EXPM, D^{3/2}=-((D^3)^{1/2}+(-D^3)^{1/2})/sqrt(2) \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
filename=sprintf('Comp/Ex61-EXPM-mixedDer-multi-t-delx=%d.png',delta_x);

print(1,filename, '-dpng','-S600,400')



% Plotting of v(x,t)
clf
hold on 
listoftimes=1:(Nt-1)/5:Nt;
cmap=colormap();
for ti=1:size(listoftimes)(2)
    tau=listoftimes(ti);
    colour_i=round(1+(tau*63/Nt));
    plot(x, v_expm(:,tau),'Color', color=cmap(colour_i,:),'linewidth',2)
    axis([0 1 0 1])
end
xlabel('x')
ylabel('v')
h=get(gcf, "currentaxes");
set(h, "fontsize", 12, "linewidth", 1);
cbar=colorbar();
caxis([0,T])
ylabel(cbar,'t','Rotation',0)

title(sprintf('EXPM, D^2, Periodic BCs, \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
filename=sprintf('Comp/Comp-EXPM-2D3over4-PBCs-multi-t-delx=%d.png',delta_x);
print(1,filename, '-dpng','-S600,400')
