1;


%% Exact Solution for Problem 4.1
function u=HeatEquationExact(x,t)
    u=exp(-(pi^2)*t).*sin(pi*x);
endfunction



%%%%%%%%%%%%%% Initial Values
% Here we define our global values, the total time, the number of time steps
% and the space interval size
T=0.1;  
Nt=6;
delta_x=1/6;
% Nt=21;
% delta_x=0.1;
delta_t=T/(Nt-1);

% We use these to generate the number of space steps, then define a 'continuous' and 
% discrete vector for space and time
Nx=int32(1+(1/delta_x));
x_continuous=linspace(0,1);
t_continuous=linspace(0,T)
t=linspace(0,T,Nt);
x=linspace(0,1,Nx);


%%%%%%%%%%%%%% Stability check
% Here we check the stability for the EF method; not this only warns when the 
% code is unstable, it will still run
stability=delta_t/(delta_x)^2
if stability>0.5
    printf("EF Unstable, stability ratio: %.1f\n", stability)
end





%%%%%%%%%%%%%% Defining the finite-differnce Matrices
% Here we define the finite-differnce matrix we want to use in this problem. This consists of 
% -2 along the diagonal and 1's on the two 1-off diagonals. We also use fied boundary values, 
% so the first and final rows need to be zerod. For later convience the Periodic condisiton is
% also included
der_a=-2;   % Diagonal
der_b=1;    % Right-1-off Diagonal
der_c=1;    % Left-1-off Diagonal
der_d=0;    % Left-2-off Diagonal
der_f=0;    % left-2-off Diagonal
D2 = diag(der_a*ones(1,Nx)) + diag(der_b*ones(1,Nx-1),1) +diag(der_f*ones(1,Nx-2),2)+diag(der_d*ones(1,Nx-2),-2)+ diag(der_c*ones(1,Nx-1),-1);
% D2(1,Nx)=der_c;     % P.B.C.s
% D2(Nx,1)=der_b;     % P.B.C.s
D2(1,:)=0;          % F.B.C.s
D2(Nx,:)=0;         % F.B.C.s
D2=D2 / delta_x^2; 




%%%%%%%%%%%%%% Initial Conditions
% We initialise our three integration method solution matrices. 
u_expm=zeros(Nx,Nt);
u_eulfor=zeros(Nx,Nt);
u_eulback=zeros(Nx,Nt);
u_initial=sin(pi*x);
u_eulfor(:,1)=u_initial;
u_eulback(:,1)=u_initial;




%%%%%%%%%%%%%% Integration Steps
% EXPM - Solve the ODE system analytically via the EXPM
for i=1:Nt
    u_expm(:,i)=(expm(t(i).*D2)*u_initial')';
end

% Euler Forward - Solve using the EF integration techniques
for i=1:Nt-1
    u_eulfor(:,i+1) = ((eye(Nx)+delta_t.*D2)*u_eulfor(:,i));
endfor 

% Euler Backward - Solve using the EB integration techniques
for i=1:Nt-1
    u_eulback(:,i+1) = ((eye(Nx)-delta_t.*D2)\u_eulback(:,i));
endfor 


%%%%%%%%%%%%%% Plotting
% Some of this is not used, and is keep as archival for future usage if required
% Also note the plots require a folder Ex41 to be present.
listoftimes=1:round((Nt-1)/5):Nt
listofspace=1:round((Nx-1)/5):Nx
h=get(gcf, "currentaxes");
set(h, "fontsize", 12, "linewidth", 1);
clf
hold on
labels=[''];
for ti=1:size(listoftimes)(2)
    tau=listoftimes(ti);
    filename=sprintf("Ex41/Ex41eulfor-delx=%d.png", delta_x );
    
    % Plot EF solution
    plot(x, u_eulfor(:,tau))
    labels=[labels;sprintf("t=%d",t(tau))];


    title(sprintf('Euler Forward: \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt))
    legend(labels=labels)
    axis([0 1 0 1])

    saveas(1, filename)
end
% Similar Plots bellow
% clf
% hold on
% labels=[''];
% for ti=1:size(listoftimes)(2)
%     tau=listoftimes(ti);
%     filename=sprintf("Ex41/Ex41eulback-delx=%d.png", delta_x );
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
%     filename=sprintf("Ex41/Ex41expm-delx=%d.png", delta_x);
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


%%%%%%%%%%%%%% Colour Plotting
%%% EXPM Colorplot
clf
hold on 
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
cbar=colorbar()
caxis([0,T])
ylabel(cbar,'t','Rotation',0)
title(sprintf('EXPM, \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
filename=sprintf('Ex41/Ex41-EXPM-multi-t-delx=%d.png',delta_x);
print(1,filename, '-dpng','-S1200,800')

% Similar plots below
%%% EF Colorplot
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
% cbar=colorbar()
% caxis([0,T])
% ylabel(cbar,'t','Rotation',0)
% title(sprintf('EF, \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
% filename=sprintf('Ex41/Ex41-EF-multi-t-delx=%d.png',delta_x);
% print(1,filename, '-dpng','-S1200,800')
% 

%%% EB Colorplot
% clf
% hold on 
% cmap=colormap();
% for ti=1:size(listoftimes)(2)
%     tau=listoftimes(ti);
%     colour_i=round(1+(tau*63/Nt));
% 
%     plot(x, u_eulback(:,tau),'Color', color=cmap(colour_i,:),'linewidth',2)
% 
% end
% xlabel('x')
% ylabel('u')
% axis([0 1 0 1])
% h=get(gcf, "currentaxes");
% set(h, "fontsize", 12, "linewidth", 1);
% cbar=colorbar()
% caxis([0,T])
% ylabel(cbar,'t','Rotation',0)
% title(sprintf('EB, \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
% filename=sprintf('Ex41/Ex41-EB-multi-t-delx=%d.png',delta_x);
% print(1,filename, '-dpng','-S1200,800')
% 

%%% Fix Time
labels=['']
for ti=2:size(listoftimes)(2)
    clf
    hold on
    labels=[''];
    tau=listofspace(ti)
    plot(x, u_eulfor(:,tau),'linewidth',2)
    labels=[labels;"EF"];

    plot(x, u_eulback(:,tau),'linewidth',2)
    labels=[labels;"EB"];

    plot(x, u_expm(:,tau),'linewidth',2)
    labels=[labels;"EXPM"];

    plot(x_continuous ,HeatEquationExact(x_continuous,t(tau)),":",'linewidth',2)
    labels=[labels;"EXACT"];

    filename=sprintf("Ex41/Ex41-fixtime-delx=%d-x=%.2f.png", delta_x, x(tau));
    title(sprintf('Fix Time, \\Delta x=%d, T=%d, Nt=%d, t=%.02f', delta_x, T, Nt, t(tau)))
    legend(labels=labels)
    xlabel('x')
    ylabel('u')

    print(1,filename, '-dpng','-S1200,800')
end

%%% Fix Space
labels=['']
for ti=2:size(listoftimes)(2)
    clf
    hold on
    labels=[''];
    tau=listofspace(ti)
    plot(t, u_eulfor(tau,:),'linewidth',2)
    labels=[labels;"EF"];

    plot(t, u_eulback(tau,:),'linewidth',2)
    labels=[labels;"EB"];

    plot(t, u_expm(tau,:),'linewidth',2)
    labels=[labels;"EXPM"];

    plot(t_continuous ,HeatEquationExact(x(tau),t_continuous),":",'linewidth',2)
    labels=[labels;"EXACT"];

    filename=sprintf("Ex41/Ex41-fixspace-delx=%d-t=%.2f.png", delta_x, x(tau));
    title(sprintf('Fix Space, \\Delta x=%d, T=%d, Nt=%d, x=%.02f', delta_x, T, Nt, x(tau)))
    legend(labels=labels)
    xlabel('t')
    ylabel('u')

    print(1,filename, '-dpng','-S1200,800')
end
