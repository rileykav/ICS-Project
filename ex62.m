
% -------------- FUNCTIONS --------------

    % Defining g(u, delta_t, gamma) function to be called in the recursive relation
    1;

    function g_vec=g_funct(u, delta_t, gamma_val)

        g_vec = u + delta_t * gamma_val*(u.*(1-u));

    endfunction

% -------------- MAIN CODE --------------

    % Defining parameters of the simulation: tmax = maximum time, delta_x = space step-size
    % d_val = diffusion coefficient, gamma_val = reaction coefficient

    t_max = 1;
    delta_x = 0.01;
    d_val = 0.001;
    gamma_val = 5;

    % Find number of steps in space N and create vector x with values in space in [0,1] 

    N = 1 + 1/delta_x
    x = linspace(0,1,N);

    % Define matrix D3 that computes third order time deivative with OBC and compute the
    % square roots with sqrtm function

    a = 0.5;
    b = -1;
    c = 1;
    D3= diag(a*ones(1,N-2),2) + diag(b*ones(1,N-1),1) + diag(c*ones(1,N-1),-1) + diag((-a)*ones(1,N-2),-2) + diag((-a)*ones(1,2), N-2) + diag((a)*ones(1,2), 2-N) + diag(c*ones(1,1), N-1) + diag((-c)*ones(1,1), 1-N)
    D3 = D3/(delta_x^3)
    D3sq = sqrtm(D3)
    D3sqmin = sqrtm(-D3)

    % Define matrix D (=D2) that computes second order time deivative with OBC

    a = -2;
    b = 1;
    c = 1;
    D2 = diag(a*ones(1,N)) + diag(b*ones(1,N-1),1) + diag(c*ones(1,N-1),-1);
    D2(1,N)=1;
    D2(N,1)=1;
    D2 = D2/(delta_x)^2;

    % Impose initial condition

    u_initial = zeros(N,1);

    for i=1:N
        u_initial(i) = exp(-50*(x(i)-1/2)^2);
    end

    % Fix number of time-steps t_interval and compute vector t with values of time in [0,t_max]

    t_interval = 500;
    t = linspace(0,t_max,t_interval)
    delta_t = t_max/(t_interval-1)

    % Implement recursive relation in time to find the solution u_approx

    u_approx=zeros(N,size(t)(2));

    u_approx(:,1) = u_initial;

    for j=1:(size(t)(2)-1);

        % Solution with D3 and OBC

            % A = -d_val*(D3sq+D3sqmin)/sqrt(2);

            % A = eye(N) - delta_t*A;

        % Solution with D2 and OBC

            A = eye(N) - delta_t*d_val*(D2);

        u_approx(:,j+1) = A\g_funct(u_approx(:,j), delta_t, gamma_val);

    end

% -------------- PLOTS --------------

    % Plot for a bunch of different times with colour map

    listoftimes = [];

    time_intervals = 10;

    delta_tint = floor(t_interval/time_intervals)

    listoftimes(1) = delta_tint;


    for i=2:time_intervals;
        listoftimes(i) = listoftimes(i-1) + delta_tint;
    end

    listoftimes
    size(listoftimes)(2)

    clf
    hold on 

        cmap=colormap();
    for ti=1:size(listoftimes)(2)
        tau=listoftimes(ti)
        colour_i=round(1+(tau*63/(t_interval)));
        plot(x, u_approx(:,tau),'Color', color=cmap(colour_i,:),'linewidth',2)
        axis([0 1 0 1])
    end

    xlabel('x');
    ylabel('u');
    h=get(gcf, "currentaxes");
    set(h, "fontsize", 12, "linewidth", 1);
    cbar=colorbar();
    caxis([0,t_max]);
    ylabel(cbar,'t','Rotation',0);
    title(sprintf('IMEX \\Deltax=%d, T=%d, Nt=%d', delta_x, t_max, t_interval));
    filename=sprintf('plots_6.2/Ex62multit_Nt=%d_D2.png',t_interval);
    % title(sprintf('EXPM, D^{3/2}=-(-D^3)^{1/2}) \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
    % filename=sprintf('Ex61/Ex61-EXPM-minusDer-multi-t-delx=%d.png',delta_x);
    % title(sprintf('EXPM, D^{3/2}=-((D^3)^{1/2}+(-D^3)^{1/2})/sqrt(2) \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt));
    % filename=sprintf('Ex61/Ex61-EXPM-mixedDer-multi-t-delx=%d.png',delta_x);
    print(1,filename, '-dpng','-S600,400')


    % Print in time for a bunch of space points

    clf
    hold on 

    cmap=colormap();
    listofspace=1
    for xi=1:6

        colour_i=round((listofspace*63/(N/2)))
        plot(t, u_approx(listofspace,:),'Color', color=cmap(colour_i,:),'linewidth',2)

        axis([0 1 0 1])

        listofspace = listofspace + (N-1)/10
    end

    xlabel('t');
    ylabel('u');
    h=get(gcf, "currentaxes");
    set(h, "fontsize", 12, "linewidth", 1);
    cbar=colorbar();
    caxis([0,0.5]);
    ylabel(cbar,'t','Rotation',0);
    title(sprintf('IMEX \\Deltax=%d, T=%d, Nt=%d', delta_x, t_max, t_interval));
    filename=sprintf('plots_6.2/Ex62multix_Nt=%d_first1_D2.png',t_interval);
    print(1,filename, '-dpng','-S600,400')


    clf
    hold on 

    cmap=colormap();

    listofspace = 51

    for xi=1:6

        colour_i=round(((listofspace-50)*63/(N/2)))
        plot(t, u_approx(listofspace,:),'Color', color=cmap(colour_i,:),'linewidth',2)

        axis([0 1 0 1])

        listofspace = listofspace + (N-1)/10
    end

    xlabel('t');
    ylabel('u');
    h=get(gcf, "currentaxes");
    set(h, "fontsize", 12, "linewidth", 1);
    cbar=colorbar();
    caxis([0.5,1]);
    ylabel(cbar,'x','Rotation',0);
    title(sprintf('IMEX \\Deltax=%d, T=%d, Nt=%d', delta_x, t_max, t_interval));
    filename=sprintf('plots_6.2/Ex62multix_Nt=%d_second_D2.png',t_interval);
    print(1,filename, '-dpng','-S600,400')


