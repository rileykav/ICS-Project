
1;
function exact_solution=Heat_eq(t,x)

    exact_solution = exp(-(pi^2)*t).*sin(pi*x);

endfunction

x_fixed = [0.5,1,1.5,2,2.5]/100
x_intervals=100;

t_interval = 100;
t_values = linspace(0,1,t_interval);

u_exact = zeros(t_interval);

for i=1:size(x_fixed)(2)

    u_exact = Heat_eq(t_values, x_fixed(i));
    % plot(t_values, u_exact)
    % hold on;

end


t_fixed = [0.5,1,1.5,2,2.5]/100

x_interval = 10;
x_values = linspace(0,1,x_interval);

u_exact = zeros(x_interval);

for i=1:size(t_fixed)(2)

    u_exact = Heat_eq(t_fixed(i), x_values);
    % plot(x_values, u_exact)
    % hold on;

end


% Approximations

t_max = 0.1
delta_x = 0.1

N = 1 + 1/delta_x
a = -2;
b = 1;
c = 1;
D = diag(a*ones(1,N)) + diag(b*ones(1,N-1),1) + diag(c*ones(1,N-1),-1);
D(1,:)=0;
D(N,:)=0;

t_interval = 11;
t = linspace(0,t_max,t_interval)
delta_t = t_max/(t_interval-1)

D = D/(delta_x^2)

u_initial = zeros(N,1);
x_interval = linspace(0,1,N);

u_initial = sin(pi*x_interval)

% Solution exponentiating D2 matrix

u_approx1=zeros(N,t_interval);

for i=1:t_interval

    u_approx1(:,i) = (expm(t(i).*D)*(u_initial)')';

    (expm(t(i).*D));

    u_initial;

end


clf
% hold on
plot(x_interval, u_approx1(:,3))

% Euler - Forward

u_approx2=zeros(N,t_interval);

    u_approx2(:,1) = u_initial;

    u_approx2;

    for j=1:(t_interval)

        j;

        u_approx2(:,j);

        u_approx2(:,j+1) = (eye(N)+delta_t.*D)*u_approx2(:,j);

        u_approx2(4,j);

    end

    u_approx2

%     for i=1:N

%         u_approx2(i,10)
%     end
% clf

%     plot(x_interval, u_approx2(:,1))

% Euler - Backwards

u_approx3=zeros(N,t_interval);

    u_approx3(:,1) = u_initial;

    u_approx3;

    for j=1:(t_interval)

        j;

        u_approx3(:,j);

        A = eye(N) - delta_t*D;

        u_approx3(:,j+1) = A\u_approx3(:,j);

    end
t(10)

% u_approx3;
% clf
% hold on
% plot(x_interval, u_approx1(:,10))
% plot(x_interval, u_approx2(:,10))
% plot(x_interval, u_approx3(:,10))
% 
% x_interval=linspace(0,1,N);
% plot(x_interval, Heat_eq(t(10),x_interval))
% 
% labels=["Exp"; "EF"; "EB"; "Exact"];
% 
% legend(labels=labels)
% 


listoftimes=([0,1,2,3,4,5,6,7,8,9,10])*((t_interval-1)/10)+1
clf
hold on
labels=[''];
for ti=1:size(listoftimes)(2)
    tau=listoftimes(ti);
    filename=sprintf("Ex41/Monica-ex41eulback.png");

    plot(x_interval, u_approx2(:,tau))
    labels=[labels;sprintf("time(%d)=%d",tau,t(tau))];


    title(sprintf('Monica Euler Forward: \\Delta x=%d, T=%d, Nt=%d', delta_x, T, Nt))
    legend(labels=labels)
    axis([0 1 0 1])

    saveas(1, filename)
end

