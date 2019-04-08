%% Linear Multi-step methods
% 3-Apr: Adams-Bashforth 2
% 5-Apr: add Trapezoidal Rule


%%
clear all
close all
clc

%% Define sample problem
% The most interesting part to vary here is ratio of wn and stepsize h.
wn = 1;
zeta = 0.1;
h = 0.4;

%%
% The equation to solve is
% $$ \ddot{x} + 2 \zeta \omega_n \dot{x} + \omega_n^2 x = 0 $$
f = @(z) [z(2); -wn^2*z(1)-2*zeta*wn*z(2)];

%%
% With initial conditions
z0 = [1;0];

%%
% Let's also define the true solution
wd = wn*sqrt(1-zeta^2);
y = @(t) exp(-zeta*wn*t).*( z0(1)*cos(wd*t) + (z0(2) + zeta*wn*z0(1))/wd*sin(wd*t) );

%%
tf = 10;
time = 0:h:tf;
nstep = length(time);

Z = zeros( 2, nstep);

%% Integrate the solution
% Let's first take an Euler step since AB2 is not self-starting
tic
Z(:,1) = z0;
Z(:,2) = Z(:,1) + h*f(z0);

%% 
% Now we can continue with AB2
fm1 = f(Z(:,1));
for i = 2:nstep-1
    
    fi = f(Z(:,i));
    Z(:,i+1) = Z(:,i) + h*( 3/2*fi - 1/2*fm1 );
    
    fm1 = fi;
    
end
toc

%% AM2/Trap
Z_am2 = zeros(size(Z));
Z_am2(:,1) = z0;
tol = 1e-6;
tic
for i = 1:nstep-1
    
    iter = 0;
    err = 1;
    if i == 1
        zest = Z_am2(:,i) + h*f(Z_am2(:,i));
    else
        zest = Z_am2(:,i) + h*(3/2*f(Z_am2(:,i)) - 1/2*f(Z_am2(:,i-1)) );
    end
    while err >= tol
        iter = iter + 1;
        zest1 = Z_am2(:,i) + h/2*( f(zest) + f(Z_am2(:,i)) );
    
        err = norm(zest1 - zest);
        zest = zest1;
        
    end
    Z_am2(:,i+1) = zest;
    iter
    
end
toc

%% Visualize the results 
% Now we can make a plot to compare the true solution and our
% approximations
plot(time, y(time), '-', ...
    'LineWidth', 2, 'DisplayName', 'True solution')
hold on
plot(time, Z(1,:), 'o', ...
    'LineWidth', 2, 'DisplayName', 'AB2')
xlabel('time')
ylabel('x(t)')
legend('show')


plot( time, Z_am2(1,:),'*', 'DisplayName', 'AM2', 'LineWidth', 3)


%%
% Convergence testing example
% z1 = [-0.328569181097212;
%    0.197477255580434];
% z2 = [-0.334745911140231;
%    0.188330026330040];
% z3 = [-0.336317094322793;
%    0.186084168432383];
% 
% z1(1) - z2(1)
% z2(1) - z3(1)


%%