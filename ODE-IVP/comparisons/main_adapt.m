%% RK Adapt. Comparisons
% 

%%
clear all
close all
clc

%% Define sample problem
tol = 0.001;
h0 = 0.01;
tf = 100;

%%
% The most interesting part to vary here is ratio of wn and stepsize h.
wn = 1;
zeta = 0.1;

% The equation to solve is
% $$ \ddot{x} + 2 \zeta \omega_n \dot{x} + \omega_n^2 x = 0 $$
f = @(t,z) [z(2); -wn^2*z(1)-2*zeta*wn*z(2)];

%%
% With initial conditions
z0 = [1;0];

%%
% Let's also define the true solution
wd = wn*sqrt(1-zeta^2);
y = @(t) exp(-zeta*wn*t).*( z0(1)*cos(wd*t) + (z0(2) + zeta*wn*z0(1))/wd*sin(wd*t) );

%% Visualize the results 
% Now we can make a plot to compare the true solution and our
% approximations
fig1 = figure();
ax1 = axes('Parent', fig1);
time = linspace(0,tf,300);
plot(ax1,time, y(time), '-', ...
    'LineWidth', 2, 'DisplayName', 'True solution')
hold(ax1,'on');
xlabel(ax1,'time')
ylabel(ax1,'x(t)')
legend(ax1,'show')

fig2 = figure();
ax2 = axes('Parent', fig2);

%% RK2 with simple time adaptivity
tic
[Z1, t1, nfail1] = rk2_adapt(f, 0, h0, tf, z0, tol);
toc

plot( ax1, t1, Z1(1,:), 'd', ...
    'LineWidth', 3, 'DisplayName', 'RK2 simple' )

plot( ax2, t1(1:end-1), diff(t1), ...
    'LineWidth', 3, ...
    'DisplayName', sprintf('1: %3d steps, %3d fails', length(t1)-1, nfail1) )
hold(ax2, 'on')
xlabel(ax2,'time t')
ylabel(ax2,'h')
legend(ax2,'show')

%% RK2 with method 2
tic
[Z2, t2, nfail2] = rk2_adapt_method2(f, 0, h0, tf, z0, tol);
toc

plot( ax1, t2, Z2(1,:), '*', ...
    'LineWidth', 3, 'DisplayName', 'RK2 method 2' )

plot( ax2, t2(1:end-1), diff(t2), ...
    'LineWidth', 3, ...
    'DisplayName', sprintf('2: %3d steps, %3d fails', length(t2)-1, nfail2) )


%% RK2 with method 3
tic
[Z3, t3, nfail3] = rk2_adapt_method3(f, 0, h0, tf, z0);
toc

plot( ax1, t3, Z3(1,:), 's', ...
    'LineWidth', 3, 'DisplayName', 'RK2 method 3' )

plot( ax2, t3(1:end-1), diff(t3), ...
    'LineWidth', 3, ...
    'DisplayName', sprintf('3: %3d steps, %3d fails', length(t3)-1, nfail3) )

