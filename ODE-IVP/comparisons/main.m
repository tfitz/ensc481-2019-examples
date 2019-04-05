%% Linear Multi-step methods
% 3-Apr: Adams-Bashforth 2
%


%%
clear all
close all
clc

%% Define sample problem
% The most interesting part to vary here is ratio of wn and stepsize h.
wn = 1;
zeta = 0.1;
h = 0.2;

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


