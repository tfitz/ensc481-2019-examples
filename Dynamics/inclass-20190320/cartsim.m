%% cartsim.m
% simulating the pendulum on a cart

clear all
close all
clc


%% Define the parameters
g = 9.81;
m = 2;
M = 5;
l = 0.5;
J_G = 1/12*m*l^2;

%%
% Initial conditions
tf = 10;
z0 = [0, 30*pi/180, 0, 0]';

%% Define the actuation
% No forcing:
f = @(t,z) 0;
torque = @(t,z) 0;

%% Solve the system
sol = ode45( @(t,z) cart(t,z, m, g, l, M, J_G, f, torque), ...
    [0,tf], z0);

%% Now let's plot
time = sol.x;
states = sol.y;

plot(time, states)
xlabel('time')

%%
% Now, let's plot the Hamiltonian
figure
H = hamiltonian(sol.y, m, g, l, M, J_G);
plot(time, H/H(1)- 1)
xlabel('time')
ylabel('Normalized total energy H(t)/H(t=0) - 1')

%% now let's make an animation
fig = figure();
ax = axes('Parent', fig);
hold(ax,'on');
ax.DataAspectRatio = [1,1,1];
xlabel(ax,'x');
ylabel(ax,'y');
xlim(ax,[-2,2]);
ylim(ax,[-2,2]);

f_cart_x = @(z) [-1, 1, 1, -1] + z(1);
h_cart = patch(ax,'XData', f_cart_x(z0), ...
        'YData', [-1,-1,1,1]*l/10,...
        'FaceAlpha', 0.3 );

f_com_x = @(z) l*sin(z(2)) + z(1);
f_com_y = @(z) l*cos(z(2));
h_com = plot(ax, f_com_x(z0), f_com_y(z0), ...
            'Marker', 'o', 'MarkerFaceColor', 'r', ...
            'MarkerSize', 12);
    

nframes = 300;
time = linspace(0,tf,nframes);

for j = 1:nframes
    
    z = deval(sol, time(j), 1:2);
    
    h_cart.XData = f_cart_x(z);
    h_com.XData = f_com_x(z);
    h_com.YData = f_com_y(z);
    
    pause(0.05);
    
end








