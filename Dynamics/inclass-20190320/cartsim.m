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

