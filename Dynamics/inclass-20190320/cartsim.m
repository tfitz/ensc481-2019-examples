%% cartsim.m

clear all
close all
clc


%%

g = 9.81;
m = 2;
M = 5;
l = 0.5;
J_G = 1/12*m*l^2;

f = @(t,z) 0;
torque = @(t,z) 0;

%%
tf = 10;
z0 = [0, 30*pi/180, 0, 0]';

sol = ode45( @(t,z) cart(t,z, m, g, l, M, J_G, f, torque), ...
    [0,tf], z0);

%% Now let's plot
time = sol.x;
states = sol.y;


plot(time, states)


%%
% Now, let's plot the Hamiltonian
figure
H = hamiltonian(sol.y, m, g, l, M, J_G);
plot(time, H/H(1))

%% now let's make an animation





