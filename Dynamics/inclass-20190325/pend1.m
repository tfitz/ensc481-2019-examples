%% Penulum as a constrained system
% Example built inclass, on 25 Mar 2019

clear all
close all
clc


%% Define some parameters
m = 1;
g = 9.81;
l = 2;

%%
% Setup initial conditions, this is where the length matters
z0 = [l*cos(0); l*sin(0); 0; 0];


%%
% Solve using ode45
sol = ode45(@(t,z) constrained_ode(t,z,g), [0,10], z0);

%% State-space function to solve
% this was derived in the jupyter notebook using Lagrange's equation
function dz = constrained_ode(t,z,g) 

x = z(1);
y = z(2);
xd = z(3);
yd = z(4);

M = [y, -x; x, y];

dz = [xd; yd; M\[x*g; -xd^2 - yd^2] ];

end