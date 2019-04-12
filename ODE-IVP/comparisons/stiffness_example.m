%% Stiff_example
% Inspired by article on the Mathworks
% [blog](https://www.mathworks.com/company/newsletters/articles/stiff-differential-equations.html)

%%
clear all
close all
clc

%%
delta = 0.01;
f = @(t,y) y^2*(1-y);
y0 = delta;
tf = 2/delta;

opts = odeset('RelTol',1.e-4);


%% ode45

ode45(f,[0 tf],delta,opts);

pause

%% ode23s
ode23s(f,[0 tf],delta,opts);

