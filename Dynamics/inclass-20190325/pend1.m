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
theta0 = 30*pi/180;
dtheta0 = 0;
tf= 10;
z0 = [...
    l*cos(theta0);          % x(0)
    l*sin(theta0);          % y(0)
    -l*dtheta0*sin(theta0); % dx(0)
    l*dtheta0*cos(theta0)]; % dy(0)


%%
% Solve using ode45
tic
sol = ode45(@(t,z) constrained_ode(t,z,g), [0,tf], z0);
toc

%%
% solve the regular version
reg_ode = @(t,z) [z(2); -g/l*cos(z(1))];
tic
sol2 = ode45( reg_ode, [0,tf], [theta0; dtheta0]);
toc

%% Visualizations
% Let's compare the results,just remember they are in different
% coordinates.
time = linspace(0,tf,300);

figure
x = deval(sol, time, 1);
y = deval(sol, time, 2);
plot( time, x, 'LineWidth', 2, ...
    'DisplayName', 'CL: x')
hold on
plot( time, y, 'LineWidth', 2, ...
    'DisplayName', 'CL: y')

theta = deval(sol2, time,1);
plot( time, l*cos(theta), '--', ...
    'LineWidth', 2, ...
    'DisplayName', 'L: x(\theta)')
plot( time, l*sin(theta), '--', ...
    'LineWidth', 2,...
    'DisplayName', 'L: y(\theta)')

%%
% solve using an Implicit Runge-Kutta method w/ fixed step size
flag_irk = 1;
if flag_irk
    tic
    sol_irk = irk(@(t,z) constrained_ode(t,z,g),0,tf,z0,0.01);
    toc
    
    ti = sol_irk.x;
    xi = sol_irk.y(1,:);
    yi = sol_irk.y(2,:);
    
    plot( ti, xi, ':', ...
        'LineWidth', 2, ...
        'DisplayName', 'irk x')
    plot( ti, yi, ':', ...
        'LineWidth', 2, ...
        'DisplayName', 'irk y')
    
end

%%
% finish up the plot
legend('show')
xlabel('time t')
ylabel('position')

%%
% We notice that the positions are drifting.  As time goes forward, the
% length of the pendulum is changing *gasp*.  That's not super-the-best.
% Let's look at the plot of length of the pendulum in time
figure
plot( time, x.^2 + y.^2 - l^2, 'LineWidth', 2, 'DisplayName', 'ode45' )
hold on
xlabel('time t')
ylabel('Contrained equation h(x,y)')

if flag_irk
    plot( ti, xi.^2 + yi.^2 - l^2, ':', 'LineWidth', 2, 'DisplayName', 'irk' )
    legend('show');
end

%%
% What about energy?
figure
xd = deval(sol, time, 3);
yd = deval(sol, time, 4);
H = 1/2*m*(xd.^2 + yd.^2) + m*g*y;
plot( time, H/H(1), 'LineWidth', 2, 'DisplayName', 'CL, ode45' )
hold on

% compute energy in terms of theta
theta = deval(sol2, time, 1);
thetad = deval(sol2, time, 2);
Hr = 1/2*m*l^2*thetad.^2 + m*g*l*sin(theta);
plot( time, Hr/H(1), '--', 'LineWidth', 2, 'DisplayName', 'L, ode45')

legend('show');
xlabel('time t')
ylabel('Normalized Total energy H(t)/H(0)')

if flag_irk
   xdi = sol_irk.y(3,:);
   ydi = sol_irk.y(4,:);
   Hi = 1/2*m*(xdi.^2 + ydi.^2) + m*g*yi;
   plot( ti, Hi/H(1), ':', 'LineWidth', 2, 'DisplayName', 'irk')
   legend('show');
end

%%
% Ouch.  That doesn't look.  Drift is real.  This is a topic we'll come
% back to later on when we discuss implicit methods.


%% Animation
% Let's end on a happy note, and make an animation.
fig = figure();
ax = axes('Parent', fig);
hold(ax, 'on');
ax.DataAspectRatio = [1,1,1];
xlabel(ax, 'x');
ylabel(ax, 'y');
xlim(ax, [-2.5,2.5]);
ylim(ax, [-3,1]);

% plot the origin
h_base = patch(ax, 'XData', 0.2*[-1,1,1,-1], 'YData', 0.1*[-1,-1,1,1]  );

% draw a circle
beta = linspace(0,2*pi,100);
h_circle = plot(ax, l*cos(beta), l*sin(beta), ...
    ':', 'color', [1,1,1]*0.3, 'LineWidth', 3 );


% draw the pendulum
h_pendulum = line(ax, 'XData', [0, z0(1)] , 'YData', [0,z0(2)],...
    'LineWidth', 3, 'color', 'm');

nframes = 300;
time = linspace(0,tf,nframes);
for i = 1:nframes
    z = deval(sol, time(i), 1:2);
    h_pendulum.XData = [0,z(1)];
    h_pendulum.YData = [0,z(2)];
    
    pause(0.1)
    
end


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