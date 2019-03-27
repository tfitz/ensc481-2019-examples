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
nframes = 300;


%% Define the actuation
% No forcing:
f = @(t,z) 0;
torque = @(t,z) 0;

% feedback
K = design_controller(m, g, l, M, J_G);
% f = @(t,z) -K(1,:)*z;
% torque = @(t,z) -K(2,:)*z;


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
hold(ax, 'on');
set(ax, 'DataAspectRatio', [1,1,1]);
xlabel(ax,'x');
ylabel(ax, 'y');
h_title = title(ax, 'Title', 'FontSize', 12, 'FontWeight', 'Bold');

% make the cart
f_cart_x = @(z) z(1) + 0.5/2*[-1, +1, +1, -1];
h_cart = patch(ax, 'XData', f_cart_x(z0), 'YData', [-1,-1,1,1]*0.25/2 );
h_cart.FaceColor = [6, 39, 79]/255;
h_cart.FaceAlpha = 0.9;

% make the rod
[xpts,ypts] = rot_pendulum(z0, l);
h_rod = patch(ax, 'XData', xpts, 'YData', ypts);
h_rod.FaceColor = [0, 125, 138]/255;
h_rod.FaceAlpha = 0.5;

% make the center of mass of the potato
f_com_x = @(z) z(1) + l*sin(z(2));
f_com_y = @(z) l*cos(z(2));
h_com = plot(ax, f_com_x(z0), f_com_y(z0), ...
    'Marker', 'o', ...
    'MarkerFaceColor', [204, 0, 51]/255, ...
    'MarkerSize', 12);

%%
% find the bounding box of the simulation
xmax = max(sol.y(1,:));
xmin = min(sol.y(1,:));
ax.XLim = [-1.1,1.1];
ax.YLim = [-l,l]*2.1;

%%
% Now let's march through the results
time = linspace(0,tf,nframes);
for j = 1:nframes
    z = deval(sol, time(j), 1:2);
    
    % update the title
    h_title.String = sprintf('Time = %.2f', time(j) );
    
    % update positions
    h_cart.XData = f_cart_x(z);
    [h_rod.XData, h_rod.YData] = rot_pendulum(z,l);
    h_com.XData = f_com_x(z);
    h_com.YData = f_com_y(z);
    
    pause(0.05)
    
end


%% points of the pendulum
%
function [xpts,ypts] = rot_pendulum(z, l)

x = z(1);
theta = z(2);
R = [ sin(theta), -cos(theta);
    cos(theta), sin(theta)];
xpts = 2*l*[0, 1, 1, 0];
ypts = l/20*[-1, -1, 1, 1];

X = zeros(2,length(xpts));

for i = 1:length(xpts)
    X(:,i) = R*[ xpts(i); ypts(i) ];
end

xpts = X(1,:) + x;
ypts = X(2,:);

end

%% Controller design
function K = design_controller(m, g, l, M, J_G)
% rank is 4:
% CM = [0 0 (J_G + l.^2.*m)./(J_G.*M + J_G.*m + M.*l.^2.*m) -l.*m./(J_G.*M + J_G.*m + M.*l.^2.*m) 0 0 g.*l.^3.*m.^3./(J_G.*M + J_G.*m + M.*l.^2.*m).^2 -g.*l.^2.*m.^2.*(M + m)./(J_G.*M + J_G.*m + M.*l.^2.*m).^2; 0 0 -l.*m./(J_G.*M + J_G.*m + M.*l.^2.*m) (M + m)./(J_G.*M + J_G.*m + M.*l.^2.*m) 0 0 -g.*l.^2.*m.^2.*(M + m)./(J_G.*M + J_G.*m + M.*l.^2.*m).^2 g.*l.*m.*(M + m).^2./(J_G.*M + J_G.*m + M.*l.^2.*m).^2; 1 0 0 0 0 0 0 0; 0 1 0 0 -g.*l.^2.*m.^2./(J_G.*M + J_G.*m + M.*l.^2.*m) g.*l.*m.*(M + m)./(J_G.*M + J_G.*m + M.*l.^2.*m) 0 0];

A = [0 0 (J_G + l.^2.*m)./(J_G.*M + J_G.*m + M.*l.^2.*m) -l.*m./(J_G.*M + J_G.*m + M.*l.^2.*m); 0 0 -l.*m./(J_G.*M + J_G.*m + M.*l.^2.*m) (M + m)./(J_G.*M + J_G.*m + M.*l.^2.*m); 0 0 0 0; 0 2*l.^2.*m.^2.*(J_G.*M.*g.*l.*m + J_G.*g.*l.*m.^2 + M.*g.*l.^3.*m.^2)./(J_G.*M + J_G.*m + M.*l.^2.*m).^2 - (-J_G.*M.*g.*l.*m - J_G.*g.*l.*m.^2 - M.*g.*l.^3.*m.^2 + 2*g.*l.^3.*m.^3)./(J_G.*M + J_G.*m + M.*l.^2.*m) 0 0];
B = [0 0; 0 0; 1 0; 0 1];

p = [-3,-4,-2+2j, -2-2j];
K = place(A,B,p);

end



