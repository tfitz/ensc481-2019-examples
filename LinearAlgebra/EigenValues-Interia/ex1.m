%%
clear all
close all
clc

%%
[I, x, y, z, IG] = make_ellipsoid();

%% Visualize the body
fig = figure();
ax = axes('Parent', fig);
set(ax, 'DataAspectRatio', [1,1,1]);
hold(ax,'on');

surf(ax,x,y,z)
arrow3(zeros(3,3), eye(3)*4, '-1.5', 1.5 )

xlabel(ax,'x')
ylabel(ax,'y')
zlabel(ax,'z')
xlim(ax,[-4,4]);
ylim(ax,[-4,4]);
zlim(ax,[-4,4]);

%% Compute eigenvalues
[V,D] = eig(I);

%% Plot the eigenvectors
arrow3(zeros(3,3), V'*4, 'b-1.5', 1.5 )


