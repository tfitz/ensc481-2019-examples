function [I, x, y, z, IG] = make_ellipsoid()

%% Define Ellipsoid parameters
m = 1;
a = 1;
b = 1.5;
c = 3;

%% Build Interia Tensor in body frame
IG = 1/5*m*diag([b^2+c^2, a^2+c^2, a^2+b^2]);

%% Build rotation Matrix
% u: the axis that the body will rotate about
u = [1; -1; 1];
u = u/norm(u,2);
% theta: the angle of the rotation
theta = 30*pi/180;
% rotation matrix based on the Rodrigues formula
R = cos(theta)*eye(3) + sin(theta)*crossmat(u) + (1-cos(theta))*u*u';

%% Rotate the Inertia Tensor
I = R*IG*R';

%% Build the Ellipsoid
[X,Y,Z] = ellipsoid(0,0,0,a,b,c);

%%
% Rotate the Ellipsoid
% x = R*X
x = zeros(size(X));
y = zeros(size(Y));
z = zeros(size(Z));

n = numel(x);
for i = 1:n
    temp = R*[ X(i); Y(i); Z(i) ];
    x(i) = temp(1);
    y(i) = temp(2);
    z(i) = temp(3);
end


end

function vx = crossmat(v)
vx = [...
    0, -v(3), v(2);
    v(3), 0, -v(1);
    -v(2), v(1), 0 ...
    ];
end