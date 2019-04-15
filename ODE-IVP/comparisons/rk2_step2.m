function [z,y] = rk2_step2(f, tn, xn, h)

alpha = 2/3;

c2 = alpha;
a21 = alpha;
b1 = 1 - 1/2/alpha;
b2 = 1/2/alpha;

k1 = f( tn, xn );
k2 = f( tn + c2*h, xn + h*a21*k1 );

% RK2
z = xn + h*(b1*k1 + b2*k2);

% Euler: lower order est
y = xn + h*k1;

