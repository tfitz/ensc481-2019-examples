function [y,flag] = irk_step(f, tn, xn, h)



c1 = 1/2 - sqrt(3)/6;
c2 = 1/2 + sqrt(3)/6;
a11 = 1/4;
a12 = 1/4 - sqrt(3)/6;
a21 = 1/4 + sqrt(3)/6;
a22 = 1/4;

b1 = 1/2;
b2 = 1/2;

tol = 1e-4;
nstates = length(xn);
iter_max = 40*nstates;

%%
% Let's solve with fixed-point iteration
iter = 0;
k1 = f( tn + c1*h, xn );
k2 = f( tn + c2*h, xn );

flag = 0;

while flag == 0
    iter = iter + 1;
    
    k1s = f( tn + c1*h, xn + h*(a11*k1 + a12*k2 ) );
    k2s = f( tn + c2*h, xn + h*(a21*k1 + a22*k2 ) );
    
    err = norm( [k1s; k2s] - [k1; k2], 2);
    k1 = k1s;
    k2 = k2s;
    
    if err <= tol
        flag = 1;
    elseif iter >= iter_max
        flag = -1;
    end
    
end

%% build y
y = xn + h*( b1*k1 + b2*k2 );