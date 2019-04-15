function [x,t, fail_count] = rk2_adapt_method2(f, ti, h0, tf, x0, tol)

%% TODO
% - include minimum h check

%%
% Let's just start out time at ti
t = [ti];

%%
% populate the initial condition
x(:,1) = x0;

%% Define some parameters for the time marching
% In a more general code, these should be changable from the function call.
max_num_steps = 3000;
h = h0;
j = 1;
fail_count = 0;


%% Integrate
%
err0 = 1;
alpha = 1/2;
beta =  0.08;
for i = 1:max_num_steps
   
    
    [z,err] = rk2_step( f, t(j), x(:,j), h );
    
    if err <= tol
        j = j+1;
        x(:,j) = z;
        t(j) = t(j-1) + h;
    else
        fail_count = fail_count + 1;
    end
    
    % Hairer and Wanner (1996), Solving Ordinary Diff Eq Vol II 
    % Equation (2.43c)
    h = h*(tol/err)^alpha * (err0/tol)^beta;
    err0 = err;
    
    if t(j) >= tf
        break
    end
        
end
