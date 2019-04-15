function [x, t, fail_count] = rk2_adapt_method3(f, ti, h0, tf, x0)

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
for i = 1:max_num_steps
   
    
    [z,ztilde] = rk2_step2( f, t(j), x(:,j), h );
    
    err = compute_err( x(:,j), z, ztilde, 1e-3, 1e-3);
    
    if err <= 1
        j = j+1;
        x(:,j) = z;
        t(j) = t(j-1) + h;
    else
        fail_count = fail_count + 1;
    end
    
    % https://en.wikipedia.org/wiki/Adaptive_stepsize
    h = 0.9*h*min(max(1/err, 0.3),2);
    
    if t(j) >= tf
        break
    end
        
end
