function [x,t, fail_count] = rk2_adapt(f, ti, h0, tf, x0, tol)

%% TODO
% - include minimum h check
% - improve step size control
% - count the number of failures

%%
% Since we don't know how many time steps we're going to take, we need to
% just keep making the arrays larger as we go.  This works in Matlab, since
% Matlab is very lax on sizing.  If we were in another language we may have
% to be very careful how we do this (perhaps allocate large chunks of
% arrays and then fill them etc).  Changing array sizes is a bit sloppy,
% and can get very expensive as the size of the array increases.  If we
% needed to compute a large number of time steps (or had a large number of
% states) we likely wouldn't want to save each timestep.

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
   
    
    [z,err] = rk2_step( f, t(j), x(:,j), h );
    
    if err >= tol
        h = h/2;
        fail_count = fail_count + 1;
    else
        j = j+1;
        x(:,j) = z;
        t(j) = t(j-1) + h;
        h = h*1.1;   
    end
    
    if t(j) >= tf
        break
    end
        
end
