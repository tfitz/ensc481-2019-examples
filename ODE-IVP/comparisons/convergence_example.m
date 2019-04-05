%% Convergence Example
% using AB2

%%
clear all
close all
clc

%% Define the range of step sizes
h_range = 2.^-(1:20);

%% Define sample problem
% The most interesting part to vary here is ratio of wn and stepsize h.
wn = 1;
zeta = 0.1;

%%
% The equation to solve is
% $$ \ddot{x} + 2 \zeta \omega_n \dot{x} + \omega_n^2 x = 0 $$
f = @(z) [z(2); -wn^2*z(1)-2*zeta*wn*z(2)];

%%
% With initial conditions
z0 = [1;0];

%%
% Let's also define the true solution
wd = wn*sqrt(1-zeta^2);
y = @(t) exp(-zeta*wn*t).*( z0(1)*cos(wd*t) + (z0(2) + zeta*wn*z0(1))/wd*sin(wd*t) );


%% Solve for each h

Zend = zeros(length(z0),length(h_range));

for j = 1:length(h_range)
    h = h_range(j);
    %%
    tf = 10; % be sure this is a multiple of all the h's
    time = 0:h:tf;
    nstep = length(time);
    
    %% Integrate the solution
    % Let's first take an Euler step since AB2 is not self-starting
    %Z(:,1) = z0;
    Z = z0;
    Z = Z + h*f(Z);
    
    %%
    % Now we can continue with AB2
    fm1 = f(Z);
    for i = 2:nstep-1
        
        fi = f(Z);
        Z = Z + h*( 3/2*fi - 1/2*fm1 );
        
        fm1 = fi;
        
    end
    
    %%
    % Save the Z(end)
    Zend(:,j) = Z;
    
end

%% Visualize the results
plot( h_range, Zend )

%% Compute the true error of x(t)
Ztrue = y(tf);
for j = length(h_range):-1:1
   err(j) = norm( Zend(1,j) - Ztrue ); 
end

figure
loglog( h_range, err, '-o', ...
    'LineWidth', 2, ...
    'DisplayName', 'Absolute Error' )
hold on
plot( h_range, h_range.^2, '--', ...
    'LineWidth', 2, ...
    'DisplayName', 'h^2')
xlabel('Step size h')
