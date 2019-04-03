%% Linear Multi-step methods
% 3-Apr: Adams-Bashforth 2
%


%%
clear all
close all
clc

%% Define sample problem

wn = 1;
zeta = 0.1;
f = @(z) [z(2); -wn^2*z(1)-2*zeta*wn*z(2)];

z0 = [1;0];

%%
h = 0.1;
tf = 10;
time = 0:h:tf;
nstep = length(time);

Z = zeros( 2, nstep);

%% Euler
Z(:,1) = z0;
Z(:,2) = Z(:,1) + h*f(z0);


%% Integrate AB2
fm1 = f(Z(:,1));
for i = 2:nstep-1
    
    fi = f(Z(:,i));
    Z(:,i+1) = Z(:,i) + h*( 3/2*fi - 1/2*fm1 );
    
    fm1 = fi;
end

plot(time, Z(1,:), 'o-', 'DisplayName', 'AB2')
xlabel('time')
ylabel('x(t)')
legend('show')


