%% Main 1: controller design

clear all
close all
clc

%% Define parameters
m1 = 1;
m2 = 1;
k1 = 20;
k2 = 10;
c1 = 0.4;
c2 = 0.2;

%% Build arrays
A = [0 0 1 0; 0 0 0 1; 
    (-k1 - k2)./m1 k2./m1 (-c1 - c2)./m1 c2./m1; 
    k2./m2 -k2./m2 c2./m2 -c2./m2];
B = [0;0;0;1/m2];
C = [0, 1, 0, 0];

%% Build C.C.F. Transformations
CMx = [B, A*B, A^2*B, A^3*B];

eigA = eig(A);
coeffA = poly(eigA)
Az = diag([1,1,1],1);
Az(end,:) = -coeffA(end:-1:2) % C.C.F of A
Bz = [0;0;0;1];

CMz = [Bz, Az*Bz, Az^2*Bz, Az^3*Bz];

P = CMx/CMz;

%% Desired poles
p1 = [-2-2j, -2+2j, -4-4j, -4+4j];
coeffD = poly(p1);

temp = flip( coeffD - coeffA );
Kz = temp(1:end-1)

Kx = Kz/P

%%
% Or we just could have used Place:
K = place( A,B, p1)

%% Observer
p2= 5*p1;
L = place( A', C', p2)'

