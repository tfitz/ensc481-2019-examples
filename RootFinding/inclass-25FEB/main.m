%% In class example
clear all
close all
clc

%% Let's make a plot
p = [0.7, 0.1, -2.5, 0, 1];


x = linspace(-2,2,100);

plot(x,polyval(p,x))
hold on
set(gca,'DataAspectRatio', [1,1,1])

%%
y = sqrt( 1 - x.^2/4);
plot( x, y, x, -y)

%% Build f(z)

f = @(z) [ z(1)^2/4 + z(2)^2 - 1;
          polyval(p,z(1)) - z(2)];
      
pd = polyder(p);      

J = @(z) [ z(1)/2, 2*z(2);
           polyval(pd,z(1)), -1];

       
%% Make an initial guess
z0 = [ 1.75; 0.5];
tol = 1e-8;
max_iter = 200;
iter = 0;
flag = 0;

fn = f(z0);
z = z0;
while flag == 0
   iter = iter + 1; 
   
   Jn = J(z);
   Dz = -Jn\fn;
   
   z = z + Dz;
   
   fn = f(z);
   
   err = norm(fn);
   
   if err <= tol
       flag = 1;
       
   elseif iter >= max_iter
       flag = -1;
   end 
    
end


plot( z(1), z(2), 'marker', 'o', 'MarkerSize', 12)



%%
z0 = [ 0.5; 0.5];
[z1,flag,iter] = myLittleNewton(z0,f,J);
plot( z1(1), z1(2), 'marker', '*', 'MarkerSize', 12)




