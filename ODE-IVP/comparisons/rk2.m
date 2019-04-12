function [x,t] = rk2(f,ti,h,tf,x0)


t = ti:h:tf;

nstates = length(x0);
nsteps = length(t);

x = zeros(nstates, nsteps);

x(:,1) = x0;


for i = 2:nsteps
   
    x(:,i) = rk2_step(f, t(i-1), x(:,i-1),h);
    
end
