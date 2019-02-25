function [z,flag,iter] = myLittleNewton(z0, f, J)

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